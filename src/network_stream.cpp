#include "network_stream.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>

#include "edge_format.hpp"
#include "parallel.hpp"
#include "tom.hpp"

namespace {

// gcen2.network._row_blocks
std::vector<std::pair<int, int>> row_blocks(int n, int64_t max_pairs) {
    std::vector<std::pair<int, int>> blocks;
    int start = 0;
    int64_t acc = 0;
    for (int i = 0; i < n; ++i) {
        const int64_t pairs = static_cast<int64_t>(n) - i - 1;
        if (acc && acc + pairs > max_pairs) {
            blocks.emplace_back(start, i);
            start = i;
            acc = 0;
        }
        acc += pairs;
    }
    if (start < n) blocks.emplace_back(start, n);
    return blocks;
}

// p-values (>= 0, +inf or NaN) as integers in numpy's sort order: the bits of
// a non-negative double order like the double; NaN (any sign or payload)
// goes last as numpy puts it, -0.0 with 0.0.
constexpr uint64_t kNanKey = 0x7ff8000000000000ULL;

inline uint64_t p_key(double p) {
    if (std::isnan(p)) return kNanKey;
    if (p == 0.0) return 0;
    uint64_t k;
    std::memcpy(&k, &p, sizeof k);
    return k;
}

inline double key_value(uint64_t k) {
    double p;
    std::memcpy(&p, &k, sizeof p);
    return p;
}

// gcen2.network._bh_lookup over the p-value keys: BH values of the sorted
// p-values with m tests in all, looked up for any of those p-values.
class BhLookup {
public:
    BhLookup(std::vector<uint64_t> keys, int64_t m, int n_threads) : keys_(std::move(keys)) {
        gcen_parallel::parallel_sort(keys_, std::less<uint64_t>(), n_threads);
        const size_t k = keys_.size();
        adjusted_.resize(k);
        const double total = static_cast<double>(m);
        for (size_t i = 0; i < k; ++i)  // sorted_p * n_total / rank, as numpy
            adjusted_[i] = key_value(keys_[i]) * total / static_cast<double>(i + 1);
        // np.minimum.accumulate from the right (NaN propagates), then capped at 1
        for (size_t i = k; i-- > 1;) {
            const double a = adjusted_[i], b = adjusted_[i - 1];
            if (std::isnan(a)) adjusted_[i - 1] = a;
            else if (!std::isnan(b) && a < b) adjusted_[i - 1] = a;
        }
        for (double &v : adjusted_)
            if (v > 1.0) v = 1.0;
    }

    double operator()(double p) const {
        const auto it = std::lower_bound(keys_.begin(), keys_.end(), p_key(p));
        const size_t pos = static_cast<size_t>(it - keys_.begin());
        return adjusted_[std::min(pos, adjusted_.size() - 1)];
    }

private:
    std::vector<uint64_t> keys_;
    std::vector<double> adjusted_;
};

struct Filters {
    double k_cor, k_p;
    bool k_pos;
    const NetworkOptions &o;

    explicit Filters(const NetworkOptions &opt) : o(opt) {
        // gcen2.network._kernel_filter
        if (o.all_edges) { k_cor = 0.0; k_p = std::numeric_limits<double>::infinity(); k_pos = false; }
        else if (o.fdr) { k_cor = 0.0; k_p = o.pval; k_pos = false; }
        else { k_cor = o.cor; k_p = o.pval; k_pos = o.positive_only; }
    }
    // gcen2.network._final_mask without the FDR term
    bool keep(const CorEdge &x) const {
        return o.all_edges || (std::fabs(x.corr) >= o.cor && x.pval <= o.pval &&
                               (!o.positive_only || x.corr >= 0.0));
    }
};

// BH values of `edges`; with all_edges every edge is kept, otherwise only
// those with q <= pval (edges and q compacted in place).
void apply_fdr(std::vector<CorEdge> &edges, std::vector<double> &q, const BhLookup &bh,
               const NetworkOptions &o, int n_threads) {
    const size_t n = edges.size();
    q.resize(n);
    const size_t T = gcen_parallel::workers_for(n, n_threads, size_t(1) << 14);
    gcen_parallel::run_parallel(T, [&](size_t t) {
        for (size_t k = n * t / T; k < n * (t + 1) / T; ++k) q[k] = bh(edges[k].pval);
    });
    if (o.all_edges) return;
    size_t w = 0;
    for (size_t k = 0; k < n; ++k) {
        if (q[k] <= o.pval) {
            edges[w] = edges[k];
            q[w++] = q[k];
        }
    }
    edges.resize(w);
    q.resize(w);
}

}  // namespace

void stream_network(int n_genes, const NetworkOptions &opt, int n_threads,
                    const Correlate &correlate, const Emit &emit_final) {
    const Filters f(opt);
    std::vector<double> t;
    // TOM of the final edges, then the tom_min filter (q compacted alongside)
    auto emit = [&](std::vector<CorEdge> &edges, std::vector<double> *q) {
        if (opt.tom) {
            opt.tom->score(edges, t);
            if (!opt.all_edges) {
                size_t w = 0;
                for (size_t k = 0; k < edges.size(); ++k) {
                    if (t[k] >= opt.tom_min) {
                        edges[w] = edges[k];
                        if (q) (*q)[w] = (*q)[k];
                        t[w++] = t[k];
                    }
                }
                edges.resize(w);
                t.resize(w);
                if (q) q->resize(w);
            }
        }
        if (!edges.empty()) emit_final(edges, q, opt.tom ? &t : nullptr);
    };
    const auto blocks = row_blocks(n_genes, std::max<int64_t>(1, opt.max_pairs_per_block));
    auto filtered = [&](int b, int e) {
        std::vector<CorEdge> edges = correlate(b, e, f.k_cor, f.k_p, f.k_pos, true);
        if (!opt.all_edges) {
            size_t w = 0;
            for (const CorEdge &x : edges)
                if (f.keep(x)) edges[w++] = x;
            edges.resize(w);
        }
        return edges;
    };

    if (!opt.fdr) {
        for (const auto &blk : blocks) {
            std::vector<CorEdge> edges = filtered(blk.first, blk.second);
            if (!edges.empty()) emit(edges, nullptr);
        }
        return;
    }

    // pass 1: the p-values of all tests BH ranks against (in any order), and
    // the edges that pass the other filters while they fit in fdr_keep_limit
    // (only those are sorted: most tests usually fail the |r| cutoff)
    std::vector<uint64_t> keys;
    std::vector<CorEdge> kept;
    std::vector<size_t> kept_end;  // end of each block in `kept`
    bool keeping = true;
    for (const auto &blk : blocks) {
        const std::vector<CorEdge> edges = correlate(blk.first, blk.second, f.k_cor, f.k_p,
                                                     f.k_pos, false);
        const size_t start = kept.size();
        for (const CorEdge &x : edges) {
            keys.push_back(p_key(x.pval));
            if (keeping && f.keep(x)) kept.push_back(x);
        }
        std::sort(kept.begin() + start, kept.end(), [](const CorEdge &x, const CorEdge &y) {
            return x.i != y.i ? x.i < y.i : x.j < y.j;
        });
        if (keeping && kept.size() > opt.fdr_keep_limit) {
            keeping = false;
            std::vector<CorEdge>().swap(kept);
            std::vector<size_t>().swap(kept_end);
        }
        if (keeping) kept_end.push_back(kept.size());
    }
    const int64_t n = n_genes;
    const BhLookup bh(std::move(keys), n * (n - 1) / 2, n_threads);

    std::vector<CorEdge> edges;
    std::vector<double> q;
    for (size_t b = 0; b < blocks.size(); ++b) {
        if (keeping) {
            const size_t lo = b ? kept_end[b - 1] : 0, hi = kept_end[b];
            edges.assign(kept.begin() + lo, kept.begin() + hi);
        } else {  // pass 2
            edges = filtered(blocks[b].first, blocks[b].second);
        }
        apply_fdr(edges, q, bh, opt, n_threads);
        if (!edges.empty()) emit(edges, &q);
    }
}

Correlate make_correlate(const std::string &method, const double *data, int n_genes,
                         int n_samples, WorkerPool &pool, const std::string &kernel) {
    const bool gemm = use_gemm_kernel(kernel, n_samples);
    if ((method == "pearson" || method == "spearman") && gemm) {
        auto rows = std::make_shared<GemmRows>(
            prepare_gemm_rows(prepare_rows(data, n_genes, n_samples, method == "spearman"),
                              n_samples));
        return [rows, &pool](int b, int e, double c, double p, bool pos, bool srt) {
            return correlate_rows_gemm(*rows, pool.size(), c, p, pos, b, e, &pool, srt);
        };
    }
    if (method == "pearson" || method == "spearman") {
        auto rows = std::make_shared<std::vector<std::vector<double>>>(
            prepare_rows(data, n_genes, n_samples, method == "spearman"));
        return [rows, n_samples, &pool](int b, int e, double c, double p, bool pos, bool srt) {
            return correlate_rows(*rows, n_samples, pool.size(), c, p, pos, b, e, &pool, srt);
        };
    }
    NonzeroMode mode;
    if (method == "pearson_nonzero") mode = NonzeroMode::BothNonzero;
    else if (method == "pearson_pairwise_complete") mode = NonzeroMode::AnyNonzero;
    else
        throw std::invalid_argument(
            "method must be 'pearson', 'spearman', 'pearson_nonzero' or "
            "'pearson_pairwise_complete'.");
    auto rows = std::make_shared<NonzeroRows>(prepare_nonzero_rows(data, n_genes, n_samples, mode));
    return [rows, &pool](int b, int e, double c, double p, bool pos, bool srt) {
        return correlate_rows_nonzero(*rows, pool.size(), c, p, pos, b, e, &pool, srt);
    };
}

// WGCNA bicor (maxPOutliers = 1, pearsonFallback = "individual"): x - median
// weighted by (1 - u^2)^2 for |u| < 1, u = (x - median) / (9 mad), mad the raw
// median absolute deviation; a gene with mad 0 falls back to x - mean.
// Correlating these rows without centring them again gives bicor.
static std::vector<double> bicor_row(const double *x, int n) {
    auto median = [](std::vector<double> v) {
        const size_t h = v.size() / 2;
        std::nth_element(v.begin(), v.begin() + h, v.end());
        const double hi = v[h];
        if (v.size() % 2) return hi;
        return (*std::max_element(v.begin(), v.begin() + h) + hi) / 2.0;
    };
    std::vector<double> v(x, x + n), out(n);
    const double med = median(v);
    for (int t = 0; t < n; ++t) v[t] = std::fabs(x[t] - med);
    const double mad = median(v);
    if (mad == 0) {
        double mean = 0;
        for (int t = 0; t < n; ++t) mean += x[t];
        mean /= n;
        for (int t = 0; t < n; ++t) out[t] = x[t] - mean;
        return out;
    }
    for (int t = 0; t < n; ++t) {
        const double u = (x[t] - med) / (9.0 * mad);
        const double w = std::fabs(u) < 1 ? (1 - u * u) * (1 - u * u) : 0.0;
        out[t] = (x[t] - med) * w;
    }
    return out;
}

std::vector<std::vector<double>> tom_rows(const std::string &method, const double *data,
                                          int n_genes, int n_samples) {
    std::vector<std::vector<double>> rows;
    if (method == "bicor") {
        rows.resize(n_genes);
        for (int i = 0; i < n_genes; ++i) {
            const double *x = data + static_cast<size_t>(i) * n_samples;
            if (std::all_of(x, x + n_samples, [](double v) { return std::isfinite(v); }))
                rows[i] = bicor_row(x, n_samples);
            else
                rows[i].assign(n_samples, 0.0);  // set to NaN below
        }
    } else {
        rows = prepare_rows(data, n_genes, n_samples, method == "spearman");
    }
    // a gene with a non-finite value gets adjacency 0 (its ranks would hide it)
    for (int i = 0; i < n_genes; ++i) {
        const double *x = data + static_cast<size_t>(i) * n_samples;
        if (!std::all_of(x, x + n_samples, [](double v) { return std::isfinite(v); }))
            std::fill(rows[i].begin(), rows[i].end(), std::numeric_limits<double>::quiet_NaN());
    }
    return rows;
}

std::unique_ptr<TomScorer> make_tom_scorer(const TomSpec &spec, const std::string &method,
                                           const double *data, int n_genes, int n_samples,
                                           bool positive_only, WorkerPool &pool) {
    if (!spec.enabled) return nullptr;
    if (method != "pearson" && method != "spearman")
        throw std::invalid_argument("TOM needs method 'pearson' or 'spearman'.");
    const auto rows = tom_rows(method, data, n_genes, n_samples);
    return std::unique_ptr<TomScorer>(new TomScorer(
        rows, n_samples, positive_only, spec.power,
        spec.powers.empty() ? default_soft_powers() : spec.powers, pool));
}

TomInfo tom_info(const TomScorer &scorer) {
    TomInfo info;
    info.power = scorer.power();
    info.estimated = scorer.estimated();
    info.fits = scorer.fits();
    return info;
}

EdgeFileWriter::EdgeFileWriter(std::FILE *fh, std::string path,
                               const std::vector<std::string> &genes, bool with_q,
                               WorkerPool &pool, bool with_t)
    : fh_(fh), path_(std::move(path)), genes_(genes), with_q_(with_q), pool_(pool),
      with_t_(with_t), has_edge_(genes.size(), 0) {
    put(std::string("gene_a\tgene_b\tcorrelation\tp_value") + (with_q ? "\tfdr" : "") +
        (with_t ? "\ttom" : "") + "\n");
}

EdgeFileWriter::~EdgeFileWriter() {
    if (fh_) std::fclose(fh_);
}

void EdgeFileWriter::put(const std::string &s) {
    if (!s.empty() && std::fwrite(s.data(), 1, s.size(), fh_) != s.size()) {
        std::fclose(fh_);
        fh_ = nullptr;
        throw std::runtime_error("error writing " + path_);
    }
}

void EdgeFileWriter::write(const std::vector<CorEdge> &edges, const std::vector<double> *q,
                           const std::vector<double> *t) {
    if (!fh_) throw std::runtime_error("error writing " + path_);
    // format and write in slices: with all edges kept, a whole block of text
    // would be ~60 bytes per pair (about 120 MB per 2M-pair block)
    constexpr size_t kSlice = 1 << 18;
    for (size_t off = 0; off < edges.size(); off += kSlice) {
        const size_t len = std::min(kSlice, edges.size() - off);
        // formatted on the pool's threads (no thread start-up per block)
        const int nt = len < 1024 ? 1 : pool_.size();
        put(format_edge_lines_with(
            genes_, len, (with_q_ ? kEdgeFdr : 0) | (with_t_ ? kEdgeTom : 0), nt,
            [&](size_t k) {
                const CorEdge &x = edges[off + k];
                return EdgeRow{x.i, x.j, x.corr, x.pval, q ? (*q)[off + k] : 0.0,
                               t ? (*t)[off + k] : 0.0};
            },
            [&](int k, auto &work) {
                if (k == 1) return work(0);
                std::atomic<int> next{0};
                pool_.run([&] {
                    for (int t; (t = next.fetch_add(1)) < k;) work(t);
                });
            }));
    }
    for (const CorEdge &x : edges) has_edge_[x.i] = has_edge_[x.j] = 1;
    edges_ += static_cast<int64_t>(edges.size());
}

void EdgeFileWriter::close() {
    if (!fh_) throw std::runtime_error("error writing " + path_);
    const int rc = std::fclose(fh_);
    fh_ = nullptr;
    if (rc != 0) throw std::runtime_error("error writing " + path_);
}

int64_t EdgeFileWriter::node_count() const {
    int64_t c = 0;
    for (const uint8_t h : has_edge_) c += h;
    return c;
}
