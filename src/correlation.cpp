/**
 * correlation.cpp
 *
 * Pairwise Pearson / Spearman correlation computation for GCEN2.
 * Ported from GCEN v1 src/network_build.cpp + src/util/pearson.hpp.
 *
 * Threading model: two interleaved sections per thread (identical to GCEN v1)
 * so that the O(n²) pair space is partitioned evenly across threads.
 * Each thread accumulates its results locally, then merges under a mutex.
 */

#include "correlation.hpp"

#include "tom_gemm.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <atomic>
#include <mutex>
#include <thread>
#include <vector>

// ---------------------------------------------------------------------------
// Pearson correlation (identical to GCEN v1 pearson.hpp)
// ---------------------------------------------------------------------------
static double pearson(const double *x, const double *y, int n) {
    constexpr double TINY = std::numeric_limits<double>::min();
    double ax = 0.0, ay = 0.0;
    for (int i = 0; i < n; ++i) {
        ax += x[i];
        ay += y[i];
    }
    ax /= n;
    ay /= n;

    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (int i = 0; i < n; ++i) {
        double xt = x[i] - ax;
        double yt = y[i] - ay;
        sxx += xt * xt;
        syy += yt * yt;
        sxy += xt * yt;
    }
    return sxy / std::sqrt(sxx * syy + TINY);
}

// ---------------------------------------------------------------------------
// Regularized Incomplete Beta Function (from GCEN v1 pearson.hpp)
// Lewis Van Winkle, https://github.com/codeplea/incbeta
// ---------------------------------------------------------------------------
static double incbeta(double a, double b, double x) {
    const double TINY = 1.0e-30;
    const double STOP = 1.0e-8;

    if (x < 0.0 || x > 1.0) return std::numeric_limits<double>::infinity();
    if (x > (a + 1.0) / (a + b + 2.0))
        return 1.0 - incbeta(b, a, 1.0 - x);

    const double lbeta_ab = lgamma(a) + lgamma(b) - lgamma(a + b);
    const double front = std::exp(std::log(x) * a + std::log(1.0 - x) * b - lbeta_ab) / a;

    double f = 1.0, c = 1.0, d = 0.0;
    for (int i = 0; i <= 200; ++i) {
        int m = i / 2;
        double numerator;
        if (i == 0) {
            numerator = 1.0;
        } else if (i % 2 == 0) {
            numerator = (m * (b - m) * x) / ((a + 2.0 * m - 1.0) * (a + 2.0 * m));
        } else {
            numerator = -((a + m) * (a + b + m) * x) / ((a + 2.0 * m) * (a + 2.0 * m + 1));
        }
        d = 1.0 + numerator * d;
        if (std::fabs(d) < TINY) d = TINY;
        d = 1.0 / d;
        c = 1.0 + numerator / c;
        if (std::fabs(c) < TINY) c = TINY;
        const double cd = c * d;
        f *= cd;
        if (std::fabs(1.0 - cd) < STOP)
            return front * (f - 1.0);
    }
    return std::numeric_limits<double>::infinity();  // 1.0 / 0.0 in v1 (MSVC rejects it)
}

// Two-sided p-value from correlation r and sample size n (identical to GCEN v1)
static double get_p_value(double r, double n) {
    constexpr double TINY = std::numeric_limits<double>::min();
    double df = n - 2.0;
    double t = r * std::sqrt(df / ((1.0 - r + TINY) * (1.0 + r + TINY)));
    return incbeta(0.5 * df, 0.5, df / (df + t * t));
}

// A lower bound on |r| for p <= max_p with at most n samples (0: none). The
// p-value falls as |r| grows and rises as n shrinks, so a pair below the
// bound cannot pass the p filter: the kernels skip its p-value (an
// incomplete beta function, most of the cost of a pair) when the caller
// filters on p only (FDR) or with a smaller |r| cutoff. Bisection gives the
// largest |r| still above max_p; the margin covers rounding in get_p_value
// near that point, so the filtered pairs are exactly those of the p filter.
static double min_abs_cor_for_p(double max_p, int n) {
    if (!(max_p < 1.0) || !(max_p >= 0.0) || n < 3) return 0.0;
    double lo = 0.0, hi = 1.0;
    for (int it = 0; it < 64; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (get_p_value(mid, n) <= max_p) hi = mid;
        else lo = mid;
    }
    return std::max(0.0, lo * (1.0 - 1e-6) - 1e-9);
}


// ---------------------------------------------------------------------------
// Rank transformation for Spearman (identical to GCEN v1 base.hpp get_rank)
// Average rank for ties.
// ---------------------------------------------------------------------------
static std::vector<double> get_rank(const double *data, int n) {
    std::vector<int> index(n);
    std::vector<double> rank(n, 0.0);
    for (int i = 0; i < n; ++i) index[i] = i;

    std::sort(index.begin(), index.end(),
              [&](int a, int b) { return data[a] > data[b]; });

    int sumranks = 0, dupcount = 0;
    for (int i = 0; i < n; ++i) {
        sumranks += i;
        dupcount++;
        if (i == n - 1 || data[index[i]] != data[index[i + 1]]) {
            double averank = static_cast<double>(sumranks) / dupcount + 1.0;
            for (int j = i - dupcount + 1; j <= i; ++j)
                rank[index[j]] = averank;
            sumranks = 0;
            dupcount = 0;
        }
    }
    return rank;
}

// ---------------------------------------------------------------------------
// Thread worker — mirrors GCEN v1 thread_task_with_fdr logic
// (collects results locally, merges under mutex)
// ---------------------------------------------------------------------------

static std::mutex g_mutex;

struct CorFilter {
    double min_abs_cor;
    double max_p;
    bool positive_only;
    bool filter_cor;
    bool filter_p;
};

static CorFilter make_filter(double min_abs_cor, double max_p, bool positive_only,
                             int n_samples) {
    const double c = std::max(min_abs_cor, min_abs_cor_for_p(max_p, n_samples));
    return CorFilter{c, max_p, positive_only, c > 0.0, !std::isinf(max_p)};
}

// Compute one pair and append it to `out` if it passes the filter.
// |r| and sign are checked before the (expensive) p-value evaluation.
static inline void eval_pair(int i, int j, int sample_num,
                             const std::vector<std::vector<double>> &matrix,
                             const CorFilter &f,
                             std::vector<CorEdge> &out) {
    double corr = pearson(matrix[i].data(), matrix[j].data(), sample_num);
    if (f.filter_cor && !(std::fabs(corr) >= f.min_abs_cor)) return;
    if (f.positive_only && !(corr >= 0.0)) return;
    double pval = get_p_value(corr, sample_num);
    if (f.filter_p && !(pval <= f.max_p)) return;
    out.push_back({i, j, corr, pval});
}

// A unit of work: pairs (i, j) for j in [j_begin, j_end).
struct PairChunk {
    int i, j_begin, j_end;
};

// Pairs per work unit. Rows are split into column chunks so that a row range
// of only a few long rows (small blocks on a large matrix) still keeps every
// thread busy; units are handed out dynamically.
constexpr int kChunkPairs = 256;

template <class Eval>
static void thread_worker(std::atomic<size_t> &next, const std::vector<PairChunk> &units,
                          const Eval &eval, std::vector<CorEdge> &global_results) {
    std::vector<CorEdge> local;
    for (;;) {
        const size_t u = next.fetch_add(1);
        if (u >= units.size()) break;
        const PairChunk &c = units[u];
        for (int j = c.j_begin; j < c.j_end; ++j) eval(c.i, j, local);
    }
    std::lock_guard<std::mutex> lock(g_mutex);
    global_results.insert(global_results.end(), local.begin(), local.end());
}

// Runs job on the pool, or on min(n_threads, n_units) new threads.
static void run_job(WorkerPool *pool, int n_threads, size_t n_units,
                    const std::function<void()> &job) {
    if (pool) {
        pool->run(job);
        return;
    }
    const int n_run = static_cast<int>(std::min<size_t>(n_threads, std::max<size_t>(1, n_units)));
    std::vector<std::thread> threads;
    threads.reserve(n_run);
    for (int t = 0; t < n_run; ++t) threads.emplace_back(job);
    for (auto &th : threads) th.join();
}

// Deterministic output order (threads finish in arbitrary order).
static void sort_edges(std::vector<CorEdge> &edges) {
    std::sort(edges.begin(), edges.end(), [](const CorEdge &x, const CorEdge &y) {
        return x.i != y.i ? x.i < y.i : x.j < y.j;
    });
}

// Pairs (i, j > i) of rows [row_begin, row_end) evaluated by eval(i, j, out)
// on the threads (or the pool), sorted by (i, j).
template <class Eval>
static std::vector<CorEdge> run_pairs(int n_genes, int n_threads, int row_begin, int row_end,
                                      WorkerPool *pool, const Eval &eval, bool sorted = true) {
    if (n_threads < 1) n_threads = 1;
    if (row_end < 0 || row_end > n_genes) row_end = n_genes;
    if (row_begin < 0) row_begin = 0;
    std::vector<PairChunk> units;
    for (int i = row_begin; i < row_end; ++i)
        for (int jb = i + 1; jb < n_genes; jb += kChunkPairs)
            units.push_back({i, jb, std::min(n_genes, jb + kChunkPairs)});

    std::vector<CorEdge> global_results;
    std::atomic<size_t> next{0};
    run_job(pool, n_threads, units.size(),
            [&] { thread_worker(next, units, eval, global_results); });
    if (sorted) sort_edges(global_results);
    return global_results;
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------
std::vector<std::vector<double>> prepare_rows(const double *data, int n_genes, int n_samples,
                                              bool spearman) {
    // Build per-gene expression vectors (row-major input)
    std::vector<std::vector<double>> matrix(n_genes, std::vector<double>(n_samples));
    for (int i = 0; i < n_genes; ++i) {
        if (spearman) {
            // Pre-compute ranks for Spearman (same approach as GCEN v1)
            std::vector<double> ranks = get_rank(data + static_cast<size_t>(i) * n_samples,
                                                 n_samples);
            matrix[i] = std::move(ranks);
        } else {
            for (int j = 0; j < n_samples; ++j)
                matrix[i][j] = data[static_cast<size_t>(i) * n_samples + j];
        }
    }
    return matrix;
}

WorkerPool::WorkerPool(int n_threads) {
    for (int t = 1; t < std::max(1, n_threads); ++t) threads_.emplace_back(&WorkerPool::loop, this);
}

WorkerPool::~WorkerPool() {
    {
        std::lock_guard<std::mutex> lk(mu_);
        stop_ = true;
    }
    start_cv_.notify_all();
    for (auto &th : threads_) th.join();
}

void WorkerPool::loop() {
    std::uint64_t seen = 0;
    for (;;) {
        const std::function<void()> *job;
        {
            std::unique_lock<std::mutex> lk(mu_);
            start_cv_.wait(lk, [&] { return stop_ || generation_ != seen; });
            if (stop_) return;
            seen = generation_;
            job = job_;
        }
        (*job)();
        {
            std::lock_guard<std::mutex> lk(mu_);
            if (--pending_ == 0) done_cv_.notify_one();
        }
    }
}

void WorkerPool::run(const std::function<void()> &job) {
    {
        std::lock_guard<std::mutex> lk(mu_);
        job_ = &job;
        pending_ = static_cast<int>(threads_.size());
        ++generation_;
    }
    start_cv_.notify_all();
    job();
    std::unique_lock<std::mutex> lk(mu_);
    done_cv_.wait(lk, [&] { return pending_ == 0; });
}

std::vector<CorEdge> correlate_rows(const std::vector<std::vector<double>> &matrix,
                                    int n_samples, int n_threads, double min_abs_cor,
                                    double max_p, bool positive_only, int row_begin,
                                    int row_end, WorkerPool *pool, bool sorted) {
    const CorFilter filter = make_filter(min_abs_cor, max_p, positive_only, n_samples);
    return run_pairs(static_cast<int>(matrix.size()), n_threads, row_begin, row_end, pool,
                     [&](int i, int j, std::vector<CorEdge> &out) {
                         eval_pair(i, j, n_samples, matrix, filter, out);
                     }, sorted);
}

// ---------------------------------------------------------------------------
// Product-kernel path (see correlation.hpp).
// ---------------------------------------------------------------------------
bool use_gemm_kernel(const std::string &kernel, int n_samples) {
    if (kernel == "auto") return n_samples >= kGemmMinSamples;
    if (kernel == "gemm") return true;
    if (kernel == "pairwise") return false;
    throw std::invalid_argument("kernel must be 'auto', 'pairwise' or 'gemm'.");
}

GemmRows prepare_gemm_rows(const std::vector<std::vector<double>> &rows, int n_samples) {
    GemmRows g;
    g.n_genes = static_cast<int>(rows.size());
    g.n_samples = n_samples;
    g.z.assign(static_cast<size_t>(g.n_genes) * n_samples, 0.0);
    for (int i = 0; i < g.n_genes; ++i) {
        const double *x = rows[i].data();
        double a = 0.0;
        for (int k = 0; k < n_samples; ++k) a += x[k];
        a /= n_samples;
        double ss = 0.0;
        for (int k = 0; k < n_samples; ++k) ss += (x[k] - a) * (x[k] - a);
        double *z = g.z.data() + static_cast<size_t>(i) * n_samples;
        if (std::isnan(ss)) {  // a non-finite value: NaN correlations, as pearson()
            std::fill(z, z + n_samples, std::numeric_limits<double>::quiet_NaN());
        } else if (ss > 0.0) {
            const double s = 1.0 / std::sqrt(ss);
            for (int k = 0; k < n_samples; ++k) z[k] = (x[k] - a) * s;
        }  // constant row: zeros, r = 0 as pearson()
    }
    return g;
}

namespace {
constexpr int kGemmTileRows = 96;   // a multiple of the kernel's 6-row blocks
constexpr int kGemmTileCols = 256;
}  // namespace

std::vector<CorEdge> correlate_rows_gemm(const GemmRows &rows, int n_threads,
                                         double min_abs_cor, double max_p, bool positive_only,
                                         int row_begin, int row_end, WorkerPool *pool,
                                         bool sorted) {
    const int n = rows.n_genes, ns = rows.n_samples;
    const CorFilter f = make_filter(min_abs_cor, max_p, positive_only, ns);
    if (n_threads < 1) n_threads = 1;
    if (row_end < 0 || row_end > n) row_end = n;
    if (row_begin < 0) row_begin = 0;
    struct Tile { int r0, R, c0, C; };
    std::vector<Tile> units;
    for (int r0 = row_begin; r0 < row_end; r0 += kGemmTileRows) {
        const int R = std::min(kGemmTileRows, row_end - r0);
        for (int c0 = r0 + 1; c0 < n; c0 += kGemmTileCols)
            units.push_back({r0, R, c0, std::min(kGemmTileCols, n - c0)});
    }
    std::vector<CorEdge> global_results;
    std::atomic<size_t> next{0};
    const std::function<void()> job = [&] {
        std::vector<CorEdge> local;
        std::vector<double> buf;
        for (;;) {
            const size_t u = next.fetch_add(1);
            if (u >= units.size()) break;
            const Tile &t = units[u];
            buf.resize(static_cast<size_t>(t.R) * t.C);
            gemm_nt(rows.z.data(), static_cast<int64_t>(ns), t.r0, t.R, t.c0, t.C, ns,
                    buf.data());
            for (int a = 0; a < t.R; ++a) {
                const int i = t.r0 + a;
                const double *b = buf.data() + static_cast<size_t>(a) * t.C;
                for (int c = std::max(0, i + 1 - t.c0); c < t.C; ++c) {
                    double corr = b[c];
                    if (corr > 1.0) corr = 1.0;
                    else if (corr < -1.0) corr = -1.0;
                    if (f.filter_cor && !(std::fabs(corr) >= f.min_abs_cor)) continue;
                    if (f.positive_only && !(corr >= 0.0)) continue;
                    const double pval = get_p_value(corr, ns);
                    if (f.filter_p && !(pval <= f.max_p)) continue;
                    local.push_back({i, t.c0 + c, corr, pval});
                }
            }
        }
        std::lock_guard<std::mutex> lock(g_mutex);
        global_results.insert(global_results.end(), local.begin(), local.end());
    };
    run_job(pool, n_threads, units.size(), job);
    if (sorted) sort_edges(global_results);
    return global_results;
}

// ---------------------------------------------------------------------------
// Dropout-aware Pearson (see correlation.hpp).
//
// AnyNonzero: samples where at least one of the two genes is non-zero. Rows
// are centred once on their mean over all n samples (a); with n0 samples zero
// in both genes and S the others (n_S = n - n0), the centred moments on S
// follow exactly from the full-sample ones, since those n0 samples have
// x - a = -a:
//   sum_S (x - a)           = n0 * a
//   sum_S (x - a)^2         = ss_x - n0 * a^2
//   sum_S (x - a)(y - b)    = dot(x - a, y - b) - n0 * a * b
// and Sxy = sum_S (x - a)(y - b) - [sum_S (x - a)][sum_S (y - b)] / n_S (same
// for Sxx, Syy). With n0 = 0 this is exactly pearson() above.
//
// BothNonzero: samples where both genes are non-zero. Values are shifted by
// the gene's mean over its non-zero samples; the loop runs over the non-zero
// samples of the sparser gene, reading the other gene's shifted values (NaN
// where it is zero), and accumulates the shifted moments.
// ---------------------------------------------------------------------------
NonzeroRows prepare_nonzero_rows(const double *data, int n_genes, int n_samples,
                                 NonzeroMode mode) {
    NonzeroRows r;
    r.mode = mode;
    r.n_samples = n_samples;
    r.mean.resize(n_genes);
    r.values.assign(n_genes, std::vector<double>(n_samples));
    if (mode == NonzeroMode::AnyNonzero) {
        r.ss.resize(n_genes);
        r.zero.assign(n_genes, std::vector<uint64_t>((static_cast<size_t>(n_samples) + 63) / 64, 0));
    } else {
        r.nz.resize(n_genes);
    }
    for (int i = 0; i < n_genes; ++i) {
        const double *x = data + static_cast<size_t>(i) * n_samples;
        double *v = r.values[i].data();
        if (mode == NonzeroMode::AnyNonzero) {
            double a = 0.0;
            for (int k = 0; k < n_samples; ++k) a += x[k];
            a /= n_samples;
            double ss = 0.0;
            for (int k = 0; k < n_samples; ++k) {
                const double xt = x[k] - a;
                v[k] = xt;
                ss += xt * xt;
                if (x[k] == 0.0) r.zero[i][k / 64] |= uint64_t(1) << (k % 64);
            }
            r.mean[i] = a;
            r.ss[i] = ss;
        } else {
            double a = 0.0;
            int cnt = 0;
            for (int k = 0; k < n_samples; ++k)
                if (x[k] != 0.0) {
                    a += x[k];
                    ++cnt;
                    r.nz[i].push_back(k);
                }
            a = cnt ? a / cnt : 0.0;
            for (int k = 0; k < n_samples; ++k)
                v[k] = x[k] != 0.0 ? x[k] - a : std::numeric_limits<double>::quiet_NaN();
            r.mean[i] = a;
        }
    }
    return r;
}

static inline int popcount64(uint64_t v) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(v);
#else
    v = v - ((v >> 1) & 0x5555555555555555ULL);
    v = (v & 0x3333333333333333ULL) + ((v >> 2) & 0x3333333333333333ULL);
    v = (v + (v >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    return static_cast<int>((v * 0x0101010101010101ULL) >> 56);
#endif
}

// A gene is constant on the subset when its centred sum of squares (sxx) is
// zero up to the rounding of the shifted sum of squares (q) it came from.
static inline bool constant_on_subset(double sxx, double q) {
    return !(sxx > 1e-12 * q);
}

void nonzero_pearson(const NonzeroRows &rows, int i, int j, double &corr, double &n_eff) {
    constexpr double TINY = std::numeric_limits<double>::min();
    const int n = rows.n_samples;
    if (rows.mode == NonzeroMode::AnyNonzero) {
        const double *x = rows.values[i].data(), *y = rows.values[j].data();
        double dot = 0.0;
        for (int k = 0; k < n; ++k) dot += x[k] * y[k];
        const uint64_t *zx = rows.zero[i].data(), *zy = rows.zero[j].data();
        int n0 = 0;
        for (size_t w = 0; w < rows.zero[i].size(); ++w) n0 += popcount64(zx[w] & zy[w]);
        const int ns = n - n0;
        n_eff = ns;
        if (ns < 3) {
            corr = 0.0;
            return;
        }
        if (n0 == 0) {
            corr = dot / std::sqrt(rows.ss[i] * rows.ss[j] + TINY);
            return;
        }
        const double a = rows.mean[i], b = rows.mean[j];
        const double sx = n0 * a, sy = n0 * b;
        const double qx = rows.ss[i] - n0 * a * a, qy = rows.ss[j] - n0 * b * b;
        const double sxx = qx - sx * sx / ns;
        const double syy = qy - sy * sy / ns;
        const double sxy = dot - n0 * a * b - sx * sy / ns;
        if (constant_on_subset(sxx, qx) || constant_on_subset(syy, qy)) {
            corr = 0.0;
            return;
        }
        corr = std::max(-1.0, std::min(1.0, sxy / std::sqrt(sxx * syy + TINY)));
        return;
    }
    // BothNonzero: iterate over the sparser gene's non-zero samples
    const bool swap = rows.nz[i].size() > rows.nz[j].size();
    const int p = swap ? j : i, q = swap ? i : j;
    const double *xp = rows.values[p].data(), *xq = rows.values[q].data();
    int ns = 0;
    double sp = 0.0, sq = 0.0, spp = 0.0, sqq = 0.0, spq = 0.0;
    for (int k : rows.nz[p]) {
        const double vq = xq[k];
        if (vq != vq) continue;  // zero in q
        const double vp = xp[k];
        ++ns;
        sp += vp;
        sq += vq;
        spp += vp * vp;
        sqq += vq * vq;
        spq += vp * vq;
    }
    n_eff = ns;
    if (ns < 3) {
        corr = 0.0;
        return;
    }
    const double sxx = spp - sp * sp / ns, syy = sqq - sq * sq / ns;
    const double sxy = spq - sp * sq / ns;
    if (constant_on_subset(sxx, spp) || constant_on_subset(syy, sqq)) {
        corr = 0.0;
        return;
    }
    corr = std::max(-1.0, std::min(1.0, sxy / std::sqrt(sxx * syy + TINY)));
}

std::vector<CorEdge> correlate_rows_nonzero(const NonzeroRows &rows, int n_threads,
                                            double min_abs_cor, double max_p,
                                            bool positive_only, int row_begin, int row_end,
                                            WorkerPool *pool, bool sorted) {
    // n_eff <= n_samples: the bound for n_samples holds for every pair
    const CorFilter f = make_filter(min_abs_cor, max_p, positive_only, rows.n_samples);
    return run_pairs(static_cast<int>(rows.values.size()), n_threads, row_begin, row_end,
                     pool, [&](int i, int j, std::vector<CorEdge> &out) {
                         double corr, n_eff;
                         nonzero_pearson(rows, i, j, corr, n_eff);
                         if (f.filter_cor && !(std::fabs(corr) >= f.min_abs_cor)) return;
                         if (f.positive_only && !(corr >= 0.0)) return;
                         const double pval = n_eff < 3 ? 1.0 : get_p_value(corr, n_eff);
                         if (f.filter_p && !(pval <= f.max_p)) return;
                         out.push_back({i, j, corr, pval});
                     }, sorted);
}

std::vector<std::tuple<int, int, double, double>>
compute_pairwise_correlation(const double *data,
                             int n_genes,
                             int n_samples,
                             bool spearman,
                             int n_threads,
                             double min_abs_cor,
                             double max_p,
                             bool positive_only,
                             int row_begin,
                             int row_end) {
    const auto matrix = prepare_rows(data, n_genes, n_samples, spearman);
    std::vector<CorEdge> global_results = correlate_rows(
        matrix, n_samples, n_threads, min_abs_cor, max_p, positive_only, row_begin, row_end);

    // Convert to output format (free the intermediate buffer as we go)
    std::vector<std::tuple<int, int, double, double>> result;
    result.reserve(global_results.size());
    for (const auto &r : global_results)
        result.emplace_back(r.i, r.j, r.corr, r.pval);
    std::vector<CorEdge>().swap(global_results);

    return result;
}
