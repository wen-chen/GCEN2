#include "layer_fuse.hpp"

#include "parallel.hpp"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <exception>
#include <iterator>
#include <stdexcept>
#include <thread>
#include <unordered_map>

namespace {

using gcen_parallel::parallel_sort;
using gcen_parallel::run_parallel;
using gcen_parallel::workers_for;

struct Pair {
    int64_t code;
    double w;
};

// Unique pair codes (ascending) of one layer with their largest weight.
void dedup_layer(const FuseLayerCodes &L, int64_t n_genes, int n_threads,
                 std::vector<int64_t> &codes, std::vector<double> &w) {
    // bucket the rows by lo (counting sort), then sort each bucket by hi
    std::vector<int64_t> start(static_cast<size_t>(n_genes) + 1, 0);
    for (size_t r = 0; r < L.n; ++r) {
        const int64_t a = L.a[r], b = L.b[r];
        if (a < 0 || a >= n_genes || b < 0 || b >= n_genes)
            throw std::invalid_argument("gene code outside 0.." + std::to_string(n_genes));
        if (std::isnan(L.w[r])) throw std::domain_error("NaN weight");
        if (a != b) ++start[static_cast<size_t>(std::min(a, b)) + 1];
    }
    for (size_t g = 0; g < static_cast<size_t>(n_genes); ++g) start[g + 1] += start[g];
    const size_t m = static_cast<size_t>(start[static_cast<size_t>(n_genes)]);
    std::vector<Pair> rows(m);
    {
        std::vector<int64_t> fill(start.begin(), start.end() - 1);
        for (size_t r = 0; r < L.n; ++r) {
            const int64_t a = L.a[r], b = L.b[r];
            if (a == b) continue;
            const int64_t lo = std::min(a, b), hi = std::max(a, b);
            rows[static_cast<size_t>(fill[static_cast<size_t>(lo)]++)] =
                Pair{lo * n_genes + hi, L.w[r]};
        }
    }
    // sort the buckets and collapse repeats, in gene ranges per thread
    const size_t T = workers_for(m, n_threads);
    std::vector<size_t> gb(T + 1);  // gene boundaries with ~equal row counts
    gb[0] = 0;
    gb[T] = static_cast<size_t>(n_genes);
    for (size_t t = 1; t < T; ++t) {
        const int64_t target = static_cast<int64_t>(m * t / T);
        gb[t] = static_cast<size_t>(std::lower_bound(start.begin(), start.end(), target) -
                                    start.begin());
        gb[t] = std::min(std::max(gb[t], gb[t - 1]), static_cast<size_t>(n_genes));
    }
    std::vector<size_t> kept(T, 0);
    run_parallel(T, [&](size_t t) {
        size_t out = static_cast<size_t>(start[gb[t]]);
        for (size_t g = gb[t]; g < gb[t + 1]; ++g) {
            auto first = rows.begin() + start[g], last = rows.begin() + start[g + 1];
            std::sort(first, last, [](const Pair &x, const Pair &y) { return x.code < y.code; });
            for (auto it = first; it != last; ++it) {
                if (out > static_cast<size_t>(start[gb[t]]) && rows[out - 1].code == it->code) {
                    if (it->w > rows[out - 1].w) rows[out - 1].w = it->w;
                } else {
                    rows[out++] = *it;
                }
            }
        }
        kept[t] = out - static_cast<size_t>(start[gb[t]]);
    });
    size_t total = 0;
    for (size_t k : kept) total += k;
    codes.resize(total);
    w.resize(total);
    size_t pos = 0;
    for (size_t t = 0; t < T; ++t) {
        const size_t s = static_cast<size_t>(start[gb[t]]);
        for (size_t i = 0; i < kept[t]; ++i) {
            codes[pos] = rows[s + i].code;
            w[pos] = rows[s + i].w;
            ++pos;
        }
    }
}

void normalize_weights(std::vector<double> &w, FuseSpec::Normalize how, int n_threads) {
    if (w.empty() || how == FuseSpec::None) return;
    if (how == FuseSpec::Max) {
        const double mx = *std::max_element(w.begin(), w.end());
        if (mx > 0)
            for (double &x : w) x = x / mx;
        return;
    }
    // pandas rank(method="average", pct=True): tied values share the mean of
    // their 1-based ranks (an exact half-integer), divided by the count
    struct VI {
        double v;
        size_t i;
    };
    std::vector<VI> order(w.size());
    for (size_t i = 0; i < w.size(); ++i) order[i] = VI{w[i], i};
    parallel_sort(order, [](const VI &x, const VI &y) { return x.v < y.v; }, n_threads);
    const double n = static_cast<double>(w.size());
    for (size_t s = 0; s < order.size();) {
        size_t e = s + 1;
        while (e < order.size() && order[e].v == order[s].v) ++e;
        const double rank = (static_cast<double>(s + 1) + static_cast<double>(e)) / 2.0;
        for (size_t j = s; j < e; ++j) w[order[j].i] = rank / n;
        s = e;
    }
}

}  // namespace

FuseResult fuse_layer_codes(const std::vector<FuseLayerCodes> &layers, const FuseSpec &spec) {
    FuseResult r;
    const size_t L = layers.size();
    if (L == 0 || L > 62 || spec.alpha.size() != L) {
        r.reason = "need 1 to 62 layers, one alpha each";
        return r;
    }
    const int64_t n_genes = spec.n_genes;
    std::vector<std::vector<int64_t>> codes(L);
    std::vector<std::vector<double>> weights(L);
    try {
        for (size_t k = 0; k < L; ++k) {
            dedup_layer(layers[k], n_genes, spec.n_threads, codes[k], weights[k]);
            normalize_weights(weights[k], spec.normalize, spec.n_threads);
        }
    } catch (const std::domain_error &e) {
        r.reason = e.what();
        return r;
    }

    // union of the pair codes, merged layer by layer
    std::vector<int64_t> uni;
    for (size_t k = 0; k < L; ++k) {
        std::vector<int64_t> merged;
        merged.reserve(uni.size() + codes[k].size());
        std::set_union(uni.begin(), uni.end(), codes[k].begin(), codes[k].end(),
                       std::back_inserter(merged));
        uni.swap(merged);
    }
    const size_t U = uni.size();
    std::vector<std::vector<double>> wide(L, std::vector<double>(U, 0.0));
    std::vector<int64_t> bits(U, 0);
    run_parallel(L, [&](size_t k) {
        size_t p = 0;
        for (size_t i = 0; i < codes[k].size(); ++i) {
            while (uni[p] != codes[k][i]) ++p;
            wide[k][p] = weights[k][i];
        }
    });
    for (size_t k = 0; k < L; ++k) {
        size_t p = 0;
        for (size_t i = 0; i < codes[k].size(); ++i) {
            while (uni[p] != codes[k][i]) ++p;
            bits[p] |= int64_t(1) << k;
        }
        std::vector<int64_t>().swap(codes[k]);
        std::vector<double>().swap(weights[k]);
    }

    double total = 0.0;
    for (double x : spec.alpha) total += x;
    struct SI {
        double s;
        int64_t i;
    };
    std::vector<SI> rows;
    rows.reserve(U);
    for (size_t i = 0; i < U; ++i) {
        double s = 0.0;
        for (size_t k = 0; k < L; ++k) {
            const double term = wide[k][i] * spec.alpha[k];
            s = s + term;
        }
        s = s / total;
        const int nl = static_cast<int>(std::bitset<64>(static_cast<uint64_t>(bits[i])).count());
        if (nl >= spec.min_layers && s >= spec.min_score)
            rows.push_back(SI{s, static_cast<int64_t>(i)});
    }
    parallel_sort(rows,
                  [](const SI &x, const SI &y) { return x.s > y.s || (x.s == y.s && x.i < y.i); },
                  spec.n_threads);

    const size_t n = rows.size();
    r.a.resize(n);
    r.b.resize(n);
    r.score.resize(n);
    r.n_layers.resize(n);
    r.label_idx.resize(n);
    r.wide.assign(L, std::vector<double>(n));
    run_parallel(workers_for(n, spec.n_threads), [&, T = workers_for(n, spec.n_threads)](size_t t) {
        for (size_t j = n * t / T, e = n * (t + 1) / T; j < e; ++j) {
            const size_t i = static_cast<size_t>(rows[j].i);
            r.a[j] = uni[i] / n_genes;
            r.b[j] = uni[i] % n_genes;
            r.score[j] = rows[j].s;
            r.n_layers[j] =
                static_cast<int64_t>(std::bitset<64>(static_cast<uint64_t>(bits[i])).count());
            r.label_idx[j] = bits[i];  // the pattern, replaced by its index below
            for (size_t k = 0; k < L; ++k) r.wide[k][j] = wide[k][i];
        }
    });
    // presence patterns (few) in ascending order, and each row's index
    if (L <= 16) {
        std::vector<int64_t> index_of(size_t(1) << L, -1);
        for (int64_t x : r.label_idx) index_of[static_cast<size_t>(x)] = 0;
        for (size_t p = 0; p < index_of.size(); ++p)
            if (index_of[p] == 0) {
                index_of[p] = static_cast<int64_t>(r.patterns.size());
                r.patterns.push_back(static_cast<int64_t>(p));
            }
        for (auto &x : r.label_idx) x = index_of[static_cast<size_t>(x)];
    } else {
        std::unordered_map<int64_t, int64_t> index_of;
        for (int64_t x : r.label_idx) index_of.emplace(x, 0);
        for (auto &kv : index_of) r.patterns.push_back(kv.first);
        std::sort(r.patterns.begin(), r.patterns.end());
        for (size_t p = 0; p < r.patterns.size(); ++p)
            index_of[r.patterns[p]] = static_cast<int64_t>(p);
        for (auto &x : r.label_idx) x = index_of[x];
    }
    std::vector<char> seen(static_cast<size_t>(n_genes), 0);
    for (size_t j = 0; j < n; ++j) seen[static_cast<size_t>(r.a[j])] = seen[static_cast<size_t>(r.b[j])] = 1;
    for (char c : seen) r.node_count += c;
    r.ok = true;
    return r;
}
