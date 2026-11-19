/**
 * enrichment.cpp
 *
 * Hypergeometric p-value computation for GCEN2.
 * Ported from GCEN v1 src/util/hypergeometric_p_value.hpp.
 *
 * The log-factorial table is precomputed once per batch call and reused
 * across all tests to avoid redundant log() calls (same as GCEN v1 approach).
 */

#include "enrichment.hpp"

#include <cmath>
#include <vector>

// ---------------------------------------------------------------------------
// Log-factorial helpers (identical to GCEN v1 hypergeometric_p_value.hpp)
// ---------------------------------------------------------------------------
static void init_log_facs(std::vector<double> &log_facs, int n) {
    log_facs.resize(n + 1);
    log_facs[0] = 0.0;
    for (int i = 1; i <= n; ++i)
        log_facs[i] = log_facs[i - 1] + std::log(static_cast<double>(i));
}

static double log_hypergeometric_prob(const std::vector<double> &log_facs,
                                      int a, int b, int c, int d) {
    return log_facs[a + b] + log_facs[c + d] + log_facs[a + c] + log_facs[b + d]
         - log_facs[a] - log_facs[b] - log_facs[c] - log_facs[d]
         - log_facs[a + b + c + d];
}

// Fisher's exact test (two-sided: the tables no more likely than the
// observed one) with a log-factorial table covering a + b + c + d —
// identical to GCEN v1 fisher_exact()
static double fisher_exact(const std::vector<double> &log_facs, int a, int b, int c, int d) {
    const int n = a + b + c + d;
    const double log_p_cutoff = log_hypergeometric_prob(log_facs, a, b, c, d);
    double p_fraction = 0.0;
    for (int x = 0; x <= n; ++x) {
        if (a + b - x >= 0 && a + c - x >= 0 && d - a + x >= 0) {
            const double l = log_hypergeometric_prob(log_facs, x, a + b - x,
                                                     a + c - x, d - a + x);
            if (l <= log_p_cutoff)
                p_fraction += std::exp(l - log_p_cutoff);
        }
    }
    const double p = std::exp(log_p_cutoff + std::log(p_fraction));
    return (p > 1.0) ? 1.0 : p;
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------
double hypergeometric_pvalue(int pop_n, int pop_count,
                             int study_n, int study_count) {
    int a = study_count;
    int b = study_n   - study_count;
    int c = pop_count - study_count;
    int d = pop_n - pop_count - b;
    if (d < 0) return 1.0;   // degenerate case
    std::vector<double> log_facs;
    init_log_facs(log_facs, a + b + c + d);
    return fisher_exact(log_facs, a, b, c, d);
}

void batch_hypergeometric(int pop_n,
                          const int *pop_counts,
                          int study_n,
                          const int *study_counts,
                          double *out,
                          std::size_t n) {
    // Precompute log-factorial table once for the largest N we'll encounter.
    // Worst case table size = pop_n + 1.
    std::vector<double> log_facs;
    init_log_facs(log_facs, pop_n);

    for (std::size_t k = 0; k < n; ++k) {
        int a = study_counts[k];
        int b = study_n      - study_counts[k];
        int c = pop_counts[k] - study_counts[k];
        int d = pop_n - pop_counts[k] - b;
        if (d < 0 || a < 0 || b < 0 || c < 0) {
            out[k] = 1.0;
            continue;
        }

        out[k] = fisher_exact(log_facs, a, b, c, d);  // the shared table
    }
}

// ---------------------------------------------------------------------------
// One-sided upper tail (GCEN2 enrichment test)
// ---------------------------------------------------------------------------
#include <algorithm>
#include <functional>
#include <thread>

namespace {

struct LogFac {
    std::vector<double> t;  // t[i] = log(i!) via lgamma, like scipy gammaln
    explicit LogFac(double max_n) : t(static_cast<std::size_t>(max_n) + 1) {
        for (std::size_t i = 0; i < t.size(); ++i) t[i] = std::lgamma(static_cast<double>(i) + 1.0);
    }
    double binom(double a, double b) const {
        return t[static_cast<std::size_t>(a)] - t[static_cast<std::size_t>(b)]
             - t[static_cast<std::size_t>(a - b)];
    }
};

double sf_one(const LogFac &lf, double k, double M, double n, double N) {
    double xmin = std::max(0.0, N - (M - n));
    double xmax = std::min(n, N);
    if (k <= xmin) return 1.0;
    if (k > xmax) return 0.0;
    double mode = std::floor((N + 1) * (n + 1) / (M + 2));
    bool upward = k > mode;
    double x = upward ? k : k - 1;
    double term = std::exp(lf.binom(n, x) + lf.binom(M - n, N - x) - lf.binom(M, N));
    double total = term;
    if (upward) {
        while (x < xmax) {
            double ratio = (n - x) * (N - x) / ((x + 1) * (M - n - N + x + 1));
            x += 1;
            term *= ratio;
            total += term;
            if (!(term > total * 1e-17)) break;
        }
        return std::min(total, 1.0);
    }
    while (x > xmin) {
        double ratio = x * (M - n - N + x) / ((n - x + 1) * (N - x + 1));
        x -= 1;
        term *= ratio;
        total += term;
        if (!(term > total * 1e-17)) break;
    }
    return std::min(std::max(1.0 - total, 0.0), 1.0);
}

}  // namespace

void hypergeom_sf(const double *k, const double *M, const double *n, const double *N,
                  double *out, std::size_t len, int n_threads) {
    double max_m = 0.0;
    for (std::size_t i = 0; i < len; ++i) max_m = std::max(max_m, M[i]);
    LogFac lf(max_m);
    int nt = std::max(1, std::min<int>(n_threads, static_cast<int>(len / 20000) + 1));
    auto work = [&](int t) {
        std::size_t b = len * t / nt, e = len * (t + 1) / nt;
        for (std::size_t i = b; i < e; ++i) out[i] = sf_one(lf, k[i], M[i], n[i], N[i]);
    };
    std::vector<std::thread> th;
    for (int t = 1; t < nt; ++t) th.emplace_back(work, t);
    work(0);
    for (auto &x : th) x.join();
}

// ---------------------------------------------------------------------------
// Study counts (queries x genes) @ (genes x terms)
// ---------------------------------------------------------------------------

namespace {

struct CountPart {
    std::vector<int64_t> rows, terms, counts;
};

void study_counts_range(const int64_t *q_ptr, const int64_t *q_idx, int64_t r0, int64_t r1,
                        const int64_t *g_ptr, const int32_t *g_idx, int32_t n_terms,
                        CountPart &out) {
    std::vector<int64_t> counter(static_cast<size_t>(n_terms), 0);
    std::vector<int32_t> touched;
    for (int64_t r = r0; r < r1; ++r) {
        for (int64_t p = q_ptr[r]; p < q_ptr[r + 1]; ++p) {
            const int64_t g = q_idx[p];
            for (int64_t q = g_ptr[g]; q < g_ptr[g + 1]; ++q) {
                const int32_t t = g_idx[q];
                if (counter[t]++ == 0) touched.push_back(t);
            }
        }
        for (const int32_t t : touched) {
            out.rows.push_back(r);
            out.terms.push_back(t);
            out.counts.push_back(counter[t]);
            counter[t] = 0;
        }
        touched.clear();
    }
}

}  // namespace

void study_counts(const int64_t *q_ptr, const int64_t *q_idx, int64_t n_queries,
                  const int64_t *g_ptr, const int32_t *g_idx, int32_t n_terms,
                  int n_threads, std::vector<int64_t> &rows,
                  std::vector<int64_t> &terms, std::vector<int64_t> &counts) {
    n_threads = static_cast<int>(std::max<int64_t>(1, std::min<int64_t>(n_threads, n_queries)));
    std::vector<CountPart> parts(static_cast<size_t>(n_threads));
    if (n_threads == 1) {
        study_counts_range(q_ptr, q_idx, 0, n_queries, g_ptr, g_idx, n_terms, parts[0]);
    } else {
        std::vector<std::thread> pool;
        for (int t = 0; t < n_threads; ++t) {
            const int64_t r0 = n_queries * t / n_threads;
            const int64_t r1 = n_queries * (t + 1) / n_threads;
            pool.emplace_back(study_counts_range, q_ptr, q_idx, r0, r1, g_ptr, g_idx,
                              n_terms, std::ref(parts[t]));
        }
        for (auto &th : pool) th.join();
    }
    size_t total = 0;
    for (const auto &p : parts) total += p.rows.size();
    rows.clear(); terms.clear(); counts.clear();
    rows.reserve(total); terms.reserve(total); counts.reserve(total);
    for (auto &p : parts) {
        rows.insert(rows.end(), p.rows.begin(), p.rows.end());
        terms.insert(terms.end(), p.terms.begin(), p.terms.end());
        counts.insert(counts.end(), p.counts.begin(), p.counts.end());
        p = CountPart{};
    }
}

// ---------------------------------------------------------------------------
// Study counts fused with the test: only significant entries are kept
// ---------------------------------------------------------------------------

namespace {

struct SigPart {
    std::vector<int64_t> rows, terms, counts;
    std::vector<double> p, fdr;
};

// BH adjusted p of the entries with out_pos >= 0: rp holds every tested p of
// one query, m the number of tests (untested terms have p = 1 and do not
// change the adjusted p of the others, which p.adjust caps at 1).
void bh_row(std::vector<double> &rp, const std::vector<int64_t> &out_pos, double m,
            std::vector<size_t> &order, std::vector<double> &fdr) {
    const size_t k = rp.size();
    order.resize(k);
    for (size_t i = 0; i < k; ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&rp](size_t a, size_t b) { return rp[a] < rp[b] || (rp[a] == rp[b] && a < b); });
    double run = 1.0;
    for (size_t j = k; j-- > 0;) {
        const size_t i = order[j];
        run = std::min(run, rp[i] * m / static_cast<double>(j + 1));
        if (out_pos[i] >= 0) fdr[static_cast<size_t>(out_pos[i])] = run;
    }
}

void significant_range(const int64_t *q_ptr, const int64_t *q_idx, int64_t r0, int64_t r1,
                       const int64_t *g_ptr, const int32_t *g_idx, int32_t n_terms,
                       const int64_t *study_n, const int64_t *pop_counts, int64_t pop_n,
                       double pval, double m, const LogFac &lf, SigPart &out,
                       uint8_t *has_hits) {
    std::vector<int64_t> counter(static_cast<size_t>(n_terms), 0);
    std::vector<int32_t> touched;
    std::vector<double> rp;
    std::vector<int64_t> out_pos;
    std::vector<size_t> order;
    const double M = static_cast<double>(pop_n);
    for (int64_t r = r0; r < r1; ++r) {
        for (int64_t a = q_ptr[r]; a < q_ptr[r + 1]; ++a) {
            const int64_t g = q_idx[a];
            for (int64_t b = g_ptr[g]; b < g_ptr[g + 1]; ++b) {
                const int32_t t = g_idx[b];
                if (counter[t]++ == 0) touched.push_back(t);
            }
        }
        const double N = static_cast<double>(study_n[r]);
        uint8_t hit = 0;
        for (const int32_t t : touched) {
            const int64_t c = counter[t];
            counter[t] = 0;
            if (pop_counts[t] <= 0) continue;
            hit = 1;
            const double pv = sf_one(lf, static_cast<double>(c), M,
                                     static_cast<double>(pop_counts[t]), N);
            rp.push_back(pv);
            if (pv <= pval) {
                out_pos.push_back(static_cast<int64_t>(out.p.size()));
                out.rows.push_back(r);
                out.terms.push_back(t);
                out.counts.push_back(c);
                out.p.push_back(pv);
                out.fdr.push_back(1.0);
            } else {
                out_pos.push_back(-1);
            }
        }
        if (!rp.empty()) bh_row(rp, out_pos, m, order, out.fdr);
        rp.clear();
        out_pos.clear();
        has_hits[r] = hit;
        touched.clear();
    }
}

}  // namespace

void significant_counts(const int64_t *q_ptr, const int64_t *q_idx, int64_t n_queries,
                        const int64_t *g_ptr, const int32_t *g_idx, int32_t n_terms,
                        const int64_t *study_n, const int64_t *pop_counts, int64_t pop_n,
                        double pval, int n_threads, std::vector<int64_t> &rows,
                        std::vector<int64_t> &terms, std::vector<int64_t> &counts,
                        std::vector<double> &p, std::vector<double> &fdr,
                        std::vector<uint8_t> &has_hits) {
    LogFac lf(static_cast<double>(pop_n));
    double m = 0.0;
    for (int32_t t = 0; t < n_terms; ++t) m += pop_counts[t] > 0 ? 1.0 : 0.0;
    has_hits.assign(static_cast<size_t>(n_queries), 0);
    n_threads = static_cast<int>(std::max<int64_t>(1, std::min<int64_t>(n_threads, n_queries)));
    std::vector<SigPart> parts(static_cast<size_t>(n_threads));
    std::vector<std::thread> pool;
    for (int t = 0; t < n_threads; ++t) {
        const int64_t r0 = n_queries * t / n_threads;
        const int64_t r1 = n_queries * (t + 1) / n_threads;
        auto job = [&, r0, r1, t] {
            significant_range(q_ptr, q_idx, r0, r1, g_ptr, g_idx, n_terms, study_n,
                              pop_counts, pop_n, pval, m, lf, parts[t], has_hits.data());
        };
        if (t + 1 == n_threads) job(); else pool.emplace_back(job);
    }
    for (auto &th : pool) th.join();
    size_t total = 0;
    for (const auto &q : parts) total += q.rows.size();
    rows.clear(); terms.clear(); counts.clear(); p.clear(); fdr.clear();
    rows.reserve(total); terms.reserve(total); counts.reserve(total); p.reserve(total);
    fdr.reserve(total);
    for (auto &q : parts) {
        rows.insert(rows.end(), q.rows.begin(), q.rows.end());
        terms.insert(terms.end(), q.terms.begin(), q.terms.end());
        counts.insert(counts.end(), q.counts.begin(), q.counts.end());
        p.insert(p.end(), q.p.begin(), q.p.end());
        fdr.insert(fdr.end(), q.fdr.begin(), q.fdr.end());
        q = SigPart{};
    }
}
