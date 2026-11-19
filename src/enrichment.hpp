#pragma once
#include <cstddef>

/**
 * Compute a single hypergeometric p-value (upper tail, Fisher exact test).
 * Matches the algorithm in GCEN v1 hypergeometric_p_value.hpp.
 */
double hypergeometric_pvalue(int pop_n, int pop_count, int study_n, int study_count);

/**
 * Vectorized batch hypergeometric p-value computation.
 * The log-factorial table is precomputed once and reused across all tests.
 *
 * @param pop_n       Total population size.
 * @param pop_counts  Array of population counts per term.
 * @param study_n     Study set size.
 * @param study_counts Array of study counts per term.
 * @param out         Output array (preallocated, length n).
 * @param n           Number of terms.
 */
void batch_hypergeometric(int pop_n,
                          const int *pop_counts,
                          int study_n,
                          const int *study_counts,
                          double *out,
                          std::size_t n);

/**
 * One-sided hypergeometric upper tail P(X >= k), X ~ Hypergeom(M, n, N),
 * for arrays of integer-valued (k, M, n, N). Same algorithm as the numpy
 * gcen2.annotation._hypergeom_sf_numpy: each tail is summed away from the
 * mode with the pmf ratio recurrence (stops when terms fall below 1e-17 of
 * the running sum); k <= mode uses 1 - P(X <= k - 1).
 */
void hypergeom_sf(const double *k, const double *M, const double *n, const double *N,
                  double *out, std::size_t len, int n_threads);

#include <cstdint>
#include <vector>

/**
 * Study counts of many queries at once: for query r (genes
 * q_idx[q_ptr[r]:q_ptr[r+1]], unique) and term t, the number of its genes
 * annotated with t, where gene g has the unique terms
 * g_idx[g_ptr[g]:g_ptr[g+1]]. Same as the sparse product Q @ G.
 *
 * Only non-zero counts are returned, grouped by query (ascending r); the
 * order of terms within a query is unspecified.
 */
void study_counts(const int64_t *q_ptr, const int64_t *q_idx, int64_t n_queries,
                  const int64_t *g_ptr, const int32_t *g_idx, int32_t n_terms,
                  int n_threads, std::vector<int64_t> &rows,
                  std::vector<int64_t> &terms, std::vector<int64_t> &counts);

/**
 * study_counts fused with the one-sided test: for each query r and term t
 * with count c > 0 and pop_counts[t] > 0, p = hypergeom_sf(c, pop_n,
 * pop_counts[t], study_n[r]) (same arithmetic as hypergeom_sf); only entries
 * with p <= pval are returned (grouped by query). has_hits[r] is 1 when query
 * r has any term with c > 0 and pop_counts[t] > 0. Requires
 * 0 <= pop_counts[t], study_n[r] <= pop_n.
 *
 * fdr holds the Benjamini-Hochberg adjusted p of each returned entry, over
 * the m terms with pop_counts[t] > 0 (terms without a query gene have p = 1;
 * as R's p.adjust(method = "BH") on all m p-values of the query).
 */
void significant_counts(const int64_t *q_ptr, const int64_t *q_idx, int64_t n_queries,
                        const int64_t *g_ptr, const int32_t *g_idx, int32_t n_terms,
                        const int64_t *study_n, const int64_t *pop_counts, int64_t pop_n,
                        double pval, int n_threads, std::vector<int64_t> &rows,
                        std::vector<int64_t> &terms, std::vector<int64_t> &counts,
                        std::vector<double> &p, std::vector<double> &fdr,
                        std::vector<uint8_t> &has_hits);
