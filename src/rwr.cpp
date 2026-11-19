/**
 * rwr.cpp — sparse Random Walk with Restart for GCEN2.
 *
 * Ported from GCEN v1 src/util/rwr.hpp rwr(); W is now a CSR matrix.
 * W is expected to be pre-normalised (column sums = 1) by the caller.
 */

#include "rwr.hpp"

#include <cmath>

std::pair<std::vector<double>, int>
random_walk_restart_csr(const int64_t *indptr,
                        const int32_t *indices,
                        const double *data,
                        const double *P0,
                        int32_t n,
                        double gamma,
                        int max_iter,
                        double tol) {
    std::vector<double> PT(P0, P0 + n);
    std::vector<double> PT_new(n);
    const double walk = 1.0 - gamma;

    for (int step = 0; step < max_iter; ++step) {
        double delta = 0.0;
        for (int32_t i = 0; i < n; ++i) {
            double acc = 0.0;
            for (int64_t p = indptr[i]; p < indptr[i + 1]; ++p)
                acc += data[p] * PT[indices[p]];
            const double v = walk * acc + gamma * P0[i];
            delta += std::fabs(v - PT[i]);
            PT_new[i] = v;
        }
        PT.swap(PT_new);
        if (delta < tol) return {PT, step + 1};
    }
    return {PT, max_iter};
}

std::pair<std::vector<double>, int>
random_walk_restart_counts(const int64_t *indptr,
                           const int32_t *indices,
                           const int32_t *counts,
                           const double *P0,
                           int32_t n,
                           double gamma,
                           int max_iter,
                           double tol) {
    std::vector<double> inv(static_cast<size_t>(n), 0.0);
    const int64_t nnz = indptr[n];
    for (int64_t p = 0; p < nnz; ++p) inv[indices[p]] += static_cast<double>(counts[p]);
    for (auto &v : inv) v = 1.0 / (v == 0.0 ? 1.0 : v);

    std::vector<double> PT(P0, P0 + n);
    std::vector<double> PT_new(n);
    std::vector<double> q(static_cast<size_t>(n));
    const double walk = 1.0 - gamma;

    for (int step = 0; step < max_iter; ++step) {
        // q[j] = W[., j] * PT[j] for a count of 1 (inv[j] * 1.0 == inv[j]
        // exactly), so most entries need a single gather
        for (int32_t j = 0; j < n; ++j) q[j] = inv[j] * PT[j];
        double delta = 0.0;
        for (int32_t i = 0; i < n; ++i) {
            double acc = 0.0;
            for (int64_t p = indptr[i]; p < indptr[i + 1]; ++p) {
                const int32_t j = indices[p];
                const int32_t c = counts[p];
                acc += c == 1 ? q[j] : (inv[j] * static_cast<double>(c)) * PT[j];
            }
            const double v = walk * acc + gamma * P0[i];
            delta += std::fabs(v - PT[i]);
            PT_new[i] = v;
        }
        PT.swap(PT_new);
        if (delta < tol) return {PT, step + 1};
    }
    return {PT, max_iter};
}
