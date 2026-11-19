#pragma once
#include <cstdint>
#include <utility>
#include <vector>

/**
 * Random Walk with Restart (RWR) on a sparse transition matrix.
 *
 * Iteration:   PT+1 = (1-γ)·W·PT + γ·P0
 * Convergence: ||PT+1 - PT||₁ < tol
 *
 * Same update rule and stopping test as GCEN v1 src/util/rwr.hpp, but W is a
 * CSR matrix, so memory and time per step scale with the number of edges
 * instead of n² (a dense 27k-gene matrix needs ~6 GB).
 *
 * @param indptr, indices, data  CSR form of the column-normalised W (n × n).
 * @param P0        Seed probability vector (n), should sum to 1.
 * @param gamma     Restart probability.
 * @param max_iter  Maximum number of iterations.
 * @param tol       Convergence tolerance (L1 norm).
 * @return (steady-state vector PT, number of iterations performed).
 */
std::pair<std::vector<double>, int>
random_walk_restart_csr(const int64_t *indptr,
                        const int32_t *indices,
                        const double *data,
                        const double *P0,
                        int32_t n,
                        double gamma,
                        int max_iter,
                        double tol);

/**
 * The same walk on an unnormalised adjacency with entry multiplicities
 * ``counts``: W[i, j] = counts[p] / col_sum[j] (a zero column sum counts as
 * 1), formed on the fly as (1 / col_sum[j]) * counts[p], the value
 * ``(1.0 / col_sums)[indices] * counts`` has in numpy, so the result is
 * bit-identical to passing that array as ``data`` without storing it.
 */
std::pair<std::vector<double>, int>
random_walk_restart_counts(const int64_t *indptr,
                           const int32_t *indices,
                           const int32_t *counts,
                           const double *P0,
                           int32_t n,
                           double gamma,
                           int max_iter,
                           double tol);
