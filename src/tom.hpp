#pragma once
/**
 * Topological overlap (TOM) edge scores, as WGCNA computes them
 * (TOMsimilarityFromExpr, TOMDenom = "min"):
 *
 *   a_ij = ((1 + r_ij) / 2)^power   (signed network)   or   |r_ij|^power
 *   k_i  = sum_{u != i} a_iu,   l_ij = sum_{u != i,j} a_iu a_uj
 *   TOM_ij = (l_ij + a_ij) / (min(k_i, k_j) + 1 - a_ij)
 *
 * r is the Pearson correlation of the prepared rows (ranks for Spearman),
 * computed by a double-precision matrix product. The adjacency matrix is
 * kept in single precision (n^2 floats: 0.58 GB for 12,000 genes) and l_ij
 * comes from single-precision products of fixed tiles of it, computed only
 * for the tiles that hold requested edges. The products have a fixed
 * summation order (tom_gemm.hpp), so results do not depend on the number of
 * threads, on how the edges are batched or on the instruction set. Genes with a constant
 * or non-finite row get adjacency 0.
 *
 * With power <= 0 the power is picked as WGCNA 1.74's pickSoftThreshold
 * does: the first of `powers` whose scale-free fit R^2 (scaleFreeFitIndex,
 * 10 bins, plain R^2) exceeds 0.85; if none does, the WGCNA FAQ default for
 * the number of samples (signed 18/16/14/12, unsigned 9/8/7/6 for < 20,
 * < 30, < 40, >= 40 samples).
 */
#include <cstdint>
#include <memory>
#include <vector>

#include "correlation.hpp"

struct SoftThresholdFit {
    double power, r2, slope, mean_k, median_k, max_k;
};

// TOM column requested by a caller: power <= 0 picks it from `powers`.
struct TomSpec {
    bool enabled = false;
    double power = 0.0;
    std::vector<double> powers;
    double tom_min = 0.0;
};

// What the TOM step chose, for the caller to report.
struct TomInfo {
    double power = 0.0;
    bool estimated = false;
    std::vector<SoftThresholdFit> fits;
};

// WGCNA's default powerVector: 1..10, 12..20 by 2.
std::vector<double> default_soft_powers();

class TomScorer {
public:
    TomScorer(const std::vector<std::vector<double>> &rows, int n_samples, bool signed_net,
              double power, const std::vector<double> &powers, WorkerPool &pool);
    ~TomScorer();
    TomScorer(const TomScorer &) = delete;
    TomScorer &operator=(const TomScorer &) = delete;

    double power() const { return power_; }
    // power <= 0 given: whether a power reached R^2 > 0.85 (else the FAQ power)
    bool estimated() const { return estimated_; }
    const std::vector<SoftThresholdFit> &fits() const { return fits_; }  // auto only

    // TOM of each edge (i < j); edges come in (i, j) order, also across calls.
    void score(const std::vector<CorEdge> &edges, std::vector<double> &tom);

    int n_genes() const;

    struct Impl;  // defined in tom.cpp

private:
    std::unique_ptr<Impl> impl_;
    double power_ = 0.0;
    bool estimated_ = false;
    std::vector<SoftThresholdFit> fits_;
};

/**
 * 1 - TOM of all gene pairs in double precision (adjacency and products in
 * double: average-linkage trees are sensitive to errors of ~1e-8) into out
 * (n x n, row-major, diagonal 0); the power as TomScorer picks it, but from
 * double connectivities. Memory: `out` itself (the correlations and the
 * adjacency are built in it) plus about n^2 / 16 doubles.
 */
TomInfo tom_dissimilarity(const std::vector<std::vector<double>> &rows, int n_samples,
                          bool signed_net, double power, const std::vector<double> &powers,
                          WorkerPool &pool, double *out, bool center = true,
                          bool similarity = false);

/**
 * WGCNA's rowQuantileC on the columns of `sets` (n_sets arrays of length n):
 * out[i] = the q quantile of sets[0][i], ..., sets[n_sets - 1][i], computed by
 * WGCNA's pivot() (quickselect; the interpolation, and so its rounding,
 * follow the pivots, as in WGCNA's consensus TOM).
 */
void wgcna_row_quantile(const std::vector<const double *> &sets, int64_t n, double q,
                        double *out, WorkerPool &pool);

/** pickSoftThreshold's table without building the TOM. center = false:
 *  the rows are already centred (bicor_rows), as for tom_dissimilarity. */
std::vector<SoftThresholdFit> soft_threshold_fits(const std::vector<std::vector<double>> &rows,
                                                  int n_samples, bool signed_net,
                                                  const std::vector<double> &powers,
                                                  WorkerPool &pool, bool center = true);

/** WGCNA scaleFreeFitIndex(k, nBreaks) (Rsquared.SFT, slope.SFT) and the
 *  mean, median and max of k. */
SoftThresholdFit scale_free_fit(const std::vector<double> &k, double power, int n_breaks = 10);

/** The FAQ power when no power reaches the R^2 cutoff. */
double faq_power(int n_samples, bool signed_net);
