#pragma once
/**
 * Weighted fusion of network layers given as gene codes, for
 * gcen2.multilayer.fuse / fuse_to_file: the same table as the numpy
 * reference in gcen2.multilayer._fuse_arrays_numpy, value for value.
 *
 * Per layer: undirected pairs (lo, hi) as one int64 lo * n_genes + hi,
 * self-loops dropped, repeated pairs keeping their largest weight, then the
 * weights scaled (max: w / max w when max w > 0; rank: pandas
 * ``rank(method="average", pct=True)``; none). Over the union of all pairs:
 * score = (sum_k alpha_k * w_k) / sum_k alpha_k (summed in layer order, an
 * absent pair counting as 0), rows kept when present in >= min_layers layers
 * and score >= min_score, ordered by score (descending), then pair code.
 *
 * Declines (ok = false) on NaN weights and on more than 62 layers; the
 * caller then uses the numpy reference.
 */
#include <cstdint>
#include <string>
#include <vector>

struct FuseLayerCodes {
    const int64_t *a = nullptr;
    const int64_t *b = nullptr;
    const double *w = nullptr;
    size_t n = 0;
};

struct FuseSpec {
    int64_t n_genes = 0;
    enum Normalize { Max, Rank, None } normalize = Max;
    std::vector<double> alpha;  // one per layer, >= 0 with a positive sum
    int min_layers = 1;
    double min_score = 0.0;
    int n_threads = 4;
};

struct FuseResult {
    bool ok = false;
    std::string reason;
    std::vector<int64_t> a, b;          // gene codes, a < b
    std::vector<double> score;
    std::vector<int64_t> n_layers;
    std::vector<int64_t> label_idx;     // index into patterns
    std::vector<int64_t> patterns;      // presence bit patterns, ascending
    std::vector<std::vector<double>> wide;  // per layer, one value per row
    int64_t node_count = 0;
};

// Throws std::invalid_argument for a gene code outside [0, n_genes).
FuseResult fuse_layer_codes(const std::vector<FuseLayerCodes> &layers, const FuseSpec &spec);
