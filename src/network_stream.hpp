#pragma once
/**
 * A co-expression network streamed in blocks of gene rows, with the
 * Benjamini-Hochberg FDR when asked: the C++ counterpart of
 * gcen2.network.build_to_file / build (same edges, same values, same order).
 *
 * Without FDR every block is correlated, filtered and handed to `emit` at
 * once. With FDR the BH value of a p-value needs its rank among all tested
 * pairs with p <= pval (or all pairs with all_edges), so a first pass keeps
 * those p-values (8 bytes each, sorted at the end) and, up to
 * `fdr_keep_limit`, the edges that pass the final |r| / p / sign filter
 * (24 bytes each). The BH values of the kept edges are then looked up and
 * the edges emitted without correlating again; only when more edges than
 * that pass does a second pass recompute the blocks.
 */
#include <cstdint>
#include <cstdio>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "correlation.hpp"

#include "tom.hpp"

struct NetworkOptions {
    double cor = 0.0;
    double pval = 1.0;
    bool positive_only = false;
    bool fdr = false;
    bool all_edges = false;
    int64_t max_pairs_per_block = 500000;
    size_t fdr_keep_limit = size_t(1) << 24;  // edges kept between FDR passes
    // TOM column (tom.hpp): every emitted edge gets its TOM and, unless
    // all_edges, edges with TOM < tom_min are dropped
    TomScorer *tom = nullptr;
    double tom_min = 0.0;
};

// Edges of the rows [begin, end) passing the kernel pre-filter, sorted by
// (i, j) when `sorted`, otherwise in no particular order.
using Correlate = std::function<std::vector<CorEdge>(int begin, int end, double min_abs_cor,
                                                     double max_p, bool positive_only,
                                                     bool sorted)>;
// One block of final edges, in order; q (null without FDR) gives their BH
// values, t (null without TOM) their TOM.
using Emit = std::function<void(const std::vector<CorEdge> &edges, const std::vector<double> *q,
                                const std::vector<double> *t)>;

void stream_network(int n_genes, const NetworkOptions &opt, int n_threads,
                    const Correlate &correlate, const Emit &emit);

/**
 * The correlation kernel of `method` ("pearson", "spearman",
 * "pearson_nonzero", "pearson_pairwise_complete") on a genes x samples
 * row-major matrix, prepared once (ranks, centred values) and run on `pool`.
 * `data` is only read while this is called. Pearson / Spearman use the
 * product kernel (correlate_rows_gemm) when use_gemm_kernel(kernel,
 * n_samples), else the per-pair kernel; the dropout-aware methods ignore
 * `kernel` (but still validate it).
 */
Correlate make_correlate(const std::string &method, const double *data, int n_genes,
                         int n_samples, WorkerPool &pool, const std::string &kernel = "auto");

/**
 * The TOM scorer of `spec` (nullptr when not enabled) for pearson/spearman
 * data (signed network: positive_only); other methods are rejected.
 */
std::unique_ptr<TomScorer> make_tom_scorer(const TomSpec &spec, const std::string &method,
                                           const double *data, int n_genes, int n_samples,
                                           bool positive_only, WorkerPool &pool);
TomInfo tom_info(const TomScorer &scorer);
// The rows TomScorer correlates: prepare_rows (or, method "bicor", WGCNA's
// bicor weighting, to be correlated without centring), NaN for genes with a
// non-finite value.
std::vector<std::vector<double>> tom_rows(const std::string &method, const double *data,
                                          int n_genes, int n_samples);

/** Edge-list file writer: header, GCEN v1 number formatting, node count. */
class EdgeFileWriter {
public:
    // Takes ownership of fh (already opened for writing).
    EdgeFileWriter(std::FILE *fh, std::string path, const std::vector<std::string> &genes,
                   bool with_q, WorkerPool &pool, bool with_t = false);
    ~EdgeFileWriter();
    void write(const std::vector<CorEdge> &edges, const std::vector<double> *q,
               const std::vector<double> *t = nullptr);
    void close();  // throws std::runtime_error when writing failed
    int64_t node_count() const;
    int64_t edge_count() const { return edges_; }

private:
    void put(const std::string &s);
    std::FILE *fh_;
    std::string path_;
    const std::vector<std::string> &genes_;
    bool with_q_;
    WorkerPool &pool_;
    bool with_t_;
    std::vector<uint8_t> has_edge_;
    int64_t edges_ = 0;
};
