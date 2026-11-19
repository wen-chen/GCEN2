#pragma once
#include <cstdint>
#include <string>

#include "tom.hpp"

/**
 * `gcen2 build` without numpy: read a plain expression matrix, compute the
 * correlations block by block and write the edge list, with the same result
 * as gcen2._matrix_io.read_matrix + gcen2.network.build_to_file (log
 * transforms are not handled here; FDR as in network_stream.hpp).
 *
 * Only a plain subset of matrix files is accepted: UTF-8, tab-separated,
 * values that are plain decimal numbers ([-]digits[.digits][e[+-]digits],
 * no quotes, spaces, NA strings or inf/nan), at least 3 columns and one gene.
 * Short rows are padded with NaN, as read_matrix does. For anything else, and
 * when the output file cannot be opened, `supported` is false and nothing is
 * written: the caller then uses the Python path, which handles (or reports)
 * those cases.
 */
struct BuildFileResult {
    bool supported = false;
    int64_t node_count = 0;
    int64_t edge_count = 0;
    TomInfo tom;
};

BuildFileResult build_network_file(const std::string &in_path, const std::string &out_path,
                                   bool has_header, bool spearman, double cor, double pval,
                                   bool positive_only, bool all_edges, int n_threads,
                                   int64_t max_pairs_per_block, bool fdr = false,
                                   const TomSpec &tom = TomSpec(),
                                   const std::string &kernel = "auto");
