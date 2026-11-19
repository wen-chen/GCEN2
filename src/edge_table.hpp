#pragma once
/**
 * Edge-list file to edge-list file, written as Network.to_edgelist would
 * write the table pandas reads from the input, without pandas:
 *
 *   shuffle_edge_file  gcen2.network.shuffle (numpy's random stream
 *                      reproduced in C++, pcg64.hpp)
 *   extract_edge_file  gcen2.network.extract (rows with a listed gene)
 *
 * Only inputs whose pandas reading is certain are handled; anything else
 * (bool-like or whitespace-padded values, numbers pandas might parse
 * differently, rows longer than the header, missing gene names) makes the
 * call return ok = false before anything is written, and the caller falls
 * back to the pandas pipeline. For extract the column dtypes still come
 * from every row (read_network reads the whole file).
 */
#include <cstdint>
#include <string>
#include <vector>

#include "pcg64.hpp"

struct EdgeTableSpec {
    std::string path;
    bool skip_header = false;
    int n_columns = 0;       // columns of the table (header or header-less rule)
    bool v1_layout = false;  // gene_a gene_b correlation p_value [fdr]: %f / %e
    std::string header;      // header line to write, with its newline
    std::string output;
    int n_threads = 4;
};

struct EdgeTableResult {
    bool ok = false;
    std::string reason;  // why the file was not handled (ok == false)
    int64_t node_count = 0;
    int64_t edge_count = 0;
    int64_t swaps_applied = 0;  // shuffle only
};

EdgeTableResult shuffle_edge_file(const EdgeTableSpec &spec, int64_t swap_factor, Pcg64 &rng);

// The rows whose first or second gene (unquoted) is one of `genes`.
EdgeTableResult extract_edge_file(const EdgeTableSpec &spec, const std::vector<std::string> &genes);
