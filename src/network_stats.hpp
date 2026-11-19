#pragma once
/**
 * Per-node quantities behind gcen2.models.NetworkStats, computed on the
 * symmetric CSR adjacency of read_edge_adjacency / module._unique_adjacency
 * (rows sorted ascending, unique entries, a self-loop is a diagonal entry).
 */
#include <cstdint>
#include <vector>

struct GraphStatsArrays {
    std::vector<int64_t> degree_no_loop;  // neighbours other than the node itself
    std::vector<int64_t> self_loop;       // 1 if the node has a self-loop
    std::vector<int64_t> triangles;       // triangles through the node
    int64_t edge_count = 0;               // undirected edges, self-loops included
    int64_t components = 0;               // connected components
    int64_t giant_component = 0;          // nodes in the largest component
};

GraphStatsArrays graph_stats(const int64_t *indptr, const int32_t *indices, int32_t n,
                             int n_threads);
