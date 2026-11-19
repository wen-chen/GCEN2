#pragma once
#include <cstdint>
#include <vector>

/**
 * Identify gene modules using Jaccard similarity + connected components,
 * the GCEN v1 module_identify algorithm, on a sparse (CSR) adjacency matrix.
 *
 * For nodes i != j with neighbour sets N(i), N(j) (a self-loop puts a node in
 * its own set) and a(i,j) = 1 if i and j are adjacent:
 *
 *   J(i,j) = (|N(i) ∩ N(j)| + a(i,j)) / (|N(i)| + |N(j)| - |N(i) ∩ N(j)| - a(i,j))
 *
 * Pairs with J > cutoff are linked; each connected component of those links
 * is a module. Because cutoff > 0 only pairs that are adjacent or share a
 * neighbour can qualify, so only 1- and 2-hop pairs are visited (instead of
 * all n^2/2 pairs as in v1).
 *
 * @param indptr     CSR row pointer (n + 1 entries).
 * @param indices    CSR column indices; the matrix must be symmetric, with
 *                   unique, in-range entries per row.
 * @param n          Number of nodes.
 * @param cutoff     Jaccard similarity cutoff (> 0).
 * @param n_threads  Worker threads for the similarity phase.
 * @return Per-node component label: -1 for nodes without any similarity
 *         link, otherwise 0..k-1 numbered by first appearance in node order.
 */
std::vector<int32_t> jaccard_modules_csr(const int64_t *indptr,
                                         const int32_t *indices,
                                         int32_t n,
                                         double cutoff,
                                         int n_threads);

#include <string>

/** Symmetric 0/1 adjacency read straight from an edge-list file. */
struct EdgeAdjacency {
    std::vector<std::string> names;  // sorted (byte order), node i = names[i]
    std::vector<int64_t> indptr;     // CSR, n + 1 entries
    std::vector<int32_t> indices;    // sorted, unique per row
    std::vector<int32_t> counts;     // multiplicity of each entry (if requested)
};

/**
 * Read the first two tab-separated columns of an edge list (GCEN v1
 * ``.network`` / ``gcen2 build`` output) and build the symmetric CSR
 * adjacency: duplicate and reversed edges collapse to one entry, a
 * self-loop is a diagonal entry. Same matrix and node order as
 * ``gcen2.module._unique_adjacency`` on ``read_network(path, genes_only=True)``.
 *
 * Empty lines are skipped, a trailing '\r' is removed and a field wrapped
 * in double quotes is unquoted ("" -> "). Throws std::runtime_error if the
 * file cannot be read or a line has fewer than two columns.
 *
 * @param path         File path.
 * With ``with_counts`` the number of times each entry occurs is returned in
 * ``counts`` (an undirected edge listed twice counts 2, a self-loop counts 2
 * per line, as ``Network.sparse_adjacency``).
 *
 * @param path         File path.
 * @param skip_header  Skip the first line.
 * @param with_counts  Also return the entry multiplicities.
 * @param descending   Store each row's columns in descending order.
 */
EdgeAdjacency read_edge_adjacency(const std::string &path, bool skip_header,
                                  bool with_counts = false, bool descending = false);
