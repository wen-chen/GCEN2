#pragma once
/**
 * Streaming edge-list file operations behind `gcen2 merge` and `gcen2 diff`:
 * the file is read in chunks, so a large network never has to be held in
 * memory as text.
 */
#include <cstdint>
#include <string>
#include <vector>

struct EdgeFileSpec {
    std::string path;
    bool skip_header;
    int cor_column;  // 0-based column holding the correlation
};

struct MergedEdges {
    std::string text;       // "gene_a\tgene_b\tcorrelation\tdirection" rows
    int64_t node_count = 0;
    int64_t edge_count = 0;
};

// gcen2.network.merge on files: the correlations of each ordered pair
// (gene_a, gene_b) summed over all files (as pandas groupby().sum(), with
// its Kahan compensation and missing values skipped), pairs sorted by name,
// kept when |sum| >= cor_cutoff, the sum written as Python repr().
MergedEdges merge_edge_files(const std::vector<EdgeFileSpec> &files, double cor_cutoff);

// gcen2.differential.compare on files: the undirected edges of both files
// (the last row wins for a repeated edge, NaN for a missing or NA
// correlation), joined and sorted by (gene_a, gene_b) with the names in
// sorted order and gene_a <= gene_b.
struct DiffTable {
    std::vector<std::string> names;  // sorted
    std::vector<int32_t> lo, hi;     // per edge: indices into names, lo <= hi
    std::vector<double> cor_a, cor_b;
    std::vector<uint8_t> flags;      // bit 0: in A, bit 1: in B

    // Writes the result table (header + rows, as DataFrame.to_csv with
    // sep="\t") to `path`. `z`, `p` and `q` hold the Fisher z test and its
    // FDR for the shared edges, in row order. Returns the row counts of
    // A_only, B_only, shared_changed (p, or q with use_fdr, <= pval) and
    // shared_stable.
    // The rows are formatted by up to n_threads threads (the file does not
    // depend on n_threads).
    std::vector<int64_t> write(const std::string &path, const double *z, const double *p,
                               const double *q, size_t n_shared, double pval, bool use_fdr,
                               int n_threads = 1) const;
};

// Benjamini-Hochberg adjusted p-values of the finite entries of p (NaN for
// the others), as R's p.adjust(p, "BH") and gcen2.differential._bh, bit for
// bit: ranks by (p, index), p * m / rank, cumulative minimum from the
// largest p, capped at 1. Sorted by up to n_threads threads.
std::vector<double> bh_adjust(const double *p, size_t n, int n_threads);

// Both files are read in blocks of about block_bytes, parsed by up to
// n_threads threads; the result depends on neither.
DiffTable diff_join_edge_files(const EdgeFileSpec &a, const EdgeFileSpec &b, int n_threads = 1,
                               size_t block_bytes = size_t(1) << 26);
