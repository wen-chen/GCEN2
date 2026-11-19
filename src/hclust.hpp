#pragma once
/**
 * Average-linkage hierarchical clustering with the result of R's
 * fastcluster::hclust(d, method = "average") (the clustering WGCNA uses):
 * the nearest-neighbour-chain algorithm with fastcluster's search order,
 * tie rules and Lance-Williams update (s * a + t * b in double), merges
 * stable-sorted by height and labelled as R does (singletons -1..-n,
 * clusters 1..n-1 by merge step; in each row a singleton before a cluster,
 * the smaller label first). Memory: n (n - 1) / 2 doubles.
 */
#include <cstdint>
#include <vector>

struct Dendrogram {
    std::vector<int> merge;      // (n - 1) x 2, row-major, R convention
    std::vector<double> height;  // n - 1
};

// `diss`: n x n row-major dissimilarities (only the upper triangle is read).
template <class T>
Dendrogram hclust_average(const T *diss, int64_t n);

// The same on a writable n x n double matrix, without a copy: the upper
// triangle (copied to the lower one first) is the working storage and is
// restored from the lower one at the end.
Dendrogram hclust_average_inplace(double *diss, int64_t n);
