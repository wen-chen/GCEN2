#include "hclust.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

struct Step {
    int64_t node1, node2;
    double dist;
};

// fastcluster's doubly linked list of active nodes (indices ascending).
struct Active {
    std::vector<int64_t> succ, pred;
    int64_t start = 0;
    explicit Active(int64_t n) : succ(n + 1), pred(n + 1) {
        for (int64_t i = 0; i <= n; ++i) {
            pred[i] = i - 1;
            succ[i] = i + 1;
        }
    }
    void remove(int64_t idx) {
        if (idx == start) {
            start = succ[idx];
        } else {
            succ[pred[idx]] = succ[idx];
            pred[succ[idx]] = pred[idx];
        }
        succ[idx] = 0;  // inactive
    }
};

// fastcluster's union-find for relabelling (singletons 0..n-1, clusters n..)
struct UnionFind {
    std::vector<int64_t> parent;
    int64_t next;
    explicit UnionFind(int64_t n) : parent(2 * n - 1, 0), next(n) {}
    int64_t find(int64_t idx) {
        if (parent[idx] != 0) {
            int64_t p = idx;
            idx = parent[idx];
            if (parent[idx] != 0) {
                do { idx = parent[idx]; } while (parent[idx] != 0);
                do {  // path compression
                    const int64_t q = parent[p];
                    parent[p] = idx;
                    p = q;
                } while (parent[p] != idx);
            }
        }
        return idx;
    }
    void unite(int64_t a, int64_t b) { parent[a] = parent[b] = next++; }
};

}  // namespace

// at(i, j), i < j: the working dissimilarities (overwritten).
template <class At>
Dendrogram nn_chain_average(int64_t n, At at) {
    Dendrogram out;
    if (n < 2) return out;

    Active active(n);
    std::vector<double> members(n, 1.0);
    std::vector<int64_t> chain(n);
    int64_t tip = 0;
    std::vector<Step> steps;
    steps.reserve(n - 1);
    int64_t idx1, idx2;
    double mn;
    for (int64_t step = 0; step < n - 1; ++step) {
        if (tip <= 3) {
            chain[0] = idx1 = active.start;
            tip = 1;
            idx2 = active.succ[idx1];
            mn = at(idx1, idx2);
            for (int64_t i = active.succ[idx2]; i < n; i = active.succ[i])
                if (at(idx1, i) < mn) {
                    mn = at(idx1, i);
                    idx2 = i;
                }
        } else {
            tip -= 3;
            idx1 = chain[tip - 1];
            idx2 = chain[tip];
            mn = idx1 < idx2 ? at(idx1, idx2) : at(idx2, idx1);
        }
        do {
            chain[tip] = idx2;
            for (int64_t i = active.start; i < idx2; i = active.succ[i])
                if (at(i, idx2) < mn) {
                    mn = at(i, idx2);
                    idx1 = i;
                }
            for (int64_t i = active.succ[idx2]; i < n; i = active.succ[i])
                if (at(idx2, i) < mn) {
                    mn = at(idx2, i);
                    idx1 = i;
                }
            idx2 = idx1;
            idx1 = chain[tip++];
        } while (idx2 != chain[tip - 2]);

        steps.push_back({idx1, idx2, mn});
        if (idx1 > idx2) std::swap(idx1, idx2);
        const double size1 = members[idx1], size2 = members[idx2];
        const double s = size1 / (size1 + size2), t = size2 / (size1 + size2);
        active.remove(idx1);
        // average linkage: the merged cluster lives at idx2
        int64_t i = active.start;
        for (; i < idx1; i = active.succ[i]) at(i, idx2) = s * at(i, idx1) + t * at(i, idx2);
        for (; i < idx2; i = active.succ[i]) at(i, idx2) = s * at(idx1, i) + t * at(i, idx2);
        for (i = active.succ[idx2]; i < n; i = active.succ[i])
            at(idx2, i) = s * at(idx1, i) + t * at(idx2, i);
        members[idx2] = size1 + size2;
    }

    std::stable_sort(steps.begin(), steps.end(),
                     [](const Step &a, const Step &b) { return a.dist < b.dist; });
    UnionFind uf(n);
    out.merge.resize(static_cast<size_t>(2 * (n - 1)));
    out.height.resize(static_cast<size_t>(n - 1));
    for (int64_t k = 0; k < n - 1; ++k) {
        int64_t a = uf.find(steps[k].node1), b = uf.find(steps[k].node2);
        uf.unite(a, b);
        if (a > b) std::swap(a, b);
        out.merge[2 * k] = static_cast<int>(a < n ? -a - 1 : a - n + 1);
        out.merge[2 * k + 1] = static_cast<int>(b < n ? -b - 1 : b - n + 1);
        out.height[k] = steps[k].dist;
    }
    return out;
}

template <class T>
Dendrogram hclust_average(const T *diss, int64_t n) {
    if (n < 2) return Dendrogram();
    // condensed upper triangle: D_(i, j), i < j
    std::vector<double> D(static_cast<size_t>(n * (n - 1) / 2));
    auto at = [&](int64_t i, int64_t j) -> double & {
        return D[static_cast<size_t>((2 * n - 3 - i) * i / 2 + j - 1)];
    };
    for (int64_t i = 0; i < n; ++i)
        for (int64_t j = i + 1; j < n; ++j) {
            const double v = static_cast<double>(diss[i * n + j]);
            if (std::isnan(v)) throw std::invalid_argument("dissimilarities contain NaN.");
            at(i, j) = v;
        }
    return nn_chain_average(n, at);
}

Dendrogram hclust_average_inplace(double *diss, int64_t n) {
    for (int64_t i = 0; i < n * n; ++i)
        if (std::isnan(diss[i])) throw std::invalid_argument("dissimilarities contain NaN.");
    for (int64_t i = 0; i < n; ++i)  // the lower triangle keeps the input
        for (int64_t j = i + 1; j < n; ++j) diss[j * n + i] = diss[i * n + j];
    Dendrogram d = nn_chain_average(n, [&](int64_t i, int64_t j) -> double & {
        return diss[i * n + j];
    });
    for (int64_t i = 0; i < n; ++i)
        for (int64_t j = i + 1; j < n; ++j) diss[i * n + j] = diss[j * n + i];
    return d;
}

template Dendrogram hclust_average<float>(const float *, int64_t);
template Dendrogram hclust_average<double>(const double *, int64_t);
