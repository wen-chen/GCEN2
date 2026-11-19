/**
 * module_identify.cpp — sparse Jaccard module identification for GCEN2.
 *
 * Same similarity definition and clustering as GCEN v1 module_identify.cpp
 * (verified identical on the v1 sample network), but:
 *   - input is a CSR adjacency matrix of integer node ids instead of string
 *     hash sets;
 *   - only 1- and 2-hop pairs are visited: for row i, shared-neighbour
 *     counts |N(i) ∩ N(j)| for all j > i are accumulated in a dense counter
 *     by walking i -> k -> j (this is row i of A·A);
 *   - rows are handed out to threads dynamically in small blocks;
 *   - components are found with union-find instead of recursive DFS.
 */

#include "module_identify.hpp"
#include "edge_text.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <atomic>
#include <functional>
#include <numeric>
#include <thread>
#include <utility>
#include <vector>

namespace {

struct UnionFind {
    std::vector<int32_t> parent;
    explicit UnionFind(int32_t n) : parent(n) {
        std::iota(parent.begin(), parent.end(), 0);
    }
    int32_t find(int32_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }
    void unite(int32_t a, int32_t b) {
        a = find(a);
        b = find(b);
        if (a == b) return;
        if (a < b) std::swap(a, b);
        parent[a] = b;  // keep the smaller id as root
    }
};

// Each worker unites its linked pairs in its own union-find (n ints) instead
// of collecting the pairs: modules are connected components, so merging the
// per-thread forests afterwards gives the same result without storing
// millions of pairs.
struct LinkSet {
    UnionFind uf;
    std::vector<uint8_t> linked;
    explicit LinkSet(int32_t n) : uf(n), linked(n, 0) {}
};

void similarity_worker(const int64_t *indptr, const int32_t *indices,
                       int32_t n, double cutoff, std::atomic<int32_t> &next_row,
                       LinkSet &out) {
    constexpr int32_t kBlock = 64;
    std::vector<int32_t> shared(n, 0);   // |N(i) ∩ N(j)| for the current row
    std::vector<uint8_t> adjacent(n, 0); // a(i, j) for the current row
    std::vector<int32_t> touched;
    touched.reserve(1024);

    for (;;) {
        const int32_t begin = next_row.fetch_add(kBlock);
        if (begin >= n) break;
        const int32_t end = std::min(n, begin + kBlock);

        for (int32_t i = begin; i < end; ++i) {
            const int64_t deg_i = indptr[i + 1] - indptr[i];
            for (int64_t p = indptr[i]; p < indptr[i + 1]; ++p) {
                const int32_t k = indices[p];
                if (k > i) {
                    if (!adjacent[k] && !shared[k]) touched.push_back(k);
                    adjacent[k] = 1;
                }
                for (int64_t q = indptr[k]; q < indptr[k + 1]; ++q) {
                    const int32_t j = indices[q];
                    if (j <= i) continue;
                    if (!adjacent[j] && !shared[j]) touched.push_back(j);
                    ++shared[j];
                }
            }
            for (const int32_t j : touched) {
                const double s = shared[j];
                const double a = adjacent[j];
                const double deg_j = static_cast<double>(indptr[j + 1] - indptr[j]);
                const double denom = static_cast<double>(deg_i) + deg_j - s - a;
                if (denom > 0.0 && (s + a) / denom > cutoff) {
                    out.uf.unite(i, j);
                    out.linked[i] = out.linked[j] = 1;
                }
                shared[j] = 0;
                adjacent[j] = 0;
            }
            touched.clear();
        }
    }
}

}  // namespace

std::vector<int32_t> jaccard_modules_csr(const int64_t *indptr,
                                         const int32_t *indices,
                                         int32_t n,
                                         double cutoff,
                                         int n_threads) {
    if (n <= 0) return {};
    n_threads = std::max(1, n_threads);

    std::atomic<int32_t> next_row{0};
    std::vector<LinkSet> links;
    links.reserve(n_threads);
    for (int t = 0; t < n_threads; ++t) links.emplace_back(n);
    if (n_threads == 1) {
        similarity_worker(indptr, indices, n, cutoff, next_row, links[0]);
    } else {
        std::vector<std::thread> pool;
        pool.reserve(n_threads);
        for (int t = 0; t < n_threads; ++t)
            pool.emplace_back(similarity_worker, indptr, indices, n, cutoff,
                              std::ref(next_row), std::ref(links[t]));
        for (auto &th : pool) th.join();
    }

    // merge the per-thread forests: same components as uniting every pair
    UnionFind uf(n);
    std::vector<uint8_t> linked(n, 0);
    for (auto &ls : links) {
        for (int32_t v = 0; v < n; ++v) {
            const int32_t r = ls.uf.find(v);
            if (r != v) uf.unite(v, r);
            linked[v] |= ls.linked[v];
        }
        ls = LinkSet(0);
    }

    std::vector<int32_t> label(n, -1);
    std::vector<int32_t> root_label(n, -1);
    int32_t next_label = 0;
    for (int32_t v = 0; v < n; ++v) {
        if (!linked[v]) continue;
        const int32_t r = uf.find(v);
        if (root_label[r] < 0) root_label[r] = next_label++;
        label[v] = root_label[r];
    }
    return label;
}

// ---------------------------------------------------------------------------
// Edge-list reader
// ---------------------------------------------------------------------------

namespace {

std::string unquote_field(const char *b, const char *e) {
    if (e - b >= 2 && *b == '"' && *(e - 1) == '"') {
        std::string out;
        out.reserve(static_cast<size_t>(e - b - 2));
        for (const char *p = b + 1; p < e - 1; ++p) {
            out.push_back(*p);
            if (*p == '"' && p + 1 < e - 1 && *(p + 1) == '"') ++p;
        }
        return out;
    }
    return std::string(b, e);
}

}  // namespace

EdgeAdjacency read_edge_adjacency(const std::string &path, bool skip_header,
                                  bool with_counts, bool descending) {
    std::unordered_map<std::string, int32_t> ids;
    std::vector<std::string> names;
    // Edge endpoints as node ids: uint16 while there are < 65536 nodes (the
    // usual case), int32 after that. The edge list is the largest buffer
    // next to the CSR it is turned into, so this lowers the peak memory.
    std::vector<uint16_t> ea16, eb16;
    std::vector<int32_t> ea, eb;
    bool wide = false;
    auto push_edge = [&](int32_t a, int32_t b) {
        if (!wide && names.size() > 65536) {
            ea.assign(ea16.begin(), ea16.end());
            eb.assign(eb16.begin(), eb16.end());
            ea16 = {};
            eb16 = {};
            wide = true;
        }
        if (wide) {
            ea.push_back(a);
            eb.push_back(b);
        } else {
            ea16.push_back(static_cast<uint16_t>(a));
            eb16.push_back(static_cast<uint16_t>(b));
        }
    };
    auto id_of = [&](std::string &&s) -> int32_t {
        auto it = ids.find(s);
        if (it != ids.end()) return it->second;
        const auto id = static_cast<int32_t>(names.size());
        names.push_back(s);
        ids.emplace(std::move(s), id);
        return id;
    };

    // lines read in large chunks (edge_text::for_each_line; '\r' removed)
    int64_t line_no = 0;
    bool header_pending = skip_header;
    edge_text::for_each_line(path, [&](const char *b, const char *e) {
        ++line_no;
        if (header_pending) { header_pending = false; return; }
        if (b == e) return;
        const char *t1 = static_cast<const char *>(std::memchr(b, '\t', e - b));
        if (!t1)
            throw std::runtime_error(path + ": line " + std::to_string(line_no) +
                                     " has fewer than two columns");
        const char *t2 = static_cast<const char *>(std::memchr(t1 + 1, '\t', e - t1 - 1));
        if (!t2) t2 = e;
        const int32_t a = id_of(unquote_field(b, t1));
        push_edge(a, id_of(unquote_field(t1 + 1, t2)));
    });
    ids = {};

    const auto n = static_cast<int32_t>(names.size());
    // rank nodes by name so the order matches the sorted pandas categories
    std::vector<int32_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int32_t x, int32_t y) { return names[x] < names[y]; });
    std::vector<int32_t> rank(n);
    for (int32_t r = 0; r < n; ++r) rank[order[r]] = r;

    EdgeAdjacency out;
    out.names.reserve(n);
    for (int32_t r = 0; r < n; ++r) out.names.push_back(std::move(names[order[r]]));
    names = {};
    order = {};

    // both directions; duplicates removed per row afterwards
    std::vector<int64_t> start(static_cast<size_t>(n) + 1, 0);
    std::vector<int32_t> cols;
    auto fill = [&](auto &xa, auto &xb) {
        for (size_t k = 0; k < xa.size(); ++k) {
            ++start[rank[xa[k]] + 1];
            ++start[rank[xb[k]] + 1];
        }
        for (int32_t i = 0; i < n; ++i) start[i + 1] += start[i];
        cols.resize(static_cast<size_t>(start[n]));
        std::vector<int64_t> pos(start.begin(), start.end() - 1);
        for (size_t k = 0; k < xa.size(); ++k) {
            const int32_t a = rank[xa[k]], b = rank[xb[k]];
            cols[pos[a]++] = b;
            cols[pos[b]++] = a;
        }
        std::decay_t<decltype(xa)>().swap(xa);
        std::decay_t<decltype(xb)>().swap(xb);
    };
    if (wide) fill(ea, eb); else fill(ea16, eb16);

    out.indptr.assign(static_cast<size_t>(n) + 1, 0);
    if (with_counts) out.counts.reserve(cols.size());
    int64_t w = 0;
    for (int32_t i = 0; i < n; ++i) {
        auto *b = cols.data() + start[i];
        auto *e = cols.data() + start[i + 1];
        if (descending) std::sort(b, e, std::greater<int32_t>());
        else std::sort(b, e);
        // compact runs of equal columns; w <= start[i], so in place is safe
        for (auto *p = b; p < e;) {
            auto *q = p + 1;
            while (q < e && *q == *p) ++q;
            if (with_counts) out.counts.push_back(static_cast<int32_t>(q - p));
            cols[w++] = *p;
            p = q;
        }
        out.indptr[i + 1] = w;
    }
    cols.resize(static_cast<size_t>(w));  // no shrink_to_fit: it would copy the array
    out.indices = std::move(cols);
    return out;
}
