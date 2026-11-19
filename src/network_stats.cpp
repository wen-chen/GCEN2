#include "network_stats.hpp"

#include <algorithm>
#include <atomic>
#include <thread>

GraphStatsArrays graph_stats(const int64_t *indptr, const int32_t *indices, int32_t n,
                             int n_threads) {
    GraphStatsArrays s;
    s.degree_no_loop.assign(n, 0);
    s.self_loop.assign(n, 0);
    s.triangles.assign(n, 0);
    int64_t loops = 0;
    for (int32_t i = 0; i < n; ++i) {
        const int32_t *b = indices + indptr[i], *e = indices + indptr[i + 1];
        const bool loop = std::binary_search(b, e, i);
        s.self_loop[i] = loop;
        s.degree_no_loop[i] = (e - b) - loop;
        loops += loop;
    }
    s.edge_count = (indptr[n] - loops) / 2 + loops;

    // Triangles: orient every edge from the lower to the higher (degree, id)
    // rank, so each triangle is found once from its lowest-ranked node and a
    // hub only scans its few higher-ranked neighbours. For node v the
    // out-neighbours are marked; for each out-neighbour u, the marked
    // entries of out(u) close a triangle (v, u, w).
    std::vector<int32_t> order(n);
    for (int32_t i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](int32_t a, int32_t b) {
        return s.degree_no_loop[a] != s.degree_no_loop[b]
            ? s.degree_no_loop[a] < s.degree_no_loop[b] : a < b;
    });
    std::vector<int32_t> rank(n);
    for (int32_t r = 0; r < n; ++r) rank[order[r]] = r;
    std::vector<int64_t> optr(static_cast<size_t>(n) + 1, 0);
    for (int32_t i = 0; i < n; ++i) {
        int64_t c = 0;
        for (int64_t k = indptr[i]; k < indptr[i + 1]; ++k) c += rank[indices[k]] > rank[i];
        optr[i + 1] = optr[i] + c;
    }
    std::vector<int32_t> out(static_cast<size_t>(optr[n]));
    for (int32_t i = 0; i < n; ++i) {
        int64_t w = optr[i];
        for (int64_t k = indptr[i]; k < indptr[i + 1]; ++k)
            if (rank[indices[k]] > rank[i]) out[w++] = indices[k];
    }

    n_threads = std::max(1, std::min(n_threads, static_cast<int>(std::max<int32_t>(n, 1))));
    std::vector<std::vector<int64_t>> part(n_threads, std::vector<int64_t>(n, 0));
    std::atomic<int32_t> next{0};
    auto work = [&](int t) {
        std::vector<int64_t> &tri = part[t];
        std::vector<uint8_t> mark(n, 0);
        for (;;) {
            const int32_t begin = next.fetch_add(64);
            if (begin >= n) break;
            const int32_t end = std::min(n, begin + 64);
            for (int32_t v = begin; v < end; ++v) {
                const int32_t *vb = out.data() + optr[v], *ve = out.data() + optr[v + 1];
                if (ve - vb < 2) continue;
                for (const int32_t *p = vb; p < ve; ++p) mark[*p] = 1;
                for (const int32_t *p = vb; p < ve; ++p) {
                    const int32_t u = *p;
                    int64_t found = 0;
                    for (int64_t k = optr[u]; k < optr[u + 1]; ++k) {
                        const int32_t w = out[k];
                        if (mark[w]) {
                            ++tri[w];
                            ++found;
                        }
                    }
                    tri[u] += found;
                    tri[v] += found;
                }
                for (const int32_t *p = vb; p < ve; ++p) mark[*p] = 0;
            }
        }
    };
    std::vector<std::thread> threads;
    for (int t = 1; t < n_threads; ++t) threads.emplace_back(work, t);
    work(0);
    for (auto &th : threads) th.join();
    for (const auto &p : part)
        for (int32_t i = 0; i < n; ++i) s.triangles[i] += p[i];

    // connected components (iterative DFS)
    std::vector<int32_t> label(n, -1), stack;
    for (int32_t r = 0; r < n; ++r) {
        if (label[r] >= 0) continue;
        int64_t size = 0;
        label[r] = static_cast<int32_t>(s.components);
        stack.push_back(r);
        while (!stack.empty()) {
            const int32_t v = stack.back();
            stack.pop_back();
            ++size;
            for (int64_t k = indptr[v]; k < indptr[v + 1]; ++k) {
                const int32_t w = indices[k];
                if (label[w] < 0) {
                    label[w] = label[r];
                    stack.push_back(w);
                }
            }
        }
        ++s.components;
        s.giant_component = std::max(s.giant_component, size);
    }
    return s;
}
