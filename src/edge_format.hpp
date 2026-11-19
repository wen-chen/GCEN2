#pragma once
/**
 * Edge-list lines with GCEN v1 number formatting: correlation "%f",
 * p-value / FDR "%e" (v1 network_build: std::to_string / double_to_string).
 * Shared by the Python EdgeFormatter and build_network_file.
 */
#include <algorithm>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>

struct EdgeRow {
    int i, j;
    double r, p, q;
    double t = 0.0;  // TOM
};

// Optional columns after the p-value, as a bit mask (a bool `with_q` is 1).
enum EdgeExtra : int { kEdgeFdr = 1, kEdgeTom = 2 };

// get(k) returns the EdgeRow of edge k; `extra` (EdgeExtra bits) selects the
// optional columns fdr (q) and tom (t), each printed "%e". The n
// lines are cut into nt parts; par(nt, work) must call work(t) for every
// t < nt (in any order or concurrently) and return when all are done.
template <class Get, class Par>
std::string format_edge_lines_with(const std::vector<std::string> &names, std::size_t n,
                                   int extra, int nt, Get get, Par par) {
    nt = std::max(1, nt);
    std::vector<std::string> parts(nt);
    auto work = [&](int t) {
        const std::size_t b = n * t / nt, e = n * (t + 1) / nt;
        std::string &s = parts[t];
        s.reserve((e - b) * 64);
        char buf[96];
        for (std::size_t k = b; k < e; ++k) {
            const EdgeRow row = get(k);
            s += names[row.i];
            s += '\t';
            s += names[row.j];
            int len;
            switch (extra) {
                case kEdgeFdr:
                    len = std::snprintf(buf, sizeof buf, "\t%f\t%e\t%e\n", row.r, row.p, row.q);
                    break;
                case kEdgeTom:
                    len = std::snprintf(buf, sizeof buf, "\t%f\t%e\t%e\n", row.r, row.p, row.t);
                    break;
                case kEdgeFdr | kEdgeTom:
                    len = std::snprintf(buf, sizeof buf, "\t%f\t%e\t%e\t%e\n", row.r, row.p,
                                        row.q, row.t);
                    break;
                default:
                    len = std::snprintf(buf, sizeof buf, "\t%f\t%e\n", row.r, row.p);
            }
            s.append(buf, static_cast<std::size_t>(len));
        }
    };
    par(nt, work);
    std::size_t total = 0;
    for (auto &s : parts) total += s.size();
    std::string out;
    out.reserve(total);
    for (auto &s : parts) {
        out += s;
        std::string().swap(s);
    }
    return out;
}

// As above with fresh threads (one per 50k lines, at most n_threads).
template <class Get>
std::string format_edge_lines(const std::vector<std::string> &names, std::size_t n, int extra,
                              int n_threads, Get get) {
    const int nt = std::max(1, std::min<int>(n_threads, static_cast<int>(n / 50000) + 1));
    return format_edge_lines_with(names, n, extra, nt, get, [](int k, auto &work) {
        std::vector<std::thread> threads;
        for (int t = 1; t < k; ++t) threads.emplace_back(work, t);
        work(0);
        for (auto &th : threads) th.join();
    });
}
