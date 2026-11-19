#pragma once
// Small fork-join helpers shared by the kernels: parallel jobs, sorting,
// and writing a table formatted by several threads.
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <string>
#include <thread>
#include <vector>

namespace gcen_parallel {

// Runs job(t) for t in [0, n) on up to n threads (the caller's included).
template <class Job>
void run_parallel(size_t n, Job job) {
    std::vector<std::exception_ptr> errors(n);
    auto guarded = [&](size_t t) {
        try {
            job(t);
        } catch (...) {
            errors[t] = std::current_exception();
        }
    };
    std::vector<std::thread> threads;
    for (size_t t = 1; t < n; ++t) threads.emplace_back(guarded, t);
    if (n) guarded(0);
    for (auto &th : threads) th.join();
    for (auto &e : errors)
        if (e) std::rethrow_exception(e);
}

inline size_t workers_for(size_t n_items, int n_threads, size_t min_per_worker = 1 << 16) {
    const size_t t = static_cast<size_t>(std::max(1, n_threads));
    return std::max<size_t>(1, std::min(t, n_items / min_per_worker));
}

// Chunks sorted in parallel, then merged pairwise (in parallel per round).
// With `stable`, equal elements keep their order (std::merge takes the
// first range's element of a tie), as std::stable_sort.
template <class Elem, class Cmp>
void parallel_sort(std::vector<Elem> &v, Cmp cmp, int n_threads, bool stable = false) {
    const size_t n = v.size();
    const size_t T = workers_for(n, n_threads);
    auto sort = [&](auto b, auto e) {
        if (stable) std::stable_sort(b, e, cmp);
        else std::sort(b, e, cmp);
    };
    if (T <= 1) {
        sort(v.begin(), v.end());
        return;
    }
    std::vector<size_t> bounds(T + 1);
    for (size_t t = 0; t <= T; ++t) bounds[t] = n * t / T;
    run_parallel(T, [&](size_t t) { sort(v.begin() + bounds[t], v.begin() + bounds[t + 1]); });
    std::vector<Elem> tmp(n);
    std::vector<Elem> *src = &v, *dst = &tmp;
    while (bounds.size() > 2) {
        const size_t pairs = (bounds.size() - 1) / 2;
        const bool odd = (bounds.size() - 1) % 2 == 1;
        std::vector<size_t> next;
        for (size_t p = 0; p < pairs; ++p) next.push_back(bounds[2 * p]);
        if (odd) next.push_back(bounds[bounds.size() - 2]);
        next.push_back(n);
        run_parallel(pairs + (odd ? 1 : 0), [&](size_t p) {
            if (p == pairs) {  // the odd chunk out, copied
                std::copy(src->begin() + bounds[2 * p], src->end(), dst->begin() + bounds[2 * p]);
                return;
            }
            const size_t lo = bounds[2 * p], mid = bounds[2 * p + 1], hi = bounds[2 * p + 2];
            std::merge(src->begin() + lo, src->begin() + mid, src->begin() + mid,
                       src->begin() + hi, dst->begin() + lo, cmp);
        });
        bounds.swap(next);
        std::swap(src, dst);
    }
    if (src != &v) v.swap(*src);
}

// Writes rows [0, n_rows) to fh: format(out, begin, end) appends the text of
// rows [begin, end) to an empty `out`; rounds of up to n_threads blocks of
// block_rows rows are formatted in parallel and written in order, so the
// file does not depend on the thread count. Exceptions of format propagate
// (the caller closes fh). Returns false on a write error.
template <class Format>
bool write_blocks(std::FILE *fh, size_t n_rows, int n_threads, Format format,
                  size_t block_rows = size_t(1) << 15) {
    const size_t n_blocks = (n_rows + block_rows - 1) / block_rows;
    const size_t workers = std::max<size_t>(
        1, std::min<size_t>(static_cast<size_t>(std::max(1, n_threads)), n_blocks));
    std::vector<std::string> bufs(workers);
    for (size_t first = 0; first < n_blocks; first += workers) {
        const size_t count = std::min(workers, n_blocks - first);
        run_parallel(count, [&](size_t w) {
            const size_t begin = (first + w) * block_rows;
            bufs[w].clear();
            format(bufs[w], begin, std::min(n_rows, begin + block_rows));
        });
        for (size_t w = 0; w < count; ++w)
            if (std::fwrite(bufs[w].data(), 1, bufs[w].size(), fh) != bufs[w].size())
                return false;
    }
    return true;
}

}  // namespace gcen_parallel
