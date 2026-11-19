#pragma once
/**
 * Degree-preserving edge swaps for gcen2.network.shuffle. The random edge
 * indices are numpy's stream (drawn by numpy and passed to run(), or
 * reproduced in C++ by shuffle()), so a given seed gives the same network
 * as the pure-Python loop.
 */
#include <cstdint>
#include <vector>

#include "huge_buffer.hpp"
#include "pcg64.hpp"

#if defined(__GNUC__) || defined(__clang__)
#define GCEN_PREFETCH(p) __builtin_prefetch(p)
#else
#define GCEN_PREFETCH(p) ((void)0)
#endif

/** Set of keys below Key's maximum (which marks an empty slot): open
 *  addressing with linear probing; erase shifts the following entries back
 *  (no tombstones). */
template <typename Key>
class KeySet {
public:
    /** Room for n keys at a load factor of at most 1/4. */
    void reserve(size_t n);
    bool contains(Key k) const;
    void insert(Key k);
    void erase(Key k);
    void prefetch(Key k) const { GCEN_PREFETCH(&t_[slot(k)]); }

private:
    static constexpr Key kEmpty = static_cast<Key>(~Key(0));
    size_t slot(Key k) const;
    void rebuild(size_t cap);
    HugeBuffer<Key> t_;
    size_t mask_ = 0, used_ = 0;
};

class EdgeSwapper {
public:
    /** Edges (a[k], b[k]) over node codes in [0, base). */
    EdgeSwapper(std::vector<int64_t> a, std::vector<int64_t> b, int64_t base);
    /** Try the swaps (i[t], j[t]) in order; returns the number applied. */
    int64_t run(const int64_t *i, const int64_t *j, int64_t m);
    /** n_swaps swaps drawn as gcen2.network.shuffle draws them from numpy:
     *  per batch of up to 2^20, rng.integers(0, n, m) for i, then for j. */
    int64_t shuffle(Pcg64 &rng, int64_t n_swaps);
    /** Current second endpoint of every edge (the first never changes). */
    std::vector<int64_t> b() const;

private:
    struct Edge {
        int32_t a, b;
    };
    uint64_t key(int64_t u, int64_t v) const {
        return static_cast<uint64_t>(u < v ? u * base_ + v : v * base_ + u);
    }
    template <typename Key>
    int64_t run_with(KeySet<Key> &set, const int64_t *i, const int64_t *j, int64_t m);
    // both endpoints in one 8-byte record: one cache line per edge looked
    // up, and 3.6M edges (29 MB) stay mostly in L3
    HugeBuffer<Edge> e_;
    int64_t base_;
    // edge keys u * base + v (u < v): 32-bit when they fit (base <= 65535),
    // which halves the table
    bool small_keys_;
    KeySet<uint32_t> keys32_;
    KeySet<uint64_t> keys64_;
};
