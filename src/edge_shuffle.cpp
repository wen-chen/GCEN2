#include "edge_shuffle.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <utility>

EdgeSwapper::EdgeSwapper(std::vector<int64_t> a, std::vector<int64_t> b, int64_t base)
    : e_(a.size()), base_(base), small_keys_(base <= 65535) {
    if (a.size() != b.size()) throw std::invalid_argument("a and b must have the same length.");
    if (base_ < 1 || base_ > (int64_t(1) << 31)) throw std::invalid_argument("bad node count.");
    for (size_t k = 0; k < a.size(); ++k) {
        if (a[k] < 0 || a[k] >= base_ || b[k] < 0 || b[k] >= base_)
            throw std::out_of_range("node code out of range.");
        e_[k] = Edge{static_cast<int32_t>(a[k]), static_cast<int32_t>(b[k])};
    }
    if (small_keys_) {
        keys32_.reserve(a.size());
        for (size_t k = 0; k < a.size(); ++k) keys32_.insert(static_cast<uint32_t>(key(a[k], b[k])));
    } else {
        keys64_.reserve(a.size());
        for (size_t k = 0; k < a.size(); ++k) keys64_.insert(key(a[k], b[k]));
    }
}

std::vector<int64_t> EdgeSwapper::b() const {
    std::vector<int64_t> v(e_.size());
    for (size_t k = 0; k < e_.size(); ++k) v[k] = e_[k].b;
    return v;
}

namespace {
inline uint64_t mix(uint64_t x) {  // splitmix64 finaliser
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
}  // namespace

template <typename Key>
void KeySet<Key>::rebuild(size_t cap) {
    HugeBuffer<Key> old(cap);
    old.swap(t_);
    t_.fill(kEmpty);
    mask_ = cap - 1;
    used_ = 0;
    for (size_t i = 0; i < old.size(); ++i)
        if (old[i] != kEmpty) insert(old[i]);
}

template <typename Key>
void KeySet<Key>::reserve(size_t n) {
    size_t cap = 16;
    while (cap < 4 * n + 4) cap <<= 1;  // load factor <= 1/4: short probe runs
    rebuild(cap);
}

template <typename Key>
size_t KeySet<Key>::slot(Key k) const { return mix(static_cast<uint64_t>(k)) & mask_; }

template <typename Key>
bool KeySet<Key>::contains(Key k) const {
    for (size_t p = slot(k);; p = (p + 1) & mask_) {
        if (t_[p] == k) return true;
        if (t_[p] == kEmpty) return false;
    }
}

template <typename Key>
void KeySet<Key>::insert(Key k) {
    if (t_.empty()) reserve(1);
    size_t p = slot(k);
    for (; t_[p] != kEmpty; p = (p + 1) & mask_)
        if (t_[p] == k) return;
    t_[p] = k;
    if (++used_ * 2 > t_.size()) rebuild(t_.size() * 2);
}

template <typename Key>
void KeySet<Key>::erase(Key k) {
    if (t_.empty()) return;
    size_t p = slot(k);
    for (; t_[p] != k; p = (p + 1) & mask_)
        if (t_[p] == kEmpty) return;
    // backward-shift: move later entries of the probe run into the hole
    // unless their home slot lies cyclically in (hole, entry]
    for (size_t q = (p + 1) & mask_; t_[q] != kEmpty; q = (q + 1) & mask_) {
        const size_t home = slot(t_[q]);
        if (((q - home) & mask_) >= ((q - p) & mask_)) {
            t_[p] = t_[q];
            p = q;
        }
    }
    t_[p] = kEmpty;
    --used_;
}

template class KeySet<uint32_t>;
template class KeySet<uint64_t>;

namespace {
// prefetch distances (attempts ahead) for edge records and key slots; on the
// 3.6M-edge sample edge 16..48 / key 4..16 are all within 3%
constexpr int64_t kEdgeAhead = 16;
constexpr int64_t kKeyAhead = 8;
}  // namespace

// Same steps as the Python loop in gcen2.network.shuffle.
int64_t EdgeSwapper::run(const int64_t *i, const int64_t *j, int64_t m) {
    const int64_t n = static_cast<int64_t>(e_.size());
    for (int64_t t = 0; t < m; ++t)
        if (i[t] < 0 || i[t] >= n || j[t] < 0 || j[t] >= n)
            throw std::out_of_range("edge index out of range.");
    return small_keys_ ? run_with(keys32_, i, j, m) : run_with(keys64_, i, j, m);
}

template <typename Key>
int64_t EdgeSwapper::run_with(KeySet<Key> &set, const int64_t *i, const int64_t *j, int64_t m) {
    auto k_of = [&](int64_t u, int64_t v) { return static_cast<Key>(key(u, v)); };
    int64_t applied = 0;
    for (int64_t t = 0; t < m; ++t) {
        if (t + kEdgeAhead < m) {  // the swap targets are known ahead: prefetch them
            GCEN_PREFETCH(&e_[i[t + kEdgeAhead]]);
            GCEN_PREFETCH(&e_[j[t + kEdgeAhead]]);
        }
        if (t + kKeyAhead < m) {  // and the key slots they would touch (if unchanged)
            const Edge ex = e_[i[t + kKeyAhead]], ey = e_[j[t + kKeyAhead]];
            set.prefetch(k_of(ex.a, ey.b));
            set.prefetch(k_of(ey.a, ex.b));
            set.prefetch(k_of(ex.a, ex.b));  // erased if the swap is applied
            set.prefetch(k_of(ey.a, ey.b));
        }
        const int64_t x = i[t], y = j[t];
        if (x == y) continue;
        const int64_t a1 = e_[x].a, b1 = e_[x].b, a2 = e_[y].a, b2 = e_[y].b;
        if (a1 == b2 || a2 == b1) continue;
        const Key k1 = k_of(a1, b2), k2 = k_of(a2, b1);
        if (k1 == k2 || set.contains(k1) || set.contains(k2)) continue;
        set.erase(k_of(a1, b1));
        set.erase(k_of(a2, b2));
        e_[x].b = static_cast<int32_t>(b2);
        e_[y].b = static_cast<int32_t>(b1);
        set.insert(k1);
        set.insert(k2);
        ++applied;
    }
    return applied;
}

int64_t EdgeSwapper::shuffle(Pcg64 &rng, int64_t n_swaps) {
    const int64_t n = static_cast<int64_t>(e_.size());
    if (n < 2 || n_swaps <= 0) return 0;
    constexpr int64_t kBatch = int64_t(1) << 20;
    std::vector<int64_t> i(static_cast<size_t>(std::min(kBatch, n_swaps)));
    std::vector<int64_t> j(i.size());
    int64_t applied = 0;
    for (int64_t done = 0; done < n_swaps;) {
        const int64_t m = std::min(kBatch, n_swaps - done);
        rng.integers(n, i.data(), m);
        rng.integers(n, j.data(), m);
        applied += run(i.data(), j.data(), m);
        done += m;
    }
    return applied;
}
