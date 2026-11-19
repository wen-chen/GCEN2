#pragma once
// A fixed-size array that asks Linux for transparent huge pages (2 MB):
// for large tables probed at random, 4 KB pages cost a TLB miss per probe.
// Elsewhere (or if the request fails) it is an ordinary allocation.
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <new>
#include <utility>
#if defined(__linux__)
#include <sys/mman.h>
#endif

template <typename T>
class HugeBuffer {
public:
    HugeBuffer() = default;
    explicit HugeBuffer(size_t n) : n_(n) {
        if (n == 0) return;
        const size_t bytes = n * sizeof(T);
        constexpr size_t kHuge = size_t(2) << 20;
#if defined(__linux__)
        // std::aligned_alloc is missing from MSVC's library and needs macOS
        // 10.15+, and huge pages are only requested on Linux anyway
        if (bytes >= kHuge) {
            const size_t rounded = (bytes + kHuge - 1) / kHuge * kHuge;
            p_ = static_cast<T *>(std::aligned_alloc(kHuge, rounded));
#if defined(MADV_HUGEPAGE)
            if (p_) madvise(p_, rounded, MADV_HUGEPAGE);
#endif
        }
#endif
        if (!p_) p_ = static_cast<T *>(std::malloc(bytes));
        if (!p_) throw std::bad_alloc();
    }
    ~HugeBuffer() { std::free(p_); }
    HugeBuffer(const HugeBuffer &) = delete;
    HugeBuffer &operator=(const HugeBuffer &) = delete;
    HugeBuffer(HugeBuffer &&o) noexcept { swap(o); }
    HugeBuffer &operator=(HugeBuffer &&o) noexcept { swap(o); return *this; }
    void swap(HugeBuffer &o) noexcept { std::swap(p_, o.p_); std::swap(n_, o.n_); }

    T &operator[](size_t i) { return p_[i]; }
    const T &operator[](size_t i) const { return p_[i]; }
    T *data() { return p_; }
    const T *data() const { return p_; }
    size_t size() const { return n_; }
    bool empty() const { return n_ == 0; }
    void fill(const T &v) { for (size_t i = 0; i < n_; ++i) p_[i] = v; }

private:
    T *p_ = nullptr;
    size_t n_ = 0;
};
