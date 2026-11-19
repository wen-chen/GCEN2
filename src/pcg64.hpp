#pragma once
/**
 * numpy's PCG64 bit generator (PCG XSL-RR 128/64) and the bounded integers
 * of Generator.integers, so that C++ code can draw the same numbers as
 * `np.random.default_rng(seed).integers(0, high, size)`.
 *
 * The generator is started from a numpy state (bit_generator.state: the
 * 128-bit `state` and `inc`, and the buffered half of a 64-bit output,
 * `has_uint32` / `uinteger`), so seeding (SeedSequence) stays in numpy.
 */
#include <cstdint>
#include <stdexcept>

class Pcg64 {
public:
    Pcg64(uint64_t state_hi, uint64_t state_lo, uint64_t inc_hi, uint64_t inc_lo,
          bool has_uint32, uint32_t uinteger)
        : s_hi_(state_hi), s_lo_(state_lo), i_hi_(inc_hi), i_lo_(inc_lo),
          has_uint32_(has_uint32), uinteger_(uinteger) {}

    uint64_t next64() {
        step();
        // XSL-RR output: (hi ^ lo) rotated right by the top 6 bits
        const uint64_t x = s_hi_ ^ s_lo_;
        const unsigned rot = static_cast<unsigned>(s_hi_ >> 58);
        return (x >> rot) | (x << ((64u - rot) & 63u));
    }

    // numpy's pcg64_next32: low half of a 64-bit output, the high half kept
    // for the next call.
    uint32_t next32() {
        if (has_uint32_) {
            has_uint32_ = false;
            return uinteger_;
        }
        const uint64_t v = next64();
        has_uint32_ = true;
        uinteger_ = static_cast<uint32_t>(v >> 32);
        return static_cast<uint32_t>(v);
    }

    // Generator.integers(0, high, size=n) with the default int64 dtype
    // (numpy's random_bounded_uint64_fill, Lemire's method, off = 0).
    void integers(int64_t high, int64_t *out, int64_t n) {
        if (high < 1) throw std::invalid_argument("high must be >= 1.");
        const uint64_t rng = static_cast<uint64_t>(high) - 1;
        if (rng == 0) {
            for (int64_t k = 0; k < n; ++k) out[k] = 0;
        } else if (rng <= 0xFFFFFFFFULL) {
            if (rng == 0xFFFFFFFFULL) {
                for (int64_t k = 0; k < n; ++k) out[k] = next32();
            } else {
                const uint32_t r = static_cast<uint32_t>(rng);
                for (int64_t k = 0; k < n; ++k) out[k] = bounded32(r);
            }
        } else {
            for (int64_t k = 0; k < n; ++k) out[k] = static_cast<int64_t>(bounded64(rng));
        }
    }

    uint64_t state_hi() const { return s_hi_; }
    uint64_t state_lo() const { return s_lo_; }
    uint64_t inc_hi() const { return i_hi_; }
    uint64_t inc_lo() const { return i_lo_; }
    bool has_uint32() const { return has_uint32_; }
    uint32_t uinteger() const { return uinteger_; }

private:
    // state = state * M + inc (mod 2^128)
    void step() {
        constexpr uint64_t kMulHi = 2549297995355413924ULL;
        constexpr uint64_t kMulLo = 4865540595714422341ULL;
        uint64_t lo_hi, lo_lo;
        mul64(s_lo_, kMulLo, lo_hi, lo_lo);
        uint64_t hi = lo_hi + s_lo_ * kMulHi + s_hi_ * kMulLo;
        uint64_t lo = lo_lo + i_lo_;
        hi += i_hi_ + (lo < lo_lo ? 1 : 0);
        s_hi_ = hi;
        s_lo_ = lo;
    }

    // full 64 x 64 -> 128-bit product
    static void mul64(uint64_t a, uint64_t b, uint64_t &hi, uint64_t &lo) {
#if defined(__SIZEOF_INT128__)
        const unsigned __int128 m = static_cast<unsigned __int128>(a) * b;
        hi = static_cast<uint64_t>(m >> 64);
        lo = static_cast<uint64_t>(m);
#else
        const uint64_t a0 = a & 0xFFFFFFFFULL, a1 = a >> 32;
        const uint64_t b0 = b & 0xFFFFFFFFULL, b1 = b >> 32;
        const uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
        const uint64_t mid = (p00 >> 32) + (p01 & 0xFFFFFFFFULL) + (p10 & 0xFFFFFFFFULL);
        hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
        lo = (mid << 32) | (p00 & 0xFFFFFFFFULL);
#endif
    }

    // numpy's buffered_bounded_lemire_uint32: a value in [0, rng]
    uint32_t bounded32(uint32_t rng) {
        const uint32_t rng_excl = rng + 1;
        uint64_t m = static_cast<uint64_t>(next32()) * rng_excl;
        uint32_t leftover = static_cast<uint32_t>(m);
        if (leftover < rng_excl) {
            const uint32_t threshold = (UINT32_MAX - rng) % rng_excl;
            while (leftover < threshold) {
                m = static_cast<uint64_t>(next32()) * rng_excl;
                leftover = static_cast<uint32_t>(m);
            }
        }
        return static_cast<uint32_t>(m >> 32);
    }

    // numpy's bounded_lemire_uint64: a value in [0, rng], rng < 2^64 - 1
    uint64_t bounded64(uint64_t rng) {
        if (rng == UINT64_MAX) return next64();
        const uint64_t rng_excl = rng + 1;
        uint64_t hi, leftover;
        mul64(next64(), rng_excl, hi, leftover);
        if (leftover < rng_excl) {
            const uint64_t threshold = (UINT64_MAX - rng) % rng_excl;
            while (leftover < threshold) mul64(next64(), rng_excl, hi, leftover);
        }
        return hi;
    }

    uint64_t s_hi_, s_lo_, i_hi_, i_lo_;
    bool has_uint32_;
    uint32_t uinteger_;
};
