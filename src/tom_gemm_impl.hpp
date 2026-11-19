// Kernel of tom_gemm.hpp, included by tom_gemm.cpp (baseline: SSE2, NEON or
// plain C++) and tom_gemm_avx.cpp (compiled with -mavx) inside an unnamed
// namespace, after tom_gemm_includes.hpp. Both must be compiled without
// floating-point contraction.

template <class T>
struct Lanes;

template <>
struct Lanes<double> {  // 4 lanes
    static constexpr int W = 4;
#if defined(GCEN2_GEMM_AVX)
    struct V { __m256d a; };
    static V load(const double *p) { return {_mm256_loadu_pd(p)}; }
    static void store(double *p, V v) { _mm256_storeu_pd(p, v.a); }
    static V madd(V c, V x, V y) { return {_mm256_add_pd(c.a, _mm256_mul_pd(x.a, y.a))}; }
#elif defined(GCEN2_GEMM_SSE2)
    struct V { __m128d a, b; };
    static V load(const double *p) { return {_mm_loadu_pd(p), _mm_loadu_pd(p + 2)}; }
    static void store(double *p, V v) {
        _mm_storeu_pd(p, v.a);
        _mm_storeu_pd(p + 2, v.b);
    }
    static V madd(V c, V x, V y) {
        return {_mm_add_pd(c.a, _mm_mul_pd(x.a, y.a)), _mm_add_pd(c.b, _mm_mul_pd(x.b, y.b))};
    }
#elif defined(GCEN2_GEMM_NEON)
    struct V { float64x2_t a, b; };
    static V load(const double *p) { return {vld1q_f64(p), vld1q_f64(p + 2)}; }
    static void store(double *p, V v) {
        vst1q_f64(p, v.a);
        vst1q_f64(p + 2, v.b);
    }
    static V madd(V c, V x, V y) {  // vmulq + vaddq: no fused vfmaq
        return {vaddq_f64(c.a, vmulq_f64(x.a, y.a)), vaddq_f64(c.b, vmulq_f64(x.b, y.b))};
    }
#else
    struct V { double l[4]; };
    static V load(const double *p) { V v; std::memcpy(v.l, p, sizeof v.l); return v; }
    static void store(double *p, V v) { std::memcpy(p, v.l, sizeof v.l); }
    static V madd(V c, V x, V y) {
        for (int i = 0; i < 4; ++i) {
            const double t = x.l[i] * y.l[i];
            c.l[i] = c.l[i] + t;
        }
        return c;
    }
#endif
};

template <>
struct Lanes<float> {  // 8 lanes
    static constexpr int W = 8;
#if defined(GCEN2_GEMM_AVX)
    struct V { __m256 a; };
    static V load(const float *p) { return {_mm256_loadu_ps(p)}; }
    static void store(float *p, V v) { _mm256_storeu_ps(p, v.a); }
    static V madd(V c, V x, V y) { return {_mm256_add_ps(c.a, _mm256_mul_ps(x.a, y.a))}; }
#elif defined(GCEN2_GEMM_SSE2)
    struct V { __m128 a, b; };
    static V load(const float *p) { return {_mm_loadu_ps(p), _mm_loadu_ps(p + 4)}; }
    static void store(float *p, V v) {
        _mm_storeu_ps(p, v.a);
        _mm_storeu_ps(p + 4, v.b);
    }
    static V madd(V c, V x, V y) {
        return {_mm_add_ps(c.a, _mm_mul_ps(x.a, y.a)), _mm_add_ps(c.b, _mm_mul_ps(x.b, y.b))};
    }
#elif defined(GCEN2_GEMM_NEON)
    struct V { float32x4_t a, b; };
    static V load(const float *p) { return {vld1q_f32(p), vld1q_f32(p + 4)}; }
    static void store(float *p, V v) {
        vst1q_f32(p, v.a);
        vst1q_f32(p + 4, v.b);
    }
    static V madd(V c, V x, V y) {
        return {vaddq_f32(c.a, vmulq_f32(x.a, y.a)), vaddq_f32(c.b, vmulq_f32(x.b, y.b))};
    }
#else
    struct V { float l[8]; };
    static V load(const float *p) { V v; std::memcpy(v.l, p, sizeof v.l); return v; }
    static void store(float *p, V v) { std::memcpy(p, v.l, sizeof v.l); }
    static V madd(V c, V x, V y) {
        for (int i = 0; i < 8; ++i) {
            const float t = x.l[i] * y.l[i];
            c.l[i] = c.l[i] + t;
        }
        return c;
    }
#endif
};

constexpr int kMR = 6, kNR = 2;  // register block: rows x columns
constexpr int kKC = 256;         // terms per pass over the panel
constexpr int kNC = 128;         // columns per panel
// Per pass the panel's kKC-term slices of up to kRC rows and kNC columns stay
// in cache (~0.5 MB for double) while every register block of the panel
// takes its kKC terms; the lane sums of the panel live in a buffer.
constexpr int kRC = 128;

// lanes of an MR x NR block += terms [k0, k1) (k1 - k0 a multiple of W)
template <class T, int MR, int NR>
inline void micro(const T *const *a, const T *const *b, int k0, int k1, T *acc) {
    using L = Lanes<T>;
    constexpr int W = L::W;
    typename L::V c[MR][NR];
    for (int i = 0; i < MR; ++i)
        for (int j = 0; j < NR; ++j) c[i][j] = L::load(acc + (i * NR + j) * W);
    for (int k = k0; k < k1; k += W) {
        typename L::V av[MR];
        for (int i = 0; i < MR; ++i) av[i] = L::load(a[i] + k);
        for (int j = 0; j < NR; ++j) {
            const typename L::V bv = L::load(b[j] + k);
            for (int i = 0; i < MR; ++i) c[i][j] = L::madd(c[i][j], av[i], bv);
        }
    }
    for (int i = 0; i < MR; ++i)
        for (int j = 0; j < NR; ++j) L::store(acc + (i * NR + j) * W, c[i][j]);
}

template <class T>
inline T reduce_lanes(T *s) {
    if (Lanes<T>::W == 8)
        for (int l = 0; l < 4; ++l) s[l] = s[l] + s[l + 4];
    return (s[0] + s[1]) + (s[2] + s[3]);
}

template <class T>
void gemm_nt_impl(const T *A, int64_t lda, int64_t r0, int R, int64_t c0, int C, int K, T *out) {
    constexpr int W = Lanes<T>::W, MR = kMR, NR = kNR, BS = MR * NR * W;
    const int KW = K / W * W;  // the rest goes to lanes 0.. one by one
    const int nbr = (kRC + MR - 1) / MR, nbc = (kNC + NR - 1) / NR;
    std::vector<T> acc(static_cast<size_t>(nbr) * nbc * BS);
    for (int q0 = 0; q0 < R; q0 += kRC) {
        const int rc = std::min(kRC, R - q0);
        for (int p0 = 0; p0 < C; p0 += kNC) {
            const int pc = std::min(kNC, C - p0);
            std::fill(acc.begin(), acc.end(), T(0));
            for (int k0 = 0; k0 < KW; k0 += kKC) {
                const int k1 = std::min(KW, k0 + kKC);
                for (int i0 = 0; i0 < rc; i0 += MR) {
                    const int mr = std::min(MR, rc - i0);
                    const T *a[MR];  // short blocks repeat their last row (unused)
                    for (int i = 0; i < MR; ++i)
                        a[i] = A + (r0 + q0 + i0 + std::min(i, mr - 1)) * lda;
                    for (int j0 = 0; j0 < pc; j0 += NR) {
                        const int nr = std::min(NR, pc - j0);
                        const T *b[NR];
                        for (int j = 0; j < NR; ++j)
                            b[j] = A + (c0 + p0 + j0 + std::min(j, nr - 1)) * lda;
                        micro<T, MR, NR>(
                            a, b, k0, k1,
                            &acc[(static_cast<size_t>(i0 / MR) * nbc + j0 / NR) * BS]);
                    }
                }
            }
            for (int i = 0; i < rc; ++i) {
                const T *ai = A + (r0 + q0 + i) * lda;
                for (int j = 0; j < pc; ++j) {
                    const T *lanes = &acc[(static_cast<size_t>(i / MR) * nbc + j / NR) * BS +
                                          ((i % MR) * NR + j % NR) * W];
                    T s[W];
                    std::copy(lanes, lanes + W, s);
                    const T *bj = A + (c0 + p0 + j) * lda;
                    for (int k = KW; k < K; ++k) {
                        const T t = ai[k] * bj[k];
                        s[k - KW] = s[k - KW] + t;
                    }
                    out[static_cast<size_t>(q0 + i) * C + p0 + j] = reduce_lanes(s);
                }
            }
        }
    }
}
