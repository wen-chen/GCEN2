#include "tom_gemm.hpp"

#include <cstdlib>
#include <cstring>

#include "tom_gemm_includes.hpp"

namespace {
#include "tom_gemm_impl.hpp"
}  // namespace

#if defined(GCEN2_AVX_KERNEL)
void gemm_nt_avx(const double *, int64_t, int64_t, int, int64_t, int, int, double *);
void gemm_nt_avx(const float *, int64_t, int64_t, int, int64_t, int, int, float *);

namespace {
bool use_avx() {
    static const bool yes = [] {
        const char *e = std::getenv("GCEN2_SIMD");
        if (e && std::strcmp(e, "baseline") == 0) return false;
        __builtin_cpu_init();
        return __builtin_cpu_supports("avx") != 0;
    }();
    return yes;
}
}  // namespace
#endif

template <class T>
void gemm_nt(const T *A, int64_t lda, int64_t r0, int R, int64_t c0, int C, int K, T *out) {
#if defined(GCEN2_AVX_KERNEL)
    if (use_avx()) return gemm_nt_avx(A, lda, r0, R, c0, C, K, out);
#endif
    gemm_nt_impl(A, lda, r0, R, c0, C, K, out);
}

const char *gemm_nt_kernel() {
#if defined(GCEN2_AVX_KERNEL)
    if (use_avx()) return "avx";
#endif
    return "baseline";
}

template void gemm_nt<double>(const double *, int64_t, int64_t, int, int64_t, int, int, double *);
template void gemm_nt<float>(const float *, int64_t, int64_t, int, int64_t, int, int, float *);
