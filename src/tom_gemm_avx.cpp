// The AVX build of the tom_gemm kernel (this file is compiled with -mavx).
#include <cstdint>

#include "tom_gemm_includes.hpp"

namespace {
#include "tom_gemm_impl.hpp"
}  // namespace

void gemm_nt_avx(const double *A, int64_t lda, int64_t r0, int R, int64_t c0, int C, int K,
                 double *out) {
    gemm_nt_impl(A, lda, r0, R, c0, C, K, out);
}
void gemm_nt_avx(const float *A, int64_t lda, int64_t r0, int R, int64_t c0, int C, int K,
                 float *out) {
    gemm_nt_impl(A, lda, r0, R, c0, C, K, out);
}
