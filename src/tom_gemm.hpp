#pragma once
/**
 * The matrix products of the TOM (TomCore in tom.cpp): out (R x C,
 * row-major) = rows [r0, r0 + R) of A times rows [c0, c0 + C) of A
 * transposed, over the first K columns (rows of A are lda apart).
 *
 * The summation order is fixed, so results are bit-identical on every
 * instruction set and for every tiling: element (i, j) is accumulated in W
 * interleaved partial sums (W = 4 for double, 8 for float; lane l takes the
 * terms k = l mod W in increasing k), each term a separately rounded
 * product added to its lane (no fused multiply-add), and the lanes are then
 * reduced pairwise ((s0 + s1) + (s2 + s3); for 8 lanes s_l + s_(l+4) first).
 * An AVX kernel is used when the CPU has AVX (x86-64 GCC/Clang builds), else
 * SSE2 (x86-64), NEON (AArch64) or plain C++; GCEN2_SIMD=baseline forces the
 * baseline kernel at run time.
 */
#include <cstdint>

template <class T>
void gemm_nt(const T *A, int64_t lda, int64_t r0, int R, int64_t c0, int C, int K, T *out);

// "avx" or "baseline": the kernel gemm_nt uses
const char *gemm_nt_kernel();
