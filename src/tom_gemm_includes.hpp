#pragma once
// Headers of tom_gemm_impl.hpp (which is included inside a namespace).
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>
// Which vector code the kernel uses (all give the same results): AVX, SSE2
// (also MSVC x64, which does not define __SSE2__), NEON (AArch64), or plain
// C++ (GCEN2_GEMM_SCALAR forces it, for tests).
#if !defined(GCEN2_GEMM_SCALAR)
#if defined(__AVX__)
#define GCEN2_GEMM_AVX 1
#include <immintrin.h>
#elif defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
#define GCEN2_GEMM_SSE2 1
#include <emmintrin.h>
#elif defined(__aarch64__) && defined(__ARM_NEON)
#define GCEN2_GEMM_NEON 1
#include <arm_neon.h>
#endif
#endif
