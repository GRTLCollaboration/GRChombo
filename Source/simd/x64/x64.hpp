/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMD_X64_HPP_
#define SIMD_X64_HPP_

#if !defined(SIMD_HPP_)
#error "This file should only be included through simd.hpp"
#endif

#if defined(__AVX512F__)

#include "avx512.hpp"

#elif defined(__AVX__)

#include "avx.hpp"

#elif defined(__SSE2__)

#include "sse.hpp"

#else

#error "SSE2 must not be disabled on the x86-64 platform"

#endif

#endif /* SIMD_X64_HPP_ */
