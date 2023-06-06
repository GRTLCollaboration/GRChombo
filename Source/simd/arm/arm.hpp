/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMD_ARM_HPP_
#define SIMD_ARM_HPP_

#if !defined(SIMD_HPP_)
#error "This file should only be included through simd.hpp"
#endif

#if defined(__ARM_FEATURE_SVE) && defined(__ARM_FEATURE_SVE_BITS)

#include "sve.hpp"

#elif defined(__ARM_NEON)

#include "neon.hpp"

#else

#warning                                                                       \
    "Not using any intrinsics on AArch64. If targetting SVE, you might need to enable SVE bits with a compiler flag e.g. -msve-vector-bits=... with GCC"

#endif

#endif /* SIMD_ARM_HPP_ */
