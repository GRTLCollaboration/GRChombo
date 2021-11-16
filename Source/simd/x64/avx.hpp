/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMD_AVX_HPP_
#define SIMD_AVX_HPP_

#if !defined(SIMD_X64_HPP_)
#error "This file should only be included through simd_x64.hpp"
#endif

#if defined(__AVX__)

#include <immintrin.h>

template <> struct simd_traits<double>
{
    typedef __m256d data_t;
    typedef __m256d mask_t;
    static const int simd_len = 4;
};

template <> struct simd_traits<float>
{
    typedef __m256 data_t;
    typedef __m256 mask_t;
    static const int simd_len = 8;
};

template <> struct simd<double> : public simd_base<double>
{
    typedef typename simd_traits<double>::data_t data_t;
    typedef typename simd_traits<double>::mask_t mask_t;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuninitialized"
    ALWAYS_INLINE
    simd() : simd_base<double>(_mm256_setzero_pd()) {}
#pragma GCC diagnostic pop

    ALWAYS_INLINE
    simd(const double &s) : simd_base<double>(_mm256_set1_pd(s)) {}

    ALWAYS_INLINE
    simd(const data_t &v) : simd_base<double>(v) {}

    ALWAYS_INLINE
    static simd load(const double *ptr) { return _mm256_loadu_pd(ptr); }

    ALWAYS_INLINE
    static void store(double *ptr, const simd &a)
    {
        _mm256_storeu_pd(ptr, a.m_value);
    }

    ALWAYS_INLINE
    simd &operator+=(const simd &a)
    {
        m_value = _mm256_add_pd(m_value, a.m_value);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator-=(const simd &a)
    {
        m_value = _mm256_sub_pd(m_value, a.m_value);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator*=(const simd &a)
    {
        m_value = _mm256_mul_pd(m_value, a.m_value);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator/=(const simd &a)
    {
        m_value = _mm256_div_pd(m_value, a.m_value);
        return *this;
    }

    friend ALWAYS_INLINE simd simd_conditional(const mask_t cond,
                                               const simd &true_value,
                                               const simd &false_value)
    {
        return _mm256_blendv_pd(false_value, true_value, cond);
    }

    friend ALWAYS_INLINE mask_t simd_compare_lt(const simd &a, const simd &b)
    {
        return _mm256_cmp_pd(a, b, _CMP_LT_OQ);
    }

    friend ALWAYS_INLINE mask_t simd_compare_gt(const simd &a, const simd &b)
    {
        return _mm256_cmp_pd(a, b, _CMP_GT_OQ);
    }

    friend ALWAYS_INLINE simd simd_min(const simd &a, const simd &b)
    {
        return _mm256_min_pd(a, b);
    }

    friend ALWAYS_INLINE simd simd_max(const simd &a, const simd &b)
    {
        return _mm256_max_pd(a, b);
    }

    friend ALWAYS_INLINE simd simd_sqrt(const simd &a)
    {
        return _mm256_sqrt_pd(a);
    }
};

template <> struct simd<float> : public simd_base<float>
{
    typedef typename simd_traits<float>::data_t data_t;
    typedef typename simd_traits<float>::mask_t mask_t;

    ALWAYS_INLINE
    simd() : simd_base<float>(_mm256_setzero_ps()) {}

    ALWAYS_INLINE
    simd(const float &x) : simd_base<float>(_mm256_set1_ps(x)) {}

    ALWAYS_INLINE
    simd(const data_t &v) : simd_base<float>(v) {}

    ALWAYS_INLINE
    static simd load(const float *ptr) { return _mm256_loadu_ps(ptr); }

    ALWAYS_INLINE
    static void store(float *ptr, const simd &a)
    {
        _mm256_storeu_ps(ptr, a.m_value);
    }

    ALWAYS_INLINE
    simd &operator+=(const simd &a)
    {
        m_value = _mm256_add_ps(m_value, a);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator-=(const simd &a)
    {
        m_value = _mm256_sub_ps(m_value, a);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator*=(const simd &a)
    {
        m_value = _mm256_mul_ps(m_value, a);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator/=(const simd &a)
    {
        m_value = _mm256_div_ps(m_value, a);
        return *this;
    }

    friend ALWAYS_INLINE simd simd_conditional(const mask_t cond,
                                               const simd &true_value,
                                               const simd &false_value)
    {
        return _mm256_blendv_ps(false_value, true_value, cond);
    }

    friend ALWAYS_INLINE mask_t simd_compare_lt(const simd &a, const simd &b)
    {
        return _mm256_cmp_ps(a, b, _CMP_LT_OQ);
    }

    friend ALWAYS_INLINE mask_t simd_compare_gt(const simd &a, const simd &b)
    {
        return _mm256_cmp_ps(a, b, _CMP_GT_OQ);
    }

    friend ALWAYS_INLINE simd simd_min(const simd &a, const simd &b)
    {
        return _mm256_min_ps(a, b);
    }

    friend ALWAYS_INLINE simd simd_max(const simd &a, const simd &b)
    {
        return _mm256_max_ps(a, b);
    }

    friend ALWAYS_INLINE simd simd_sqrt(const simd &a)
    {
        return _mm256_sqrt_ps(a);
    }
};

#endif /* __AVX__ */

#endif /* SIMD_AVX_HPP_ */
