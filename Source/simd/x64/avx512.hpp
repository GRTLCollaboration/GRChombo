/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMD_AVX512_HPP_
#define SIMD_AVX512_HPP_

#if !defined(SIMD_X64_HPP_)
#error "This file should only be included through simd_x64.hpp"
#endif

#if defined(__AVX512F__)

#include <immintrin.h>

template <> struct simd_traits<double>
{
    typedef __m512d data_t;
    typedef __mmask8 mask_t;
    static const int simd_len = 8;
};

template <> struct simd_traits<float>
{
    typedef __m512 data_t;
    typedef __mmask16 mask_t;
    static const int simd_len = 16;
};

template <> struct simd<double> : public simd_base<double>
{
    typedef typename simd_traits<double>::data_t data_t;
    typedef typename simd_traits<double>::mask_t mask_t;

    ALWAYS_INLINE
    simd() : simd_base<double>(_mm512_setzero_pd()) {}

    ALWAYS_INLINE
    simd(double x) : simd_base<double>(_mm512_set1_pd(x)) {}

    ALWAYS_INLINE
    simd(const data_t &v) : simd_base<double>(v) {}

    ALWAYS_INLINE
    static simd load(const double *ptr) { return _mm512_loadu_pd(ptr); }

    ALWAYS_INLINE
    static void store(double *ptr, const simd &a)
    {
        _mm512_storeu_pd(ptr, a.m_value);
    }

    ALWAYS_INLINE
    simd &operator+=(const simd &a)
    {
        m_value = _mm512_add_pd(m_value, a.m_value);
        return *this;
    }
    ALWAYS_INLINE
    simd &operator-=(const simd &a)
    {
        m_value = _mm512_sub_pd(m_value, a.m_value);
        return *this;
    }
    ALWAYS_INLINE
    simd &operator*=(const simd &a)
    {
        m_value = _mm512_mul_pd(m_value, a.m_value);
        return *this;
    }
    ALWAYS_INLINE
    simd &operator/=(const simd &a)
    {
        m_value = _mm512_div_pd(m_value, a.m_value);
        return *this;
    }

    friend ALWAYS_INLINE simd simd_conditional(const mask_t cond,
                                               const simd &true_value,
                                               const simd &false_value)
    {
        return _mm512_mask_blend_pd(cond, false_value, true_value);
    }

    friend ALWAYS_INLINE mask_t simd_compare_lt(const simd &a, const simd &b)
    {
        return _mm512_cmp_pd_mask(a, b, _CMP_LT_OQ);
    }

    friend ALWAYS_INLINE mask_t simd_compare_gt(const simd &a, const simd &b)
    {
        return _mm512_cmp_pd_mask(a, b, _CMP_GT_OQ);
    }

    friend ALWAYS_INLINE simd simd_min(const simd &a, const simd &b)
    {
        return _mm512_min_pd(a, b);
    }

    friend ALWAYS_INLINE simd simd_max(const simd &a, const simd &b)
    {
        return _mm512_max_pd(a, b);
    }

#ifdef __AVX512ER__
    friend ALWAYS_INLINE simd exp2(const simd &a)
    {
        return _mm512_exp2a23_pd(a);
    }

#ifdef LOW_PRECISION
    // This approximation has really low precision. Leaving it here mostly for
    // reference
    friend ALWAYS_INLINE simd exp(const simd &a)
    {
        // e^x = 2^(x * log2 (e))
        return exp2(a * simd(1.44269504088896340736));
    }
#endif /* LOW_PRECISION */
#endif /* __AVX512ER__ */

    friend ALWAYS_INLINE simd sqrt(const simd &a) { return _mm512_sqrt_pd(a); }
};

template <> struct simd<float> : public simd_base<float>
{
    typedef typename simd_traits<float>::data_t data_t;
    typedef typename simd_traits<float>::mask_t mask_t;

    ALWAYS_INLINE
    simd() : simd_base<float>(_mm512_setzero_ps()) {}

    ALWAYS_INLINE
    simd(float x) : simd_base<float>(_mm512_set1_ps(x)) {}

    ALWAYS_INLINE
    simd(const data_t &v) : simd_base<float>(v) {}

    ALWAYS_INLINE
    static simd load(const float *ptr) { return _mm512_loadu_ps(ptr); }

    ALWAYS_INLINE
    static void store(float *ptr, const simd &a)
    {
        _mm512_storeu_ps(ptr, a.m_value);
    }

    ALWAYS_INLINE
    simd &operator+=(const simd &a)
    {
        m_value = _mm512_add_ps(m_value, a);
        return *this;
    }
    ALWAYS_INLINE
    simd &operator-=(const simd &a)
    {
        m_value = _mm512_sub_ps(m_value, a);
        return *this;
    }
    ALWAYS_INLINE
    simd &operator*=(const simd &a)
    {
        m_value = _mm512_mul_ps(m_value, a);
        return *this;
    }
    ALWAYS_INLINE
    simd &operator/=(const simd &a)
    {
        m_value = _mm512_div_ps(m_value, a);
        return *this;
    }

    friend ALWAYS_INLINE simd simd_conditional(const mask_t cond,
                                               const simd &true_value,
                                               const simd &false_value)
    {
        return _mm512_mask_blend_ps(cond, false_value, true_value);
    }

    friend ALWAYS_INLINE mask_t simd_compare_lt(const simd &a, const simd &b)
    {
        return _mm512_cmp_ps_mask(a, b, _CMP_LT_OQ);
    }

    friend ALWAYS_INLINE mask_t simd_compare_gt(const simd &a, const simd &b)
    {
        return _mm512_cmp_ps_mask(a, b, _CMP_GT_OQ);
    }

    friend ALWAYS_INLINE simd simd_min(const simd &a, const simd &b)
    {
        return _mm512_min_ps(a, b);
    }

    friend ALWAYS_INLINE simd simd_max(const simd &a, const simd &b)
    {
        return _mm512_max_ps(a, b);
    }

#ifdef __AVX512ER__
    friend ALWAYS_INLINE simd exp2(const simd &a)
    {
        return _mm512_exp2a23_ps(a);
    }

#ifdef LOW_PRECISION
    // This approximation has really low precision. Leaving it here mostly for
    // reference
    friend ALWAYS_INLINE simd exp(const simd &a)
    {
        // e^x = 2^(x * log2 (e))
        return exp2(a * simd(1.44269504088896340736f));
    }
#endif /* LOW_PRECISION */
#endif /* __AVX512ER */

    friend ALWAYS_INLINE simd sqrt(const simd &a) { return _mm512_sqrt_ps(a); }
};

#endif /* __AVX512F__ */

#endif /* SIMD_AVX512_HPP_ */
