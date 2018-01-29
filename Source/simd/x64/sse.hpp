/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMD_SSE_HPP_
#define SIMD_SSE_HPP_

#if !defined(SIMD_X64_HPP_)
#error "This file should only be included through simd_x64.hpp"
#endif

#if defined(__SSE2__)

// Need SSE4.1 for blendvpd, blendvps
#if defined(__SSE4_1__)
#include <smmintrin.h>
#else
#include <emmintrin.h>
#endif

template <> struct simd_traits<double>
{
    typedef __m128d data_t;
    typedef __m128d mask_t;
    static const int simd_len = 2;
};

template <> struct simd_traits<float>
{
    typedef __m128 data_t;
    typedef __m128 mask_t;
    static const int simd_len = 4;
};

template <> struct simd<double> : public simd_base<double>
{
    typedef typename simd_traits<double>::data_t data_t;
    typedef typename simd_traits<double>::mask_t mask_t;

    ALWAYS_INLINE
    simd() : simd_base<double>(_mm_setzero_pd()) {}

    ALWAYS_INLINE
    simd(const double &s) : simd_base<double>(_mm_set1_pd(s)) {}

    ALWAYS_INLINE
    simd(const data_t &v) : simd_base<double>(v) {}

    ALWAYS_INLINE
    static simd load(const double *ptr) { return _mm_loadu_pd(ptr); }

    ALWAYS_INLINE
    static void store(double *ptr, const simd &a)
    {
        _mm_storeu_pd(ptr, a.m_value);
    }

    ALWAYS_INLINE
    simd &operator+=(const simd &a)
    {
        m_value = _mm_add_pd(m_value, a.m_value);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator-=(const simd &a)
    {
        m_value = _mm_sub_pd(m_value, a.m_value);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator*=(const simd &a)
    {
        m_value = _mm_mul_pd(m_value, a.m_value);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator/=(const simd &a)
    {
        m_value = _mm_div_pd(m_value, a.m_value);
        return *this;
    }

    friend ALWAYS_INLINE simd simd_conditional(const mask_t cond,
                                               const simd &true_value,
                                               const simd &false_value)
    {
#if defined(__SSE4_1__)
        return _mm_blendv_pd(false_value, true_value, cond);
#else
        double _cond[2], _true[2], _false[2], _blend[2];
        _mm_storeu_pd(&_cond[0], cond);
        _mm_storeu_pd(&_true[0], true_value);
        _mm_storeu_pd(&_false[0], false_value);
        for (int i = 0; i < 2; ++i)
            _blend[i] = _cond[i] ? _true[i] : _false[i];
        return _mm_loadu_pd(&_blend[0]);
#endif
    }

    friend ALWAYS_INLINE mask_t simd_compare_lt(const simd &a, const simd &b)
    {
        return _mm_cmplt_pd(a, b);
    }

    friend ALWAYS_INLINE mask_t simd_compare_gt(const simd &a, const simd &b)
    {
        return _mm_cmpgt_pd(a, b);
    }

    friend ALWAYS_INLINE simd simd_min(const simd &a, const simd &b)
    {
        return _mm_min_pd(a, b);
    }

    friend ALWAYS_INLINE simd simd_max(const simd &a, const simd &b)
    {
        return _mm_max_pd(a, b);
    }

    friend ALWAYS_INLINE simd simd_sqrt(const simd &a)
    {
        return _mm_sqrt_pd(a);
    }
};

template <> struct simd<float> : public simd_base<float>
{
    typedef typename simd_traits<float>::data_t data_t;
    typedef typename simd_traits<float>::mask_t mask_t;

    ALWAYS_INLINE
    simd() : simd_base<float>(_mm_setzero_ps()) {}

    ALWAYS_INLINE
    simd(const float &x) : simd_base<float>(_mm_set1_ps(x)) {}

    ALWAYS_INLINE
    simd(const data_t &v) : simd_base<float>(v) {}

    ALWAYS_INLINE
    static simd load(const float *ptr) { return _mm_loadu_ps(ptr); }

    ALWAYS_INLINE
    static void store(float *ptr, const simd &a)
    {
        _mm_storeu_ps(ptr, a.m_value);
    }

    ALWAYS_INLINE
    simd &operator+=(const simd &a)
    {
        m_value = _mm_add_ps(m_value, a);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator-=(const simd &a)
    {
        m_value = _mm_sub_ps(m_value, a);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator*=(const simd &a)
    {
        m_value = _mm_mul_ps(m_value, a);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator/=(const simd &a)
    {
        m_value = _mm_div_ps(m_value, a);
        return *this;
    }

    friend ALWAYS_INLINE simd simd_conditional(const mask_t cond,
                                               const simd &true_value,
                                               const simd &false_value)
    {
#if defined(__SSE4_1__)
        return _mm_blendv_ps(false_value, true_value, cond);
#else
        float _cond[4], _true[4], _false[4], _blend[4];
        _mm_storeu_ps(&_cond[0], cond);
        _mm_storeu_ps(&_true[0], true_value);
        _mm_storeu_ps(&_false[0], false_value);
        for (int i = 0; i < 4; ++i)
            _blend[i] = _cond[i] ? _true[i] : _false[i];
        return _mm_loadu_ps(&_blend[0]);
#endif
    }

    friend ALWAYS_INLINE mask_t simd_compare_lt(const simd &a, const simd &b)
    {
        return _mm_cmplt_ps(a, b);
    }

    friend ALWAYS_INLINE mask_t simd_compare_gt(const simd &a, const simd &b)
    {
        return _mm_cmpgt_ps(a, b);
    }

    friend ALWAYS_INLINE simd simd_min(const simd &a, const simd &b)
    {
        return _mm_min_ps(a, b);
    }

    friend ALWAYS_INLINE simd simd_max(const simd &a, const simd &b)
    {
        return _mm_max_ps(a, b);
    }

    friend ALWAYS_INLINE simd simd_sqrt(const simd &a)
    {
        return _mm_sqrt_ps(a);
    }
};

#endif /* __SSE2__ */

#endif /* SIMD_SSE_HPP_ */
