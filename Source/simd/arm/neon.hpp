#ifndef SIMD_NEON_HPP_
#define SIMD_NEON_HPP_

#if !defined(SIMD_ARM_HPP_)
#error "This file should only be included through arm.hpp"
#endif

#ifdef __ARM_NEON

#include <arm_neon.h>

template <> struct simd_traits<double>
{
    typedef float64x2_t data_t;
    typedef uint64x2_t mask_t;
    static const int simd_len = 2;
};

template <> struct simd_traits<float>
{
    typedef float32x4_t data_t;
    typedef uint32x4_t mask_t;
    static const int simd_len = 4;
};

template <> struct simd<double> : public simd_base<double>
{
    typedef typename simd_traits<double>::data_t data_t;
    typedef typename simd_traits<double>::mask_t mask_t;

    ALWAYS_INLINE
    simd() : simd_base<double>(vdupq_n_f64(0.)) {}

    ALWAYS_INLINE
    simd(double x) : simd_base<double>(vdupq_n_f64(x)) {}

    ALWAYS_INLINE
    simd(const data_t &v) : simd_base<double>(v) {}

    ALWAYS_INLINE
    static simd load(const double *ptr) { return vld1q_f64(ptr); }

    ALWAYS_INLINE
    static void store(double *ptr, const simd &a) { vst1q_f64(ptr, a.m_value); }

    ALWAYS_INLINE
    simd &operator+=(const simd &a)
    {
        m_value = vaddq_f64(m_value, a.m_value);
        return *this;
    }
    ALWAYS_INLINE
    simd &operator-=(const simd &a)
    {
        m_value = vsubq_f64(m_value, a.m_value);
        return *this;
    }
    ALWAYS_INLINE
    simd &operator*=(const simd &a)
    {
        m_value = vmulq_f64(m_value, a.m_value);
        return *this;
    }
    ALWAYS_INLINE
    simd &operator/=(const simd &a)
    {
        m_value = vdivq_f64(m_value, a.m_value);
        return *this;
    }

    friend ALWAYS_INLINE simd simd_conditional(const mask_t cond,
                                               const simd &true_value,
                                               const simd &false_value)
    {
        return vbslq_f64(cond, true_value, false_value);
    }

    friend ALWAYS_INLINE mask_t simd_compare_lt(const simd &a, const simd &b)
    {
        return vcltq_f64(a, b);
    }

    friend ALWAYS_INLINE mask_t simd_compare_gt(const simd &a, const simd &b)
    {
        return vcgtq_f64(a, b);
    }

    friend ALWAYS_INLINE simd simd_min(const simd &a, const simd &b)
    {
        return vminq_f64(a, b);
    }

    friend ALWAYS_INLINE simd simd_max(const simd &a, const simd &b)
    {
        return vmaxq_f64(a, b);
    }

    friend ALWAYS_INLINE simd simd_sqrt(const simd &a) { return vsqrtq_f64(a); }
};

template <> struct simd<float> : public simd_base<float>
{
    typedef typename simd_traits<float>::data_t data_t;
    typedef typename simd_traits<float>::mask_t mask_t;

    ALWAYS_INLINE
    simd() : simd_base<float>(vdupq_n_f32(0.)) {}

    ALWAYS_INLINE
    simd(float x) : simd_base<float>(vdupq_n_f32(x)) {}

    ALWAYS_INLINE
    simd(const data_t &v) : simd_base<float>(v) {}

    ALWAYS_INLINE
    static simd load(const float *ptr) { return vld1q_f32(ptr); }

    ALWAYS_INLINE
    static void store(float *ptr, const simd &a) { vst1q_f32(ptr, a.m_value); }

    ALWAYS_INLINE
    simd &operator+=(const simd &a)
    {
        m_value = vaddq_f32(m_value, a.m_value);
        return *this;
    }
    ALWAYS_INLINE
    simd &operator-=(const simd &a)
    {
        m_value = vsubq_f32(m_value, a.m_value);
        return *this;
    }
    ALWAYS_INLINE
    simd &operator*=(const simd &a)
    {
        m_value = vmulq_f32(m_value, a.m_value);
        return *this;
    }
    ALWAYS_INLINE
    simd &operator/=(const simd &a)
    {
        m_value = vdivq_f32(m_value, a.m_value);
        return *this;
    }

    friend ALWAYS_INLINE simd simd_conditional(const mask_t cond,
                                               const simd &true_value,
                                               const simd &false_value)
    {
        return vbslq_f32(cond, true_value, false_value);
    }

    friend ALWAYS_INLINE mask_t simd_compare_lt(const simd &a, const simd &b)
    {
        return vcltq_f32(a, b);
    }

    friend ALWAYS_INLINE mask_t simd_compare_gt(const simd &a, const simd &b)
    {
        return vcgtq_f32(a, b);
    }

    friend ALWAYS_INLINE simd simd_min(const simd &a, const simd &b)
    {
        return vminq_f32(a, b);
    }

    friend ALWAYS_INLINE simd simd_max(const simd &a, const simd &b)
    {
        return vmaxq_f32(a, b);
    }

    friend ALWAYS_INLINE simd simd_sqrt(const simd &a) { return vsqrtq_f32(a); }
};

#endif /* __ARM_NEON */

#endif /* SIMD_NEON_HPP_*/
