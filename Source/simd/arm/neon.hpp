#ifndef SIMD_NEON_HPP_
#define SIMD_NEON_HPP_

#include <arm_neon.h>

template <> struct simd_traits<double>
{
    typedef float64x2_t data_t;
    typedef uint64x2_t mask_t;
    static const int simd_len = 2;
};

template <> struct simd<double> : public simd_base<double>
{
    typedef typename simd_traits<double>::data_t data_t;
    typedef typename simd_traits<double>::mask_t mask_t;

    ALWAYS_INLINE
    simd() : simd_base<double>(vdupq_n_f64(0.)){ }

    ALWAYS_INLINE
    simd(double x) : simd_base<double>(vdupq_n_f64 (x)) {}

    ALWAYS_INLINE
    simd(const data_t &v) : simd_base<double>(v) {}

    ALWAYS_INLINE
    static simd load(const double *ptr) { return vld1q_f64(ptr); }

    ALWAYS_INLINE
    static void store(double *ptr, const simd &a)
    {
         vst1q_f64(ptr, a.m_value);
    }

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
        return vbslq_f64(cond, false_value, true_value);
    }

    friend ALWAYS_INLINE mask_t simd_compare_lt(const simd &a, const simd &b)
    {
        return vcltq_f64(a, b);
    }

    friend ALWAYS_INLINE mask_t simd_compare_gt(const simd &a, const simd &b)
    {
        return  vcgtq_f64(a, b);
    }

    friend ALWAYS_INLINE simd simd_min(const simd &a, const simd &b)
    {
        return vminq_f64(a, b);
    }

    friend ALWAYS_INLINE simd simd_max(const simd &a, const simd &b)
    {
        return vmaxq_f64(a, b);
    }

};


#endif
