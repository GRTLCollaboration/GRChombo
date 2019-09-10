#ifndef SIMD_NEON_HPP_
#define SIMD_NEON_HPP_

#include <arm_neon.h>

template <> struct simd_traits<double>
{
    typedef float64x2_t data_t;
    static const int simd_len = 2;
};

template <> struct simd<double> : public simd_base<double>
{
    typedef typename simd_traits<double>::data_t data_t;

    ALWAYS_INLINE
    simd(double x) : simd_base<double>(vdupq_n_f64 (x)) {}

    ALWAYS_INLINE
    simd(const data_t &v) : simd_base<double>(v) {}

    ALWAYS_INLINE
    static simd load(const double *ptr) { return vld1q_f64(ptr); }

    ALWAYS_INLINE
    simd &operator+=(const simd &a)
    {
        m_value = vaddq_f64(m_value, a.m_value);
        return *this;
    }

    ALWAYS_INLINE
    static void store(double *ptr, const simd &a)
    {
         vst1q_f64(ptr, a.m_value);
    }
};


#endif
