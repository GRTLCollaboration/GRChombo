/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMD_SVE_HPP_
#define SIMD_SVE_HPP_

#if !defined(SIMD_ARM_HPP_)
#error "This file should only be included through arm.hpp"
#endif

#if defined(__ARM_FEATURE_SVE) && defined(__ARM_FEATURE_SVE_BITS)

#include <arm_sve.h>

// Need to use this attribute in order to use the vector length specific
// instructions. These are necessary because of our implementation of simd
typedef svfloat64_t vecd
    __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));
typedef svfloat32_t vecf
    __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));
typedef svbool_t pred
    __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));

// Note that even though SVE enables the removal of remainder loops through the
// use of predicates (aka mask_t) in almost every intrinsic, we don't use these
// in order to maintain compatbility with our existing interface.

template <> struct simd_traits<double>
{
    typedef vecd data_t;
    typedef pred mask_t;
    static const int simd_len = sizeof(vecd) / sizeof(double);
};

template <> struct simd_traits<float>
{
    typedef vecf data_t;
    typedef pred mask_t;
    static const int simd_len = sizeof(vecf) / sizeof(float);
};

template <> struct simd<double> : public simd_base<double>
{
    typedef typename simd_traits<double>::data_t data_t;
    typedef typename simd_traits<double>::mask_t mask_t;

    ALWAYS_INLINE
    simd() : simd_base<double>(svdup_n_f64(0.)) {}

    ALWAYS_INLINE
    simd(const double &s) : simd_base<double>(svdup_n_f64(s)) {}

    ALWAYS_INLINE
    simd(const data_t &v) : simd_base<double>(v) {}

    ALWAYS_INLINE
    static simd load(const double *ptr)
    {
        return svld1_f64(svptrue_b64(), ptr);
    }

    ALWAYS_INLINE
    static void store(double *ptr, const simd &a)
    {
        svst1_f64(svptrue_b64(), ptr, a.m_value);
    }

    ALWAYS_INLINE
    simd &operator+=(const simd &a)
    {
        m_value = svadd_f64_z(svptrue_b64(), m_value, a.m_value);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator-=(const simd &a)
    {
        m_value = svsub_f64_z(svptrue_b64(), m_value, a.m_value);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator*=(const simd &a)
    {
        m_value = svmul_f64_z(svptrue_b64(), m_value, a.m_value);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator/=(const simd &a)
    {
        m_value = svdiv_f64_z(svptrue_b64(), m_value, a.m_value);
        return *this;
    }

    friend ALWAYS_INLINE simd simd_conditional(const mask_t cond,
                                               const simd &true_value,
                                               const simd &false_value)
    {
        return svsel_f64(cond, true_value, false_value);
    }

    friend ALWAYS_INLINE mask_t simd_compare_lt(const simd &a, const simd &b)
    {
        return svcmplt_f64(svptrue_b64(), a, b);
    }

    friend ALWAYS_INLINE mask_t simd_compare_gt(const simd &a, const simd &b)
    {
        return svcmpgt_f64(svptrue_b64(), a, b);
    }

    friend ALWAYS_INLINE simd simd_min(const simd &a, const simd &b)
    {
        return svmin_f64_z(svptrue_b64(), a, b);
    }

    friend ALWAYS_INLINE simd simd_max(const simd &a, const simd &b)
    {
        return svmax_f64_z(svptrue_b64(), a, b);
    }

    friend ALWAYS_INLINE simd simd_sqrt(const simd &a)
    {
        return svsqrt_f64_z(svptrue_b64(), a);
    }
};

template <> struct simd<float> : public simd_base<float>
{
    typedef typename simd_traits<float>::data_t data_t;
    typedef typename simd_traits<float>::mask_t mask_t;

    ALWAYS_INLINE
    simd() : simd_base<float>(svdup_n_f32(0.)) {}

    ALWAYS_INLINE
    simd(const float &s) : simd_base<float>(svdup_n_f32(s)) {}

    ALWAYS_INLINE
    simd(const data_t &v) : simd_base<float>(v) {}

    ALWAYS_INLINE
    static simd load(const float *ptr) { return svld1_f32(svptrue_b32(), ptr); }

    ALWAYS_INLINE
    static void store(float *ptr, const simd &a)
    {
        svst1_f32(svptrue_b32(), ptr, a.m_value);
    }

    ALWAYS_INLINE
    simd &operator+=(const simd &a)
    {
        m_value = svadd_f32_z(svptrue_b32(), m_value, a.m_value);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator-=(const simd &a)
    {
        m_value = svsub_f32_z(svptrue_b32(), m_value, a.m_value);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator*=(const simd &a)
    {
        m_value = svmul_f32_z(svptrue_b32(), m_value, a.m_value);
        return *this;
    }

    ALWAYS_INLINE
    simd &operator/=(const simd &a)
    {
        m_value = svdiv_f32_z(svptrue_b32(), m_value, a.m_value);
        return *this;
    }

    friend ALWAYS_INLINE simd simd_conditional(const mask_t cond,
                                               const simd &true_value,
                                               const simd &false_value)
    {
        return svsel_f32(cond, true_value, false_value);
    }

    friend ALWAYS_INLINE mask_t simd_compare_lt(const simd &a, const simd &b)
    {
        return svcmplt_f32(svptrue_b32(), a, b);
    }

    friend ALWAYS_INLINE mask_t simd_compare_gt(const simd &a, const simd &b)
    {
        return svcmpgt_f32(svptrue_b32(), a, b);
    }

    friend ALWAYS_INLINE simd simd_min(const simd &a, const simd &b)
    {
        return svmin_f32_z(svptrue_b32(), a, b);
    }

    friend ALWAYS_INLINE simd simd_max(const simd &a, const simd &b)
    {
        return svmax_f32_z(svptrue_b32(), a, b);
    }

    friend ALWAYS_INLINE simd simd_sqrt(const simd &a)
    {
        return svsqrt_f32_z(svptrue_b32(), a);
    }
};

#endif /* __ARM_FEATURE_SVE && __ARM_FEATURE_SVE_BITS */

#endif /* SIMD_SVE_HPP_ */
