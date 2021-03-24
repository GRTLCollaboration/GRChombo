/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMD_HPP_
#define SIMD_HPP_

// Chombo includes
#include "MayDay.H"

// Other includes
#include "AlwaysInline.hpp"
#include <cmath>

// Chombo namespace
#include "UsingNamespace.H"

// This struct can be used to switch between simd and non-simd versions of the
// same function by overloading
struct disable_simd
{
};

/// simd_traits holds information about the current simd usage
// This is the base template,i.e. the fallback for when there is no SIMD support
// for a data type. The fallback is to just not use simd. The datatypes that
// have simd support correspond to template specialisations, which are defined
// in the architecture-specific files and simd_base.hpp
template <typename t> struct simd_traits
{
    using data_t = t;
    using mask_t = bool;
    static const int simd_len = 1;
};

/// This class represents a "vector" for a vectorised architecture
// Base template type: fallback for when there is no SIMD support for a data
// type The fallback is to just not use simd. Hence, all the member function
// definitions of simd below are pretty trivial. The datatypes that have simd
// support correspond to template specialisations, which are defined in the
// architecture-specific files and simd_base.hpp
template <typename t> struct simd
{
    static const int simd_len = simd_traits<t>::simd_len;
    t m_value;

    ALWAYS_INLINE
    simd() : m_value() {}

    ALWAYS_INLINE
    simd(const t &x) : m_value(x) {}

    ALWAYS_INLINE
    operator t &() { return m_value; }

    ALWAYS_INLINE
    operator const t &() const { return m_value; }

    ALWAYS_INLINE
    static simd load(const double *ptr) { return *ptr; }

    ALWAYS_INLINE
    static void store(double *ptr, const simd &a) { *ptr = a.m_value; }

    //   ALWAYS_INLINE
    //  simd &operator+=(const simd &a)
    //  {
    //      m_value += a.m_value;
    //      return *this;
    //  }
    //  ALWAYS_INLINE
    //  simd &operator-=(const simd &a)
    //  {
    //      m_value -= a.m_value;
    //      return *this;
    //  }
    //  ALWAYS_INLINE
    //  simd &operator*=(const simd &a)
    //  {
    //      m_value *= a.m_value;
    //      return *this;
    //  }
    //  ALWAYS_INLINE
    //  simd &operator/=(const simd &a)
    //  {
    //      m_value /= a.m_value;
    //      return *this;
    //  }

    ALWAYS_INLINE
    t operator[](int index) const { return m_value; }

    template <typename op_t> ALWAYS_INLINE simd foreach (op_t op) const
    {
        return simd(op(m_value));
    }

    template <typename op_t> ALWAYS_INLINE simd foreach (op_t op, t arg) const
    {
        return simd(op(m_value, arg));
    }

    friend simd simd_conditional(const bool cond, const simd &true_value,
                                 const simd &false_value)
    {
        return cond ? true_value : false_value;
    }

    friend ALWAYS_INLINE bool simd_compare_lt(const simd &a, const simd &b)
    {
        return a < b;
    }

    friend ALWAYS_INLINE bool simd_compare_gt(const simd &a, const simd &b)
    {
        return a > b;
    }

    friend ALWAYS_INLINE simd simd_min(const simd &a, const simd &b)
    {
        return (a <= b) ? a : b;
    }

    friend ALWAYS_INLINE simd simd_max(const simd &a, const simd &b)
    {
        return (a > b) ? a : b;
    }
};

// Define all the simd-functions whose implementation does not depend on the
// architecture
#include "simd_base.hpp"

// Define simd-functions whose implementation depends on the architecture
#if defined(__x86_64__)
#include "x64/x64.hpp"
#elif defined(__aarch64__)
#include "arm/arm.hpp"
#endif

// We have defined various simd-specific calls (simd_compare_lt,
// simd_compare_gt, min, max etc.)  For simd<t> these are defined in the various
// architecture-specific implementations.  Here, we make sure that the same
// function calls also work when simd is switched off. (e.g.
// simd_compare_lt(double, double) must work but needs to be implemented
template <typename t>
t simd_conditional(const bool cond, const t &true_value, const t &false_value)
{
    return cond ? true_value : false_value;
}

template <typename t> ALWAYS_INLINE bool simd_compare_lt(const t &a, const t &b)
{
    return a < b;
}

template <typename t> ALWAYS_INLINE bool simd_compare_gt(const t &a, const t &b)
{
    return a > b;
}

template <typename t> ALWAYS_INLINE t simd_min(const t &a, const t &b)
{
    return (a <= b) ? a : b;
}

template <typename t> ALWAYS_INLINE t simd_max(const t &a, const t &b)
{
    return (a > b) ? a : b;
}

//<-- End: Defining the simd specific calls for non-simd datatypes.

#include "simdify.hpp"

#endif /* SIMD_HPP_ */
