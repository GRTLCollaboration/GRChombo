/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMD_BASE_HPP_
#define SIMD_BASE_HPP_

#if !defined(SIMD_HPP_)
#error "This file should only be included through simd.hpp"
#endif

#include <cmath>
#include <functional>
#include <iostream>
#include <type_traits>

// This file provides functionalities for simd<t> whose implementation is
// architecture independent.  Most importantly, all simd functions which have no
// associated vector intrinsics, e.g. calls like ostream operator <<,  are
// defined in the same way for all architectures and are, therefore, in this
// file.  A more subtle example:
//*= is an operator which has an architecture dependent vector intrinsic
// associated to it  and is, therefore, not defined in this file. However, the
// binary operation a*b can be rewritten  as { a *= b; return a; } for all
// architectures and, therefore, is defined in this file.
template <typename t> struct simd_base
{
    static const int simd_len = simd_traits<t>::simd_len;
    using simd_data_t = typename simd_traits<t>::data_t;
    simd_data_t m_value;

    ALWAYS_INLINE
    simd_base() : m_value() {}

    ALWAYS_INLINE
    simd_base(const simd_data_t &x) : m_value(x) {}

    ALWAYS_INLINE
    operator simd_data_t &() { return m_value; }

    ALWAYS_INLINE
    operator const simd_data_t &() const { return m_value; }

    // Note: These binary ops allow us to write e.g. simd<double>+int
    // or 2*y etc. The "friend function" construction makes
    // sure that these non-member binary operations
    // are instantiated at the same time as simd<t>.
    // The alternative is to template
    // on the two operand types, however this gives rise to
    // ambiguity because we have defined
    // casts between t and simd<t>
    friend ALWAYS_INLINE simd<t> operator+(const simd<t> &a, const simd<t> &b)
    {
        simd<t> out(a);
        out += b;
        return out;
    }

    friend ALWAYS_INLINE simd<t> operator-(const simd<t> &a, const simd<t> &b)
    {
        simd<t> out(a);
        out -= b;
        return out;
    }

    friend ALWAYS_INLINE simd<t> operator-(const simd<t> &a)
    {
        simd<t> out(0);
        out -= a;
        return out;
    }

    friend ALWAYS_INLINE simd<t> operator*(const simd<t> &a, const simd<t> &b)
    {
        simd<t> out(a);
        out *= b;
        return out;
    }

    friend ALWAYS_INLINE simd<t> operator/(const simd<t> &a, const simd<t> &b)
    {
        simd<t> out(a);
        out /= b;
        return out;
    }

    ALWAYS_INLINE
    t operator[](int index) const { return m_value[index]; }

    template <typename op_t> ALWAYS_INLINE simd<t> foreach (op_t op) const
    {
        t in_arr[simd_traits<t>::simd_len];
        t out_arr[simd_traits<t>::simd_len];
        simd<t>::store(in_arr, m_value);

#pragma omp simd
        for (int i = 0; i < simd_traits<t>::simd_len; ++i)
        {
            out_arr[i] = op(in_arr[i]);
        }

        return simd<t>::load(out_arr);
    }

    template <typename op_t>
    ALWAYS_INLINE simd<t> foreach (op_t op, simd<t> arg) const
    {
        t in_arr[simd_traits<t>::simd_len];
        t arg_arr[simd_traits<t>::simd_len];
        t out_arr[simd_traits<t>::simd_len];
        simd<t>::store(in_arr, m_value);
        simd<t>::store(arg_arr, arg);

#pragma omp simd
        for (int i = 0; i < simd_traits<t>::simd_len; ++i)
        {
            out_arr[i] = op(in_arr[i], arg_arr[i]);
        }

        return simd<t>::load(out_arr);
    }
};

#define define_simd_overload(op)                                               \
    template <typename t> ALWAYS_INLINE simd<t> op(const simd<t> &a)           \
    {                                                                          \
        return a.foreach (([&](t x) { return op(x); }));                       \
    }

#define define_binary_simd_overload(op)                                        \
    template <typename t>                                                      \
    ALWAYS_INLINE simd<t> op(const simd<t> &a, const simd<t> &b)               \
    {                                                                          \
        return a.foreach (([&](t x, t arg) { return op(x, arg); }), b);        \
    }

/* Trascendental support:                               */
/* exp, sin, cos, log, sqrt, pow, tanh, tan, sinh, cosh */

define_simd_overload(exp) define_simd_overload(exp2) define_simd_overload(sin)
    define_simd_overload(cos) define_simd_overload(log)
        define_simd_overload(log2) define_simd_overload(sqrt)
            define_binary_simd_overload(pow) define_simd_overload(abs)
                define_simd_overload(tanh) define_simd_overload(tan)
                    define_simd_overload(sinh) define_simd_overload(cosh)
                        define_simd_overload(acos) define_simd_overload(asin)
                            define_simd_overload(atan)
                                define_binary_simd_overload(atan2)

    /* Extra pow overloads */
    template <typename t, typename t1>
    ALWAYS_INLINE simd<t> pow(const simd<t> &a, const t1 b)
{
    simd<t> simd_b(b);
    return pow(a, simd_b);
}

/* Extra atan2 overloads */
template <typename t, typename t1>
ALWAYS_INLINE simd<t> atan2(const t1 b, const simd<t> &a)
{
    simd<t> simd_b(b);
    return atan2(simd_b, a);
}

template <typename t>
ALWAYS_INLINE std::ostream &operator<<(std::ostream &os, const simd<t> &in_simd)
{
    t in_arr[simd_traits<t>::simd_len];
    simd<t>::store(in_arr, in_simd);

    os << "( ";
    for (int i = 0; i < simd_traits<t>::simd_len; ++i)
    {
        os << in_simd[i] << " ";
    }
    os << ")";
    if (os.fail())
        MayDay::Error("operator<<(std::ostream&,simd<t>&) failed");
    return os;
}

#endif /* SIMD_BASE_HPP_ */
