/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "REAL.H"

// Other includes
#include "simd.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

// Chombo namespace
#include "UsingNamespace.H"

// #define DEBUG 1

template <class t> bool similar(t a, t b, t *error, t *tolerance)
{
    *error = std::abs(a - b);
    *tolerance = 2 * std::numeric_limits<t>::epsilon() *
                 std::max(std::abs(a), std::abs(a));
    return *error <= *tolerance;
}

template <class t, class sop_t, class vop_t>
bool sv_test(const char *name, sop_t sop, vop_t vop)
{
    constexpr int simd_length = simd_traits<t>::simd_len;
    t vals[simd_length];

    for (int i = 0; i < simd_length; i++)
        vals[i] = i + 1;
    auto simd_in = simd<t>::load(vals);
    auto simd_out = vop(simd_in);

#ifdef __INTEL_COMPILER
#pragma novector
#else
#pragma omp simd safelen(1)
#endif /* __INTEL_COMPILER */
    for (int i = 0; i < simd_length; i++)
        vals[i] = sop(static_cast<t>(i + 1));

    for (int i = 0; i < simd_length; i++)
    {
        t error, tolerance;
        if (!similar(simd_out[i], vals[i], &error, &tolerance))
        {
#if DEBUG
            std::cout.precision(std::numeric_limits<t>::max_digits10 + 1);
            std::cout << name << " input=" << i + 1 << " scalar=" << vals[i]
                      << " vector=" << simd_out[i] << " error=" << error
                      << " tolerance=" << tolerance << std::endl;
#endif
            return true;
        }
    }

    return false;
}

template <class t, class op_t, class rev_op_t>
bool rv_test(const char *name, op_t op, rev_op_t rev_op)
{
    constexpr int simd_length = simd_traits<t>::simd_len;
    t vals[simd_length];

#ifdef __INTEL_COMPILER
#pragma novector
#else
#pragma omp simd safelen(1)
#endif /* __INTEL_COMPILER */
    for (int i = 0; i < simd_length; i++)
        vals[i] = op(i + 1);
    auto simd_in = simd<t>::load(vals);
    auto simd_out = rev_op(simd_in);

    for (int i = 0; i < simd_length; i++)
    {
        t error, tolerance;
        if (!similar(simd_out[i], ((t)i + 1), &error, &tolerance))
        {
#if DEBUG
            std::cout.precision(std::numeric_limits<t>::max_digits10 + 1);
            std::cout << name << " input=" << i + 1 << " input=" << i + 1
                      << " output=" << simd_out[i] << " error=" << error
                      << " tolerance=" << tolerance << std::endl;
#endif
            return true;
        }
    }

    return false;
}

#define SV_TEST_T(type, op)                                                    \
    do                                                                         \
    {                                                                          \
        const char *name = "sv::" #type "::" #op;                              \
        if (sv_test<type>(name, ([&](auto x) { return op; }),                  \
                          ([&](auto x) { return op; })))                       \
        {                                                                      \
            std::cout << name << " test FAILED" << std::endl;                  \
            error |= true;                                                     \
        }                                                                      \
    } while (0);

#define RV_TEST_T(type, op, rev_op)                                            \
    do                                                                         \
    {                                                                          \
        const char *name = "rv::" #type "::" #rev_op;                          \
        if (rv_test<type>(name, ([&](auto x) { return op; }),                  \
                          ([&](auto x) { return rev_op; })))                   \
        {                                                                      \
            std::cout << name << " test FAILED" << std::endl;                  \
            error |= true;                                                     \
        }                                                                      \
    } while (0);

#define SV_TEST(op) SV_TEST_T(Real, op);

#define RV_TEST(op, rev_op)                                                    \
    RV_TEST_T(Real, op, rev_op);                                               \
    RV_TEST_T(Real, rev_op, op);

int main()
{
    bool error = false;

    SV_TEST(exp2(x));
    SV_TEST(exp(x));
    SV_TEST(log(x));
    SV_TEST(log2(x));
    // SV_TEST(pow(x,(decltype(x))2));
    SV_TEST(pow(x, 2));
    SV_TEST(sqrt(x));
    SV_TEST(sin(x));
    SV_TEST(cos(x));
    SV_TEST(tan(x));
    SV_TEST(sinh(x));
    SV_TEST(cosh(x));
    SV_TEST(tanh(x));

    SV_TEST_T(Real, simd_min(x, (Real)0.5));
    SV_TEST_T(Real, simd_max(x, (Real)0.5));

    RV_TEST(exp(x), log(x));
    // RV_TEST(pow(x,(decltype(x))2),sqrt(x));
    RV_TEST(pow(x, 2), sqrt(x));

    if (!error)
        std::cout << "Simd functions unit test passed" << std::endl;
    else
        std::cout << "Simd functions unit test NOT passed" << std::endl;

    return error;
}
