/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMBINATORICS_HPP_
#define COMBINATORICS_HPP_

// Chombo includes
#include "CH_assert.H"

// Chombo namespace
#include "UsingNamespace.H"

namespace Combinatorics
{
// Calculate factorials
inline double factorial(int n)
{
    double out = 1.0;
    for (int i = 1; i <= n; i++)
    {
        out *= i;
    }
    return out;
}

// Calculate combinatorics
inline double n_choose_r(int n, int r)
{
    CH_assert((r >= 0) && (n >= r)); // sense check values

    double out = factorial(n) / (factorial(r) * factorial(n - r));
    return out;
}
} // namespace Combinatorics

#endif /* COMBINATORICS_HPP_ */
