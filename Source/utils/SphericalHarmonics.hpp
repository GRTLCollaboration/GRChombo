/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPHERICALHARMONICS_HPP_
#define SPHERICALHARMONICS_HPP_

#include "Combinatorics.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"

// Functions for the spin weighted spherical harmonics
// See paper arXiv:gr-qc/0610128 eqn 40
namespace SphericalHarmonics
{
template <class data_t> struct Y_lm_t
{
    data_t Real;
    data_t Im;
    data_t magnitude;
};

// Calculates the spin weight es, el, em spherical harmonic
template <class data_t>
Y_lm_t<data_t> spin_Y_lm(const data_t x, const double y, const double z,
                         const int es, const int el, const int em)
{

    CH_assert((el >= 0) && (el >= abs(em)));

    Y_lm_t<data_t> Y_lm;

    // calculate useful position quantities
    data_t r = simd_max(sqrt(x * x + y * y + z * z), 1e-6);
    data_t theta = acos(z / r);
    data_t phi = atan2(y, x);

    using namespace Combinatorics;
    double coefficient = pow(-1.0, es) * sqrt((2.0 * el + 1.0) / (4.0 * M_PI));
    coefficient *= sqrt(factorial(el + em) * factorial(el - em) /
                        factorial(el + es) / factorial(el - es));

    data_t sum = 0.0;
    int lower_limit = em + es > 0 ? em + es : 0;
    int upper_limit = el + em < el + es ? el + em : el + es;

    for (int i = lower_limit; i <= upper_limit; i++)
    {
        double temp = n_choose_r(el + es, i) * n_choose_r(el - es, i - es - em);
        sum += temp * pow(-1.0, i) *
               pow(cos(theta / 2.0), 2 * (el - i) + es + em) *
               pow(sin(theta / 2.0), 2 * i - em - es);
    }

    Y_lm.Real = coefficient * sum * cos(em * phi);
    Y_lm.Im = coefficient * sum * sin(em * phi);
    Y_lm.magnitude = sqrt(Y_lm.Real * Y_lm.Real + Y_lm.Im * Y_lm.Im);

    return Y_lm;
}
} // namespace SphericalHarmonics

#endif /* SPHERICALHARMONICS_HPP_ */
