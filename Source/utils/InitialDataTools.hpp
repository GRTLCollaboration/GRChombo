/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALDATATOOLS_HPP_
#define INITIALDATATOOLS_HPP_

#include "AlwaysInline.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

namespace InitialDataTools
{
// Convert a Tensor (with two lower indices) in spherical coords to cartesian
// coords
template <class data_t>
static Tensor<2, data_t>
    spherical_to_cartesian_LL(Tensor<2, data_t> spherical_g, data_t x, double y,
                              double z)
{
    Tensor<2, data_t> cartesian_g;

    // calculate useful position quantities
    data_t rho2 = simd_max(x * x + y * y, 1e-12);
    data_t rho = sqrt(rho2);
    data_t r2 = simd_max(x * x + y * y + z * z, 1e-12);
    data_t r = sqrt(r2);

    // And the sines and cosines of phi and theta
    data_t cos_phi = x / rho;
    data_t sin_phi = y / rho;

    // derivatives for jacobian matrix - drdx etc
    Tensor<2, data_t> jac;
    jac[0][0] = x / r;
    jac[1][0] = cos_phi * z / r2;
    jac[2][0] = -y / rho2;
    jac[0][1] = y / r;
    jac[1][1] = sin_phi * z / r2;
    jac[2][1] = x / rho2;
    jac[0][2] = z / r;
    jac[1][2] = -rho / r2;
    jac[2][2] = 0.0;

    // Convert the Tensor to cartesian coords
    FOR2(i, j)
    {
        cartesian_g[i][j] = 0;
        FOR2(k, m)
        {
            cartesian_g[i][j] += spherical_g[k][m] * jac[k][i] * jac[m][j];
        }
    }
    return cartesian_g;
}

// Convert a vector (with one upper index) in spherical coords to cartesian
// coords
template <class data_t>
Tensor<1, data_t> spherical_to_cartesian_U(Tensor<1, data_t> spherical_v,
                                           data_t x, double y, double z)
{
    Tensor<1, data_t> cartesian_v;

    // calculate useful position quantities
    data_t rho2 = simd_max(x * x + y * y, 1e-12);
    data_t rho = sqrt(rho2);
    data_t r2 = simd_max(x * x + y * y + z * z, 1e-12);
    data_t r = sqrt(r2);

    // And the sines and cosines of phi and theta
    data_t sin_theta = rho / r;
    data_t cos_phi = x / rho;
    data_t sin_phi = y / rho;

    // calculate the inverse jacobian, dxdr etc
    Tensor<2, data_t> inv_jac;
    inv_jac[0][0] = x / r;
    inv_jac[1][0] = y / r;
    inv_jac[2][0] = z / r;
    inv_jac[0][1] = z * cos_phi;
    inv_jac[1][1] = z * sin_phi;
    inv_jac[2][1] = -r * sin_theta;
    inv_jac[0][2] = -y;
    inv_jac[1][2] = x;
    inv_jac[2][2] = 0.0;

    // transform the vector to cartesian coords
    FOR1(i)
    {
        cartesian_v[i] = 0.0;
        FOR1(j) { cartesian_v[i] += inv_jac[i][j] * spherical_v[j]; }
    }
    return cartesian_v;
}
}
#endif /* INITIALDATATOOLS_HPP_ */
