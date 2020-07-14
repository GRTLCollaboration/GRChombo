/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COORDINATETRANSFORMATIONS_HPP_
#define COORDINATETRANSFORMATIONS_HPP_

#include "AlwaysInline.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

namespace CoordinateTransformations
{
// The jacobian matrix for transformation of one forms from spherical to
// cartesian and vectors from cartesian to spherical - drdx drdy drdz, dthetadx
// etc
template <class data_t>
static Tensor<2, data_t> jacobian_drdx(const data_t x, const double y,
                                       const double z)
{
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

    return jac;
}

// The jacobian matrix for transformation of vectors from spherical to cartesian
// and one forms from cartesian to spherical - dxdr dxdtheta dxdphi, dydr etc
template <class data_t>
static Tensor<2, data_t> jacobian_dxdr(const data_t x, const double y,
                                       const double z)
{
    // calculate useful position quantities
    data_t rho2 = simd_max(x * x + y * y, 1e-12);
    data_t rho = sqrt(rho2);
    data_t r2 = simd_max(x * x + y * y + z * z, 1e-12);
    data_t r = sqrt(r2);

    // And the sines and cosines of phi and theta
    data_t sin_theta = rho / r;
    data_t cos_phi = x / rho;
    data_t sin_phi = y / rho;

    // derivatives for inverse jacobian matrix - drdx etc
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

    return inv_jac;
}

// Convert a Tensor (with two lower indices) in spherical coords to cartesian
// coords
template <class data_t>
static Tensor<2, data_t>
    spherical_to_cartesian_LL(Tensor<2, data_t> spherical_g, data_t x, double y,
                              double z)
{
    // The output - g in cartesian coords
    Tensor<2, data_t> cartesian_g;

    // derivatives for jacobian matrix - drdx etc
    Tensor<2, data_t> jac = jacobian_drdx(x, y, z);

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

// Convert a Tensor (with two lower indices) in cartesian coords to spherical
// coords
template <class data_t>
static Tensor<2, data_t>
    cartesian_to_spherical_LL(Tensor<2, data_t> cartesian_g, data_t x, double y,
                              double z)
{
    // The output - g in spherical coords
    Tensor<2, data_t> spherical_g;

    // derivatives for jacobian matrix - dxdr etc
    Tensor<2, data_t> inv_jac = jacobian_dxdr(x, y, z);

    // Convert the Tensor to spherical coords
    FOR2(i, j)
    {
        spherical_g[i][j] = 0;
        FOR2(k, m)
        {
            spherical_g[i][j] +=
                cartesian_g[k][m] * inv_jac[k][i] * inv_jac[m][j];
        }
    }
    return spherical_g;
}

// Convert a vector (with one upper index) in spherical coords to cartesian
// coords
template <class data_t>
Tensor<1, data_t> spherical_to_cartesian_U(Tensor<1, data_t> spherical_v_U,
                                           data_t x, double y, double z)
{
    // The output - v in cartesian coords
    Tensor<1, data_t> cartesian_v_U;

    // derivatives for jacobian matrix - dxdr etc
    Tensor<2, data_t> inv_jac = jacobian_dxdr(x, y, z);

    // transform the vector to cartesian coords
    FOR1(i)
    {
        cartesian_v[i] = 0.0;
        FOR1(j) { cartesian_v_U[i] += inv_jac[i][j] * spherical_v_U[j]; }
    }
    return cartesian_v_U;
}

// Convert a one form (with one lower index) in spherical coords to cartesian
// coords
template <class data_t>
Tensor<1, data_t> spherical_to_cartesian_L(Tensor<1, data_t> spherical_v_L,
                                           data_t x, double y, double z)
{
    // The output - w in cartesian coords
    Tensor<1, data_t> cartesian_v_L;

    // derivatives for jacobian matrix - drdx etc
    Tensor<2, data_t> jac = jacobian_drdx(x, y, z);

    // transform the vector to cartesian coords
    FOR1(i)
    {
        cartesian_v_L[i] = 0.0;
        FOR1(j) { cartesian_v_L[i] += jac[j][i] * spherical_v_L[j]; }
    }
    return cartesian_v_L;
}

// Get the determinant of the spherical metric at a point
template <class data_t>
data_t get_det_spherical_area(Tensor<2, data_t> cartesian_g, data_t x, double y,
                              double z)
{
    // The metric in spherical coords
    const Tensor<2, data_t> spherical_g =
        cartesian_to_spherical_LL(cartesian_g, x, y, z);

    // the r components should now be zero in this adapted basis, so
    const data_t det_Sigma =
        spherical_g[1][1] * spherical_g[2][2] - spherical_g[1][2] * spherical_g[2][1];

    return det_Sigma;
}
} // namespace CoordinateTransformations
#endif /* COORDINATETRANSFORMATIONS_HPP_ */