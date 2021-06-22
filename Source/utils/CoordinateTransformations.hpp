/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COORDINATETRANSFORMATIONS_HPP_
#define COORDINATETRANSFORMATIONS_HPP_

#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"

namespace CoordinateTransformations
{

// Jacobian transformation matrix
template <class data_t>
static Tensor<2, data_t> spherical_jacobian(const data_t x, const double y,
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

// Incerse Jacobian
template <class data_t>
static Tensor<2, data_t>
inverse_spherical_jacobian(const data_t x, const double y, const double z)
{
    // calculate useful position quantities
    data_t rho2 = simd_max(x * x + y * y, 1e-12);
    data_t rho = sqrt(rho2);
    data_t r2 = simd_max(x * x + y * y + z * z, 1e-12);
    data_t r = sqrt(r2);

    // And the sines and cosines of phi and theta
    // data_t sin_theta = rho / r;
    data_t cos_phi = x / rho;
    data_t sin_phi = y / rho;

    // derivatives for inverse jacobian matrix - drdx etc
    Tensor<2, data_t> inv_jac;
    inv_jac[0][0] = x / r;
    inv_jac[1][0] = y / r;
    inv_jac[2][0] = z / r;
    inv_jac[0][1] = z * cos_phi;
    inv_jac[1][1] = z * sin_phi;
    inv_jac[2][1] = -rho;
    inv_jac[0][2] = -y;
    inv_jac[1][2] = x;
    inv_jac[2][2] = 0.0;
    return inv_jac;
}

// Convert a Tensor (with two lower indices) in spherical coords to cartesian
// coords
template <class data_t>
static Tensor<2, data_t>
spherical_to_cartesian_LL(const Tensor<2, data_t> &spherical_g, const data_t x,
                          const double y, const double z)
{
    Tensor<2, data_t> cartesian_g;

    // derivatives for jacobian matrix - drdx etc
    Tensor<2, data_t> jac = spherical_jacobian(x, y, z);

    // Convert the Tensor to cartesian coords
    FOR(i, j)
    {
        cartesian_g[i][j] = 0.;
        FOR(k, l)
        {
            cartesian_g[i][j] += spherical_g[k][l] * jac[k][i] * jac[l][j];
        }
    }
    return cartesian_g;
}

// Convert a Tensor (with two upper indices) in spherical coords to cartesian
// coords
template <class data_t>
static Tensor<2, data_t>
spherical_to_cartesian_UU(const Tensor<2, data_t> &spherical_g_UU,
                          const data_t x, const double y, const double z)
{
    Tensor<2, data_t> cartesian_g_UU;

    // derivatives for jacobian matrix - drdx etc
    Tensor<2, data_t> inv_jac = inverse_spherical_jacobian(x, y, z);

    // Convert the Tensor to cartesian coords
    FOR(i, j)
    {
        cartesian_g_UU[i][j] = 0.;
        FOR(k, l)
        {
            cartesian_g_UU[i][j] +=
                spherical_g_UU[k][l] * inv_jac[i][k] * inv_jac[j][l];
        }
    }
    return cartesian_g_UU;
}

// Convert a Tensor (with two lower indices) in cartesian coords to spherical
// coords
template <class data_t>
static Tensor<2, data_t>
cartesian_to_spherical_LL(const Tensor<2, data_t> &cartesian_g, const data_t x,
                          const double y, const double z)
{
    Tensor<2, data_t> spherical_g;

    // derivatives for inverse jacobian matrix - drdx etc
    Tensor<2, data_t> inv_jac = inverse_spherical_jacobian(x, y, z);

    // Convert the Tensor to spherical coords
    FOR(i, j)
    {
        spherical_g[i][j] = 0.;
        FOR(k, l)
        {
            spherical_g[i][j] +=
                cartesian_g[k][l] * inv_jac[k][i] * inv_jac[l][j];
        }
    }
    return spherical_g;
}

// Convert a Tensor (with two upper indices) in cartesian coords to spherical
// coords
template <class data_t>
static Tensor<2, data_t>
cartesian_to_spherical_UU(const Tensor<2, data_t> &cartesian_g_UU, data_t x,
                          double y, double z)
{
    Tensor<2, data_t> spherical_g_UU;

    // derivatives for jacobian matrix - drdx etc
    Tensor<2, data_t> jac = spherical_jacobian(x, y, z);

    // Convert the Tensor to spherical coords
    FOR(i, j)
    {
        spherical_g_UU[i][j] = 0.;
        FOR(k, l)
        {
            spherical_g_UU[i][j] +=
                cartesian_g_UU[k][l] * jac[i][k] * jac[j][l];
        }
    }
    return spherical_g_UU;
}

// Convert a vector (with one upper index) in spherical coords to cartesian
// coords
template <class data_t>
Tensor<1, data_t>
spherical_to_cartesian_U(const Tensor<1, data_t> &spherical_v_U, data_t x,
                         double y, double z)
{
    Tensor<1, data_t> cartesian_v_U;

    // derivatives for inverse jacobian matrix - drdx etc
    Tensor<2, data_t> inv_jac = inverse_spherical_jacobian(x, y, z);

    // transform the vector to cartesian coords
    FOR(i)
    {
        cartesian_v_U[i] = 0.0;
        FOR(j) { cartesian_v_U[i] += inv_jac[i][j] * spherical_v_U[j]; }
    }
    return cartesian_v_U;
}

// Convert a vector (with one lower index) in spherical coords to cartesian
// coords
template <class data_t>
Tensor<1, data_t>
spherical_to_cartesian_L(const Tensor<1, data_t> &spherical_v_L, data_t x,
                         double y, double z)
{
    Tensor<1, data_t> cartesian_v_L;

    // derivatives for jacobian matrix - drdx etc
    Tensor<2, data_t> jac = spherical_jacobian(x, y, z);

    // transform the vector to cartesian coords
    FOR(i)
    {
        cartesian_v_L[i] = 0.0;
        FOR(j) { cartesian_v_L[i] += spherical_v_L[j] * jac[j][i]; }
    }
    return cartesian_v_L;
}

// Convert a vector (with one upper index) in cartesian coords to spherical
// coords
template <class data_t>
Tensor<1, data_t>
cartesian_to_spherical_U(const Tensor<1, data_t> &cartesian_v_U, data_t x,
                         double y, double z)
{
    Tensor<1, data_t> spherical_v_U;

    // derivatives for jacobian matrix - drdx etc
    Tensor<2, data_t> jac = spherical_jacobian(x, y, z);

    // transform the vector to cartesian coords
    FOR(i)
    {
        spherical_v_U[i] = 0.0;
        FOR(j) { spherical_v_U[i] += jac[i][j] * cartesian_v_U[j]; }
    }
    return spherical_v_U;
}

// Convert a vector (with one lower index) in cartesian coords to spherical
// coords
template <class data_t>
Tensor<1, data_t>
cartesian_to_spherical_L(const Tensor<1, data_t> &cartesian_v_L, data_t x,
                         double y, double z)
{
    Tensor<1, data_t> spherical_v_L;

    // derivatives for inverse jacobian matrix - drdx etc
    Tensor<2, data_t> inv_jac = inverse_spherical_jacobian(x, y, z);

    // transform the vector to cartesian coords
    FOR(i)
    {
        spherical_v_L[i] = 0.0;
        FOR(j) { spherical_v_L[i] += cartesian_v_L[j] * inv_jac[j][i]; }
    }
    return spherical_v_L;
}

// The area element of a sphere
template <class data_t>
data_t area_element_sphere(const Tensor<2, data_t> &spherical_g)
{
    return sqrt(spherical_g[1][1] * spherical_g[2][2] -
                spherical_g[1][2] * spherical_g[2][1]);
}

} // namespace CoordinateTransformations
#endif /* COORDINATETRANSFORMATIONS_HPP_ */
