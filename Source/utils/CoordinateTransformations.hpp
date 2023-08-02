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
template <class data_t, class data2_t>
Tensor<2, data_t, 3> spherical_jacobian(const data_t x, const data2_t y,
                                        const data2_t z)
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
    Tensor<2, data_t, 3> jac;
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
template <class data_t, class data2_t>
Tensor<2, data_t, 3> inverse_spherical_jacobian(const data_t x, const data2_t y,
                                                const data2_t z)
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
    Tensor<2, data_t, 3> inv_jac;
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

// transform a vector according to a given jacobian
// jacobian can be "inverse jacobian" in opposite transformation
template <class data_t, class data2_t>
Tensor<1, data_t, 3> transform_vector(const Tensor<1, data_t, 3> &vec_U,
                                      const Tensor<2, data2_t, 3> &jacobian)
{
    Tensor<1, data_t, 3> transformed_U;
    FOR(i)
    {
        transformed_U[i] = 0.0;
        FOR(j) { transformed_U[i] += jacobian[i][j] * vec_U[j]; }
    }
    return transformed_U;
}

// transform a covector according to a given jacobian
// jacobian can be "inverse jacobian" in opposite transformation
template <class data_t, class data2_t>
Tensor<1, data_t, 3> transform_covector(const Tensor<1, data_t, 3> &vec_L,
                                        const Tensor<2, data2_t, 3> &jacobian)
{
    Tensor<1, data_t, 3> transformed_L;
    FOR(i)
    {
        transformed_L[i] = 0.0;
        FOR(j) { transformed_L[i] += vec_L[j] * jacobian[j][i]; }
    }
    return transformed_L;
}

// transform a tensor with two raised indices according to a given jacobian
// jacobian can be "inverse jacobian" in opposite transformation
template <class data_t, class data2_t>
Tensor<2, data_t, 3> transform_tensor_UU(const Tensor<2, data_t, 3> &tensor_UU,
                                         const Tensor<2, data2_t, 3> &jacobian)
{
    Tensor<2, data_t, 3> transformed_UU;
    FOR(i, j)
    {
        transformed_UU[i][j] = 0.;
        FOR(k, l)
        {
            transformed_UU[i][j] +=
                jacobian[i][k] * jacobian[j][l] * tensor_UU[k][l];
        }
    }
    return transformed_UU;
}

// transform a tensor with two lowered indices according to a given jacobian
// jacobian can be "inverse jacobian" in opposite transformation
template <class data_t, class data2_t>
Tensor<2, data_t, 3> transform_tensor_LL(const Tensor<2, data_t, 3> &tensor_LL,
                                         const Tensor<2, data2_t, 3> &jacobian)
{
    Tensor<2, data_t, 3> transformed_LL;
    FOR(i, j)
    {
        transformed_LL[i][j] = 0.;
        FOR(k, l)
        {
            transformed_LL[i][j] +=
                tensor_LL[k][l] * jacobian[k][i] * jacobian[l][j];
        }
    }
    return transformed_LL;
}

// Convert a Tensor (with two lower indices) in spherical coords to cartesian
// coords
template <class data_t, class data2_t>
Tensor<2, data_t, 3>
spherical_to_cartesian_LL(const Tensor<2, data_t, 3> &spherical_g,
                          const data_t x, const data2_t y, const data2_t z)
{
    return transform_tensor_LL(spherical_g, spherical_jacobian(x, y, z));
}

// Convert a Tensor (with two upper indices) in spherical coords to cartesian
// coords
template <class data_t, class data2_t>
Tensor<2, data_t, 3>
spherical_to_cartesian_UU(const Tensor<2, data_t, 3> &spherical_g_UU,
                          const data_t x, const data2_t y, const data2_t z)
{
    return transform_tensor_UU(spherical_g_UU,
                               inverse_spherical_jacobian(x, y, z));
}

// Convert a Tensor (with two lower indices) in cartesian coords to spherical
// coords
template <class data_t, class data2_t>
Tensor<2, data_t, 3>
cartesian_to_spherical_LL(const Tensor<2, data_t, 3> &cartesian_g,
                          const data_t x, const data2_t y, const data2_t z)
{
    return transform_tensor_LL(cartesian_g,
                               inverse_spherical_jacobian(x, y, z));
}

// Convert a Tensor (with two upper indices) in cartesian coords to spherical
// coords
template <class data_t, class data2_t>
Tensor<2, data_t, 3>
cartesian_to_spherical_UU(const Tensor<2, data_t, 3> &cartesian_g_UU, data_t x,
                          data2_t y, data2_t z)
{
    return transform_tensor_UU(cartesian_g_UU, spherical_jacobian(x, y, z));
}

// Convert a vector (with one upper index) in spherical coords to cartesian
// coords
template <class data_t, class data2_t>
Tensor<1, data_t, 3>
spherical_to_cartesian_U(const Tensor<1, data_t, 3> &spherical_v_U, data_t x,
                         data2_t y, data2_t z)
{
    return transform_vector(spherical_v_U, inverse_spherical_jacobian(x, y, z));
}

// Convert a vector (with one lower index) in spherical coords to cartesian
// coords
template <class data_t, class data2_t>
Tensor<1, data_t, 3>
spherical_to_cartesian_L(const Tensor<1, data_t, 3> &spherical_v_L, data_t x,
                         data2_t y, data2_t z)
{
    return transform_covector(spherical_v_L, spherical_jacobian(x, y, z));
}

// Convert a vector (with one upper index) in cartesian coords to spherical
// coords
template <class data_t, class data2_t>
Tensor<1, data_t, 3>
cartesian_to_spherical_U(const Tensor<1, data_t, 3> &cartesian_v_U, data_t x,
                         data2_t y, data2_t z)
{
    return transform_vector(cartesian_v_U, spherical_jacobian(x, y, z));
}

// Convert a vector (with one lower index) in cartesian coords to spherical
// coords
template <class data_t, class data2_t>
Tensor<1, data_t, 3>
cartesian_to_spherical_L(const Tensor<1, data_t, 3> &cartesian_v_L, data_t x,
                         data2_t y, data2_t z)
{
    return transform_covector(cartesian_v_L,
                              inverse_spherical_jacobian(x, y, z));
}

// The area element of a sphere
template <class data_t>
data_t area_element_sphere(const Tensor<2, data_t, 3> &spherical_g)
{
    return sqrt(spherical_g[1][1] * spherical_g[2][2] -
                spherical_g[1][2] * spherical_g[2][1]);
}

// cartesian coordinates rotation matrix
// axis must be normalized to 1
template <class data_t>
Tensor<2, data_t, 3> rotation_matrix(const Tensor<1, data_t, 3> &axis,
                                     data_t cos_angle)
{
    data_t sine = sqrt(1. - cos_angle * cos_angle);
    data_t one_minus_cos = 1. - cos_angle;

    // as in
    // https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    Tensor<2, data_t, 3> R;
    R[0][0] = axis[0] * axis[0] * one_minus_cos + cos_angle;
    R[0][1] = axis[1] * axis[0] * one_minus_cos - axis[2] * sine;
    R[0][2] = axis[2] * axis[0] * one_minus_cos + axis[1] * sine;
    R[1][0] = axis[0] * axis[1] * one_minus_cos + axis[2] * sine;
    R[1][1] = axis[1] * axis[1] * one_minus_cos + cos_angle;
    R[1][2] = axis[2] * axis[1] * one_minus_cos - axis[0] * sine;
    R[2][0] = axis[0] * axis[2] * one_minus_cos - axis[1] * sine;
    R[2][1] = axis[1] * axis[2] * one_minus_cos + axis[0] * sine;
    R[2][2] = axis[2] * axis[2] * one_minus_cos + cos_angle;

    return R;
}

// cartesian coordinates rotation matrix
template <class data_t>
Tensor<2, data_t, 3> rotation_matrix(const Tensor<1, data_t, 3> &origin,
                                     const Tensor<1, data_t, 3> &destination)
{
    // Calculate axis of rotation:
    Tensor<1, data_t, 3> axis = {
        origin[1] * destination[2] - origin[2] * destination[1],
        origin[2] * destination[0] - origin[0] * destination[2],
        origin[0] * destination[1] - origin[1] * destination[0]};

    static const double eps = 1.e-13;
    if (std::abs(axis[0]) < eps && std::abs(axis[1]) < eps &&
        std::abs(axis[2]) < eps)
        return rotation_matrix({0., 0., 0.}, 1.);

    data_t norm_inv =
        1.0 / sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);

    axis[0] *= norm_inv;
    axis[1] *= norm_inv;
    axis[2] *= norm_inv;

    // Calculate angle of rotation
    data_t norm_origin = sqrt(origin[0] * origin[0] + origin[1] * origin[1] +
                              origin[2] * origin[2]);
    data_t norm_dest =
        sqrt(destination[0] * destination[0] + destination[1] * destination[1] +
             destination[2] * destination[2]);

    data_t cos_angle = origin[0] * destination[0] + origin[1] * destination[1] +
                       origin[2] * destination[2];
    cos_angle /= (norm_origin * norm_dest);

    return rotation_matrix(axis, cos_angle);
}

} // namespace CoordinateTransformations
#endif /* COORDINATETRANSFORMATIONS_HPP_ */
