/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(KERRBH_HPP_)
#error "This file should only be included through KerrBH.hpp"
#endif

#ifndef KERRBH_IMPL_HPP_
#define KERRBH_IMPL_HPP_

#include "DimensionDefinitions.hpp"

// Computes semi-isotropic Kerr solution as detailed in Liu, Etienne and Shapiro
// 2010, arxiv gr-qc/1001.4077
template <class data_t> void KerrBH::compute(Cell<data_t> current_cell) const
{
    using namespace CoordinateTransformations;
    using namespace TensorAlgebra;

    static const Tensor<1, double> z_dir = {0., 0., 1.};
    const Tensor<1, double> spin_dir = {m_params.spin_direction[0],
                                        m_params.spin_direction[1],
                                        m_params.spin_direction[2]};

    // define the rotation matrix needed to transform standard Cartesian (x,y,z)
    // coordinates into the coordinates of the spin direction
    Tensor<2, double> R = rotation_matrix(spin_dir, z_dir);

    // set up vars for the metric and extrinsic curvature, shift and lapse in
    // spherical coords
    Tensor<2, data_t> spherical_g;
    Tensor<2, data_t> spherical_K;
    Tensor<1, data_t> spherical_shift;
    data_t kerr_lapse;

    // The cartesian variables and coords
    Vars<data_t> vars;
    Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

    // rotate point
    Tensor<1, data_t> xyz = {coords.x, coords.y, coords.z};
    xyz = transform_vector(xyz, R);

    // Compute the components in spherical coords as per 1401.1548
    compute_kerr(spherical_g, spherical_K, spherical_shift, kerr_lapse, xyz);

    // work out where we are on the grid
    data_t x = xyz[0];
    data_t y = xyz[1];
    data_t z = xyz[2];

    // Convert spherical components to cartesian components using coordinate
    // transform_tensor_UU
    Tensor<2, data_t> cartesian_h =
        spherical_to_cartesian_LL(spherical_g, x, y, z);
    Tensor<2, data_t> cartesian_K =
        spherical_to_cartesian_LL(spherical_K, x, y, z);
    Tensor<1, data_t> cartesian_shift =
        spherical_to_cartesian_U(spherical_shift, x, y, z);

    // now rotate tensors back to original coordinates
    vars.h = transform_tensor_LL(cartesian_h, R);
    vars.A = transform_tensor_LL(cartesian_K, R);
    vars.shift = transform_vector(cartesian_shift, R);

    // Convert to BSSN vars
    data_t deth = compute_determinant(vars.h);
    auto h_UU = compute_inverse_sym(vars.h);
    vars.chi = pow(deth, -1. / 3.);

    // transform extrinsic curvature into A and TrK - note h is still non
    // conformal version which is what we need here
    vars.K = compute_trace(vars.A, h_UU);
    make_trace_free(vars.A, vars.h, h_UU);

    // Make conformal
    FOR(i, j)
    {
        vars.h[i][j] *= vars.chi;
        vars.A[i][j] *= vars.chi;
    }

    // use a pre collapsed lapse, could also use analytic one
    // vars.lapse = kerr_lapse;
    vars.lapse = pow(vars.chi, 0.5);

    // Populate the variables on the grid
    // NB We stil need to set Gamma^i which is NON ZERO
    // but we do this via a separate class/compute function
    // as we need the gradients of the metric which are not yet available
    current_cell.store_vars(vars);
}

template <class data_t>
void KerrBH::compute_kerr(Tensor<2, data_t> &spherical_g,
                          Tensor<2, data_t> &spherical_K,
                          Tensor<1, data_t> &spherical_shift,
                          data_t &kerr_lapse,
                          const Tensor<1, data_t> &coords) const
{
    // Kerr black hole params - mass M and spin a
    double M = m_params.mass;
    double a = m_params.spin;

    // work out where we are on the grid
    data_t x = coords[0];
    data_t y = coords[1];
    data_t z = coords[2];

    // the radius, subject to a floor
    data_t r = sqrt(D_TERM(x * x, +y * y, +z * z));
    static const double minimum_r = 1e-6;
    r = simd_max(r, minimum_r);

    data_t r2 = r * r;

    // the radius in xy plane, subject to a floor
    data_t rho2 = simd_max(x * x + y * y, 1e-12);
    data_t rho = sqrt(rho2);

    // calculate useful position quantities
    data_t cos_theta = z / r;
    data_t sin_theta = rho / r;
    data_t cos_theta2 = cos_theta * cos_theta;
    data_t sin_theta2 = sin_theta * sin_theta;

    // calculate useful metric quantities
    double r_plus = M + sqrt(M * M - a * a);
    double r_minus = M - sqrt(M * M - a * a);

    // The Boyer-Lindquist coordinate
    data_t r_BL = r * pow(1.0 + 0.25 * r_plus / r, 2.0);

    // Other useful quantities per 1001.4077
    data_t Sigma = r_BL * r_BL + a * a * cos_theta2;
    data_t Delta = r_BL * r_BL - 2.0 * M * r_BL + a * a;
    // In the paper this is just 'A', but not to be confused with A_ij
    data_t AA = pow(r_BL * r_BL + a * a, 2.0) - Delta * a * a * sin_theta2;
    // The rr component of the conformal spatial matric
    data_t gamma_rr =
        Sigma * pow(r + 0.25 * r_plus, 2.0) / (r * r2 * (r_BL - r_minus));

    // Metric in semi isotropic Kerr-Schild coordinates, r, theta (t or th), phi
    // (p)
    FOR(i, j) { spherical_g[i][j] = 0.0; }
    spherical_g[0][0] = gamma_rr;                // gamma_rr
    spherical_g[1][1] = Sigma;                   // gamma_tt
    spherical_g[2][2] = AA / Sigma * sin_theta2; // gamma_pp

    // Extrinsic curvature
    FOR(i, j) { spherical_K[i][j] = 0.0; }

    // set non zero elements of Krtp - K_rp, K_tp
    spherical_K[0][2] =
        a * M * sin_theta2 / (Sigma * sqrt(AA * Sigma)) *
        (3.0 * pow(r_BL, 4.0) + 2 * a * a * r_BL * r_BL - pow(a, 4.0) -
         a * a * (r_BL * r_BL - a * a) * sin_theta2) *
        (1.0 + 0.25 * r_plus / r) / sqrt(r * r_BL - r * r_minus);
    spherical_K[2][0] = spherical_K[0][2];
    spherical_K[2][1] = -2.0 * pow(a, 3.0) * M * r_BL * cos_theta * sin_theta *
                        sin_theta2 / (Sigma * sqrt(AA * Sigma)) *
                        (r - 0.25 * r_plus) * sqrt(r_BL / r - r_minus / r);
    spherical_K[1][2] = spherical_K[2][1];

    // set the analytic lapse
    kerr_lapse = sqrt(Delta * Sigma / AA);

    // set the shift (only the phi component is non zero)
    spherical_shift[0] = 0.0;
    spherical_shift[1] = 0.0;
    spherical_shift[2] = -2.0 * M * a * r_BL / AA;
}

#endif /* KERRBH_IMPL_HPP_ */
