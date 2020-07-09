/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(KERRBH_HPP_)
#error "This file should only be included through KerrBH.hpp"
#endif

#ifndef KERRBH_IMPL_HPP_
#define KERRBH_IMPL_HPP_

// Computes semi-isotropic Kerr solution as detailed in Liu, Etienne and Shapiro
// 2010, arxiv gr-qc/1001.4077
template <class data_t> void KerrBH::compute(Cell<data_t> current_cell) const
{
    // The cartesian variables and coords
    Vars<data_t> vars;
    VarsTools::assign(vars,
                      0.); // Set only the non-zero components explicitly below
    Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

    // Choose the desired Kerr solution
    switch (m_kerr_solution)
    {
    case Kerr::QUASI_ISOTROPIC:
        quasi_isotropic_kerr(vars, coords);

        break;
    case Kerr::KERR_SCHILD:
        kerr_schild(vars, coords);

        break;
    case Kerr::BOWEN_YORK:
        bowen_york(vars, coords);

        break;
    default:
        MayDay::Error("KerrBH::Supplied Kerr Solution not supported.");
    }

    // use a pre collapsed lapse, one, or analytic one
    switch (m_initial_lapse)
    {
    case Lapse::ONE:
        vars.lapse = 1.;
        break;
    case Lapse::PRE_COLLAPSED:
        vars.lapse = sqrt(vars.chi);
        break;
    case Lapse::ANALYTIC:
        vars.lapse = vars.lapse;
        break;
    default:
        MayDay::Error("KerrBH::Supplied initial lapse not supported.");
    }

    // Populate the variables on the grid
    // NB We stil need to set Gamma^i if it is NON ZERO
    // but we do this via a separate class/compute function
    // as we need the gradients of the metric which are not yet available
    current_cell.store_vars(vars);
}

// Quasi isotropic solution for high spin BHs per 1001.4077
template <class data_t>
void KerrBH::quasi_isotropic_kerr(ADMConformalVars::VarsWithGauge<data_t> &vars,
                                  const Coordinates<data_t> &coords) const
{
    // Kerr black hole params - mass M and spin a
    const double M = m_params.mass;
    const double a = m_params.spin;

    // work out where we are on the grid
    const data_t x = coords.x;
    const double y = coords.y;
    const double z = coords.z;

    // the radius, subject to a floor
    const data_t r = coords.get_radius();
    const data_t r2 = r * r;

    // the radius in xy plane, subject to a floor
    data_t rho2 = x * x + y * y;
    rho2 = simd_max(rho2, 1e-12);
    const data_t rho = sqrt(rho2);

    // calculate useful position quantities
    const data_t cos_theta = z / r;
    const data_t sin_theta = rho / r;
    const data_t cos_theta2 = cos_theta * cos_theta;
    const data_t sin_theta2 = sin_theta * sin_theta;

    // calculate useful metric quantities
    const double r_plus = M + sqrt(M * M - a * a);
    const double r_minus = M - sqrt(M * M - a * a);

    // The Boyer-Lindquist coordinate
    const data_t r_BL = r * pow(1.0 + 0.25 * r_plus / r, 2.0);

    // Other useful quantities per 1001.4077
    const data_t Sigma = r_BL * r_BL + a * a * cos_theta2;
    const data_t Delta = r_BL * r_BL - 2.0 * M * r_BL + a * a;
    // In the paper this is just 'A', but not to be confused with A_ij
    const data_t AA =
        pow(r_BL * r_BL + a * a, 2.0) - Delta * a * a * sin_theta2;
    // The rr component of the conformal spatial matric
    const data_t gamma_rr =
        Sigma * pow(r + 0.25 * r_plus, 2.0) / (r * r2 * (r_BL - r_minus));

    // Metric in semi isotropic Kerr-Schild coordinates
    // r, theta (t or th), phi (p)
    Tensor<2, data_t> spherical_g;
    Tensor<2, data_t> spherical_K;
    Tensor<1, data_t> spherical_shift;
    data_t kerr_lapse;

    // assign values
    FOR2(i, j) { spherical_g[i][j] = 0.0; }
    spherical_g[0][0] = gamma_rr;                // gamma_rr
    spherical_g[1][1] = Sigma;                   // gamma_tt
    spherical_g[2][2] = AA / Sigma * sin_theta2; // gamma_pp

    // Extrinsic curvature
    FOR2(i, j) { spherical_K[i][j] = 0.0; }

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

    using namespace InitialDataTools;
    // Convert spherical components to cartesian components using coordinate
    // transforms
    vars.h = spherical_to_cartesian_LL(spherical_g, x, y, z);
    vars.A = spherical_to_cartesian_LL(spherical_K, x, y, z);
    vars.shift = spherical_to_cartesian_U(spherical_shift, x, y, z);

    using namespace TensorAlgebra;
    // Convert to BSSN vars
    data_t deth = compute_determinant(vars.h);
    auto h_UU = compute_inverse_sym(vars.h);
    vars.chi = pow(deth, -1. / 3.);

    // transform extrinsic curvature into A and TrK - note h is still non
    // conformal version which is what we need here
    vars.K = compute_trace(vars.A, h_UU);
    make_trace_free(vars.A, vars.h, h_UU);

    // Make conformal
    FOR2(i, j)
    {
        vars.h[i][j] *= vars.chi;
        vars.A[i][j] *= vars.chi;
    }

    vars.lapse = kerr_lapse;
}

// Bowen York solution with puncture
// per http://www.livingreviews.org/Articles/Volume3/2000-5cook
// Approx solution to conformal factor given by gr-qc 9710096 Gleiser et al
// Accurate to order J^4, note that J is usually set assuming that M_adm = 1
// such that a = J/Madm = J. In this case m_bare \neq 1, and needs to be found
// using M_ADM = m_bare + 0.4 * J * J * m_bare ^ (-3)
template <class data_t>
void KerrBH::bowen_york(Vars<data_t> &vars,
                        const Coordinates<data_t> &coords) const
{
    // Kerr black hole params - mass M and spin a
    // For now J is in z direction only
    const double M = m_params.mass;
    const double J_z = m_params.spin;
    const std::array<double, CH_SPACEDIM> J = {0, 0, J_z};
    const std::array<double, CH_SPACEDIM> P = m_params.boost;

    // work out where we are on the grid
    const data_t x = coords.x;
    const double y = coords.y;
    const double z = coords.z;

    // the radius, subject to a floor
    const data_t r = coords.get_radius();
    const data_t r2 = r * r;
    const data_t r3 = r * r * r;
    const data_t cos_theta_squared = z * z / r2;

    const Tensor<1, data_t> n = {x / r, y / r, z / r};

    const Tensor<3, double> epsilon = TensorAlgebra::epsilon();

    // calculate \bar A_ij (the conformal one with \bar A_ij = psi^(2) A_ij)
    FOR2(i, j)
    {
        vars.A[i][j] = 1.5 / r2 * (P[i] * n[j] + P[j] * n[i]);
        FOR1(k)
        {
            vars.A[i][j] += -1.5 / r2 *
                            (TensorAlgebra::delta(i, j) - n[i] * n[j]) * P[k] *
                            n[k];

            FOR1(l)
            {
                vars.A[i][j] += 3.0 / r3 *
                                (epsilon[k][i][l] * J[l] * n[k] * n[j] +
                                 epsilon[k][j][l] * J[l] * n[k] * n[i]);
            }
        }
    }

    // This is an approx profile for chi - as above accurate to order J^4
    // See gr-qc 9710096 Gleiser et al, we use F for their a to not confuse with
    // spin
    data_t F = 0.5 * M;

    data_t psi_0 = 1 + F / r;

    data_t psi_02 = 0.025 * pow(F, -3.0) * pow(r + F, -5.0) *
                    (pow(F, 4.0) + pow(r, 4.0) + 5.0 * F * r * pow(r + F, 2.0));

    data_t psi_22 = -0.05 * pow(r + F, -5.0) * r * r / F;

    data_t psi =
        psi_0 +
        J_z * J_z * (psi_02 + psi_22 * 0.5 * (3.0 * cos_theta_squared - 1.0));

    vars.chi = pow(psi, -4.0);

    // Now calculate the conformal BSSN Aij
    FOR2(i, j)
    {
        vars.A[i][j] = pow(vars.chi, 1.5) * vars.A[i][j];
        vars.h[i][j] = 0.0;
    }

    // Set other vars - conformally flat metric
    FOR1(i) vars.h[i][i] = 1.0;
    vars.lapse = 1.0;
}

// Kerr Schild solution, must be used with excision at centre
// vars per https://arxiv.org/abs/gr-qc/0109032
template <class data_t>
void KerrBH::kerr_schild(ADMConformalVars::VarsWithGauge<data_t> &vars,
                         const Coordinates<data_t> &coords) const
{
    // Kerr black hole params - mass M and spin a
    const double M = m_params.mass;
    const double a = m_params.spin;

    // work out where we are on the grid
    const data_t x = coords.x;
    const double y = coords.y;
    const double z = coords.z;

    // the radius, subject to a floor
    const data_t r = coords.get_radius();
    const data_t r2 = r * r;

    // find r_BL the Kerr Schild radius (also called Boyer Lindquist)
    const data_t r_BL2 =
        0.5 * (r2 - a * a + sqrt(pow(a * a - r2, 2.0) + 4.0 * z * z * a * a));
    const data_t r_BL = sqrt(r_BL2);
    const data_t H = M * r_BL * r_BL2 / (r_BL2 * r_BL2 + a * a * z * z);
    const Tensor<1, data_t> el = {(r_BL * x + a * y) / (r_BL2 + a * a),
                                  (r_BL * y - a * x) / (r_BL2 + a * a),
                                  z / r_BL};

    // Calculate the gradients in el and H
    Tensor<1, data_t> dHdx;
    Tensor<2, data_t> dldx;
    get_KS_derivs(dHdx, dldx, r_BL, H, coords);

    // populate ADM vars per
    // http://www.livingreviews.org/Articles/Volume3/2000-5cook
    vars.lapse = pow(1.0 + 2.0 * H, -0.5);
    FOR2(i, j)
    {
        vars.h[i][j] = TensorAlgebra::delta(i, j) + 2 * H * el[i] * el[j];
    }
    FOR1(i) { vars.shift[i] = 2.0 * H * pow(vars.lapse, 2.0) * el[i]; }

    // test a = 0 case for now
    //    FOR2(i, j)
    //    {
    //        vars.A[i][j] = 2.0*H*vars.lapse/r*(TensorAlgebra::delta(i,j) -
    //                                                       (2 +
    //                                                       H)*el[i]*el[j]);
    //    }

    FOR2(i, j)
    {
        vars.A[i][j] = vars.lapse * (el[i] * dHdx[j] + el[j] * dHdx[i] +
                                     H * dldx[j][i] + H * dldx[i][j]);
        FOR1(k)
        {
            vars.A[i][j] += 2.0 * vars.lapse * H * el[k] *
                            (el[i] * el[j] * dHdx[k] + H * el[i] * dldx[j][k] +
                             H * el[j] * dldx[i][k]);
        }
    }

    using namespace TensorAlgebra;
    // Convert to BSSN vars
    auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    data_t deth = compute_determinant_sym(vars.h);
    vars.chi = pow(deth, -1. / 3.);

    // transform extrinsic curvature into A and TrK - note h is still non
    // conformal version which is what we need here
    vars.K = compute_trace(vars.A, h_UU);
    make_trace_free(vars.A, vars.h, h_UU);

    // Make conformal
    FOR2(i, j)
    {
        vars.h[i][j] *= vars.chi;
        vars.A[i][j] *= vars.chi;
    }
}

// Work out the gradients of the quantities H and el appearing in the Kerr
// Schild solution
template <class data_t>
void KerrBH::get_KS_derivs(Tensor<1, data_t> &dHdx, Tensor<2, data_t> &dldx,
                           const data_t &r_BL, const data_t &H,
                           const Coordinates<data_t> &coords) const
{
    // Kerr black hole params - spin a
    const double a = m_params.spin;

    // Find where we are
    data_t x = coords.x;
    double y = coords.y;
    double z = coords.z;
    data_t r = coords.get_radius();

    data_t drdx = x / r;
    data_t drdy = y / r;
    data_t drdz = z / r;

    data_t discriminant =
        pow(a, 4.0) - 2.0 * r * r * a * a + pow(r, 4.0) + 4.0 * z * z * a * a;

    data_t dr_BLdx =
        0.25 / r_BL *
        (2 * r * drdx + 0.5 / sqrt(discriminant) *
                            (-4.0 * a * a * r * drdx + 4.0 * r * r * r * drdx));
    data_t dr_BLdy =
        0.25 / r_BL *
        (2 * r * drdy + 0.5 / sqrt(discriminant) *
                            (-4.0 * a * a * r * drdy + 4.0 * r * r * r * drdy));
    data_t dr_BLdz =
        0.25 / r_BL *
        (2 * r * drdz + 0.5 / sqrt(discriminant) *
                            (-4.0 * a * a * r * drdz + 4.0 * r * r * r * drdz +
                             8.0 * z * a * a));

    dHdx[0] = H * (dr_BLdx * 3.0 / r_BL - 4.0 * pow(r_BL, 3.0) * dr_BLdx /
                                              (pow(r_BL, 4.0) + a * a * z * z));
    dHdx[1] = H * (dr_BLdy * 3.0 / r_BL - 4.0 * pow(r_BL, 3.0) * dr_BLdy /
                                              (pow(r_BL, 4.0) + a * a * z * z));
    dHdx[2] = H * (dr_BLdz * 3.0 / r_BL -
                   (4.0 * pow(r_BL, 3.0) * dr_BLdz + 2 * a * a * z) /
                       (pow(r_BL, 4.0) + a * a * z * z));

    // note to use convention as in rest of tensors the last index is the
    // derivative index so these are d_i l_x
    dldx[0][0] = (r_BL + x * dr_BLdx) * pow(r_BL * r_BL + a * a, -1.0) -
                 2 * r_BL * dr_BLdx * (r_BL * x + a * y) *
                     pow(r_BL * r_BL + a * a, -2.0);
    dldx[0][1] = (a + x * dr_BLdy) * pow(r_BL * r_BL + a * a, -1.0) -
                 2 * r_BL * dr_BLdy * (r_BL * x + a * y) *
                     pow(r_BL * r_BL + a * a, -2.0);
    dldx[0][2] = (x * dr_BLdz) * pow(r_BL * r_BL + a * a, -1.0) -
                 2 * r_BL * dr_BLdz * (r_BL * x + a * y) *
                     pow(r_BL * r_BL + a * a, -2.0);

    // these are d_i l_y
    dldx[1][0] = (y * dr_BLdx - a) * pow(r_BL * r_BL + a * a, -1.0) -
                 2 * r_BL * dr_BLdx * (r_BL * y - a * x) *
                     pow(r_BL * r_BL + a * a, -2.0);
    dldx[1][1] = (y * dr_BLdy + r_BL) * pow(r_BL * r_BL + a * a, -1.0) -
                 2 * r_BL * dr_BLdy * (r_BL * y - a * x) *
                     pow(r_BL * r_BL + a * a, -2.0);
    dldx[1][2] = (y * dr_BLdz) * pow(r_BL * r_BL + a * a, -1.0) -
                 2 * r_BL * dr_BLdz * (r_BL * y - a * x) *
                     pow(r_BL * r_BL + a * a, -2.0);

    // these are d_i l_z
    dldx[2][0] = -z * dr_BLdx / r_BL / r_BL;
    dldx[2][1] = -z * dr_BLdy / r_BL / r_BL;
    dldx[2][2] = 1.0 / r_BL - z * dr_BLdz / r_BL / r_BL;
}
#endif /* KERRBH_IMPL_HPP_ */
