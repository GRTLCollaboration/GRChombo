/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(ISOTROPICBOOSTEDBH_HPP_)
#error "This file should only be included through IsotropicBoostedBH.hpp"
#endif

#ifndef ISOTROPICBOOSTEDBH_IMPL_HPP_
#define ISOTROPICBOOSTEDBH_IMPL_HPP_

#if CH_SPACEDIM == 3
#include "CCZ4Vars.hpp"
#elif CH_SPACEDIM == 2
#include "CCZ4CartoonVars.hpp"
#endif

#include "Coordinates.hpp"
#include "Tensor.hpp"

// Compute the value of the initial vars on the grid
template <class data_t>
void IsotropicBoostedBH::compute(Cell<data_t> current_cell) const
{
    Coordinates<data_t> coords(current_cell, m_dx, m_bh_params.center);

    const double mass = m_bh_params.mass;
    const double vel = m_bh_params.momentum[0];
    const double gamma = 1. / sqrt(1. - vel * vel);

    // coords
    data_t x = coords.x * gamma;
    double y = coords.y;
#if CH_SPACEDIM == 3
    double z = coords.z;
#elif CH_SPACEDIM == 2
    double z = 0;
#endif
    data_t r = sqrt(x * x + y * y + z * z);

    // auxiliary vars
    data_t omega = 1. - mass / (2. * r);
    data_t psi = 1. + mass / (2. * r);

    data_t omega2 = omega * omega;
    data_t psi4 = psi * psi * psi * psi;
    data_t psi6 = psi4 * psi * psi;
    data_t dx_psi = -mass * x / (2. * pow(r, 3));
    data_t dy_psi = -mass * y / (2. * pow(r, 3));
    data_t dz_psi = -mass * z / (2. * pow(r, 3));

    data_t psi_omega = sqrt(psi6 - vel * vel * omega2);
    data_t beta_x = -vel * (psi6 - omega2) / (psi_omega * psi_omega);
    data_t g_xx = gamma * gamma * (psi_omega * psi_omega) / (psi * psi);
    data_t g_yy = psi4;
    data_t g_zz = psi4;

    // build tensors
    Tensor<2, data_t, 3> KLL = {0.};
    Tensor<2, data_t, 3> gammaLL = {0.};
    Tensor<2, data_t, 3> gammaUU = {0.};

    KLL[0][0] = -gamma * gamma * vel * dx_psi *
                (2. * psi6 * (psi + 2. * omega) - 2. * vel * vel * omega2) /
                (psi4 * psi * psi_omega);
    KLL[1][1] = 2. * vel * psi * omega * dx_psi / psi_omega;
    KLL[2][2] = KLL[1][1];
    KLL[0][1] = -gamma * vel * psi * dy_psi * (3 * omega + psi) / psi_omega;
    KLL[1][0] = KLL[0][1];
    KLL[2][0] = -gamma * vel * psi * dz_psi * (3 * omega + psi) / psi_omega;
    KLL[0][2] = KLL[2][0];
    KLL[2][1] = 0.;
    KLL[1][2] = KLL[1][2];

    gammaLL[0][0] = g_xx;
    gammaLL[1][1] = g_yy;
    gammaLL[2][2] = g_zz;
    gammaUU[0][0] = 1. / g_xx;
    gammaUU[1][1] = 1. / g_yy;
    gammaUU[2][2] = 1. / g_zz;

    data_t chi = pow(g_xx * g_yy * g_zz, -1. / 3.);

#if CH_SPACEDIM == 3
    CCZ4Vars::VarsWithGauge<data_t> vars;
#elif CH_SPACEDIM == 2
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
#endif

    // Set only the non-zero components explicitly below
    VarsTools::assign(vars, 0.);

    vars.chi = chi;
    vars.shift[0] = beta_x;
    vars.lapse = sqrt(vars.chi);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vars.K += KLL[i][j] * gammaUU[i][j];
        }
    }

    FOR2(i, j) { vars.h[i][j] = vars.chi * gammaLL[i][j]; }
    FOR2(i, j)
    {
        vars.A[i][j] = vars.chi * (KLL[i][j] - vars.K * gammaLL[i][j] / 3.);
    }

#if CH_SPACEDIM == 2
    vars.hww = vars.chi * gammaLL[2][2];
    vars.Aww = vars.chi * (KLL[2][2] - one_third * vars.K * gammaLL[2][2]);
#endif

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* ISOTROPICBOOSTEDBH _IMPL_HPP_ */
