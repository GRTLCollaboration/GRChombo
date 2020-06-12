/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SPHERICALHORIZON_HPP_)
#error "This file should only be included through SphericalHorizon.hpp"
#endif

#ifndef SPHERICALHORIZON_IMPL_HPP_
#define SPHERICALHORIZON_IMPL_HPP_

template <class data_t>
void SphericalHorizon::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    Coordinates<data_t> coords(current_cell, m_dx, m_center);

    data_t out = calculate_expansion(vars, d1, coords);

    //DEBUG HACK:
    // set the potential values
    const auto matter_vars = current_cell.template load_vars<MatterVars>();
    data_t V_of_phi = 0.0;
    data_t dVdphi = 0.0;

    // compute potential
    m_potential.compute_potential(V_of_phi, dVdphi, matter_vars);

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(out, c_Omega);
    current_cell.store_vars(V_of_phi, c_VofPhi);
}

template <class data_t>
data_t
SphericalHorizon::calculate_expansion(const Vars<data_t> &vars,
                                      const Vars<Tensor<1, data_t>> &d1,
                                      const Coordinates<data_t> coords) const
{
    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    const data_t r = coords.get_radius();
    const data_t x = coords.x;
    const double y = coords.y;
    const double z = coords.z;

    // use matrix ID to find dh_UUdx, dM^-1 = -M^-1 dM M^-1
    Tensor<3, data_t> dhUUdx;
    FOR3(i, j, k)
    {
        dhUUdx[i][j][k] = 0.0;
        FOR2(l, m)
        {
            dhUUdx[i][j][k] += -h_UU[i][l] * h_UU[j][m] * d1.h[l][m][k];
        }
    }

    Tensor<2, data_t> K_tensor;
    FOR2(i, j)
    {
        K_tensor[i][j] =
            (vars.A[i][j] + 1. / 3. * vars.h[i][j] * vars.K) / vars.chi;
    }

    // Using Thornburgs method in arXiv gr-qc/0306056 in spherical symmetry

    Tensor<1, data_t> s;
    s[0] = x / r;
    s[1] = y / r;
    s[2] = z / r;

    Tensor<2, data_t> dsdx;
    dsdx[0][0] = y * y * pow(r, -3.0) + z * z * pow(r, -3.0);
    dsdx[1][1] = x * x * pow(r, -3.0) + z * z * pow(r, -3.0);
    dsdx[2][2] = x * x * pow(r, -3.0) + y * y * pow(r, -3.0);
    dsdx[0][1] = -s[0] * s[1] / r;
    dsdx[1][0] = dsdx[0][1];
    dsdx[1][2] = -s[1] * s[2] / r;
    dsdx[2][1] = dsdx[1][2];
    dsdx[0][2] = -s[0] * s[2] / r;
    dsdx[2][0] = dsdx[0][2];

    // elements of the Thornburg formula
    data_t A = 0;
    data_t C = 0;
    FOR4(i, j, k, l)
    {
        A += -vars.chi * vars.chi * h_UU[i][k] * s[k] * vars.h[j][l] * s[l] *
                 dsdx[j][i] -
             0.5 * vars.chi * h_UU[i][j] * s[j] * s[k] * s[l] *
                 (vars.chi * dhUUdx[k][l][i] + h_UU[k][l] * d1.chi[i]);

        C += vars.chi * vars.chi * h_UU[i][k] * h_UU[j][l] * K_tensor[k][l] *
             s[i] * s[j];
    }

    data_t B = 0;
    data_t D = 0;
    FOR2(i, j)
    {
        B += s[j] * (vars.chi * dhUUdx[i][j][i] + h_UU[i][j] * d1.chi[i]) +
             vars.chi * h_UU[i][j] * dsdx[j][i] -
             1.5 * d1.chi[i] * h_UU[i][j] * s[j];

        D += vars.chi * h_UU[i][j] * s[i] * s[j];
    }

    return (A / D + B) / sqrt(D) + C / D - vars.K;
}

#endif /* SPHERICALHORIZON_IMPL_HPP_ */
