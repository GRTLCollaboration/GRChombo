/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(WARPFIELD_HPP_)
#error "This file should only be included through WarpField.hpp"
#endif

#ifndef WARPFIELD_IMPL_HPP_
#define WARPFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> WarpField::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;

    // energy density rho
    out.rho = vars.rho;
    // S_i
    FOR1(i) { out.Si[i] = vars.S_vector[i]; }
    // S_ij
    FOR2(i, j) { out.Sij[i][j] = vars.S_tensor[i][j]; }
    out.S = TensorAlgebra::compute_trace(out.Sij, h_UU);
    out.S *= vars.chi;

    return out;
}

// Adds in the RHS for the matter vars
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void WarpField::add_matter_rhs(rhs_vars_t<data_t> &total_rhs,
                               const vars_t<data_t> &vars,
                               const vars_t<Tensor<1, data_t>> &d1,
                               const diff2_vars_t<Tensor<2, data_t>> &d2,
                               const vars_t<data_t> &advec) const
{
    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // calculate full spatial christoffel symbols
    const Tensor<3, data_t> chris_phys = TensorAlgebra::compute_phys_chris(
        d1.chi, vars.chi, vars.h, h_UU, chris.ULL);

    // set the lapse and shift evolution to track values
    //    total_rhs.lapse = 0.1*(1.0 - vars.lapse);
    //    total_rhs.shift[0] = -0.05 * vars.shift[0];//- m_vs * d1.shift[0][0];
    //    total_rhs.shift[1] = 0.0;
    //    total_rhs.shift[2] = 0.0;

    //    This tracks the bubble
    //    total_rhs.rho = - m_vs * d1.rho[0];

    // set the density evolution using local energy conservation
    total_rhs.rho = vars.rho * vars.lapse * vars.K + advec.rho;
    FOR2(i, j)
    {
        total_rhs.rho += -vars.chi * h_UU[i][j] *
                         (vars.lapse * d1.S_vector[j][i] // NB this is d_i S_j
                          + 2.0 * vars.S_vector[j] * d1.lapse[i]);

        FOR1(k)
        {
            total_rhs.rho += vars.lapse * vars.chi * h_UU[i][j] *
                             vars.S_vector[k] * chris_phys[k][i][j];
            FOR1(l)
            {
                total_rhs.rho +=
                    vars.lapse * vars.chi * vars.S_tensor[i][j] * h_UU[i][k] *
                    h_UU[j][l] *
                    (vars.A[k][l] + 1.0 / 3.0 * vars.h[k][l] * vars.K);
            }
        }
    }

    // set evolution of S_i
    FOR1(i)
    {
        // this was for the moving bubble
        // total_rhs.S_vector[i] = - m_vs * d1.S_vector[i][0];

        total_rhs.S_vector[i] = advec.S_vector[i] +
                                vars.lapse * vars.K * vars.S_vector[i] -
                                vars.rho * d1.lapse[i];

        FOR1(j)
        {
            total_rhs.S_vector[i] += vars.S_vector[j] * d1.shift[j][i];

            FOR1(k)
            {
                total_rhs.S_vector[i] += -vars.chi * h_UU[j][k] *
                                         (vars.S_tensor[i][k] * d1.lapse[j] +
                                          vars.lapse * d1.S_tensor[i][k][j]);
                FOR1(l)
                {
                    total_rhs.S_vector[i] +=
                        vars.chi * h_UU[j][k] * vars.lapse *
                        (chris_phys[l][j][k] * vars.S_tensor[i][l] +
                         chris_phys[l][j][i] * vars.S_tensor[k][l]);
                }
            }
        }
    }

    // Set evolution of S_ij
    // For stationary would be - m_vs * d1.S_tensor[i][j][0];
    // Here we use the damping term plus an advection along the shift.
    FOR2(i, j)
    {
        total_rhs.S_tensor[i][j] =
            advec.S_tensor[i][j] -
            vars.lapse * (1.0 * (d1.S_vector[i][j] + d1.S_vector[j][i]) +
                          1.0 * (vars.S_vector[i] * d1.lapse[j] +
                                 vars.S_vector[j] * d1.lapse[i]) +
                          m_params.a1 * vars.K * vars.S_tensor[i][j] +
                          2.0 * vars.A[i][j] / vars.chi * vars.rho +
                          2.0 * vars.h[i][j] / vars.chi * vars.K * vars.rho);

        FOR1(k)
        {
            total_rhs.S_tensor[i][j] +=
                vars.S_tensor[k][j] * d1.shift[k][i] +
                vars.S_tensor[i][k] * d1.shift[k][j] +
                2.0 * vars.lapse * chris_phys[k][i][j] * vars.S_vector[k];
            FOR1(l)
            {
                total_rhs.S_tensor[i][j] += -m_params.a2 * vars.A[i][j] *
                                            vars.S_tensor[k][l] * h_UU[k][l] *
                                            vars.lapse;
                FOR2(m, n)
                {
                    total_rhs.S_tensor[i][j] +=
                        -m_params.a3 * vars.h[i][j] * vars.S_tensor[k][l] *
                        vars.A[m][n] * h_UU[k][m] * h_UU[l][n] * vars.chi *
                        vars.chi * vars.lapse;
                }
            }
        }
    }
}

#endif /* WARPFIELD_IMPL_HPP_ */
