/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MATTERCCZ4RHS_HPP_)
#error "This file should only be included through MatterCCZ4RHS.hpp"
#endif

#ifndef MATTERCCZ4RHS_IMPL_HPP_
#define MATTERCCZ4RHS_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class matter_t, class gauge_t, class deriv_t>
MatterCCZ4RHS<matter_t, gauge_t, deriv_t>::MatterCCZ4RHS(
    matter_t a_matter, CCZ4_params_t<typename gauge_t::params_t> a_params,
    double a_dx, double a_sigma, int a_formulation, double a_G_Newton)
    : CCZ4RHS<gauge_t, deriv_t>(a_params, a_dx, a_sigma, a_formulation,
                                0.0 /*No cosmological constant*/),
      my_matter(a_matter), m_G_Newton(a_G_Newton)
{
}

template <class matter_t, class gauge_t, class deriv_t>
template <class data_t>
void MatterCCZ4RHS<matter_t, gauge_t, deriv_t>::compute(
    Cell<data_t> current_cell) const
{
    // copy data from chombo gridpoint into local variables
    const auto matter_vars = current_cell.template load_vars<Vars>();
    const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = this->m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
        this->m_deriv.template advection<Vars>(current_cell, matter_vars.shift);

    // Call CCZ4 RHS - work out RHS without matter, no dissipation
    Vars<data_t> matter_rhs;
    this->rhs_equation(matter_rhs, matter_vars, d1, d2, advec);

    // add RHS matter terms from EM Tensor
    add_emtensor_rhs(matter_rhs, matter_vars, d1);

    // add evolution of matter fields themselves
    my_matter.add_matter_rhs(matter_rhs, matter_vars, d1, d2, advec);

    // Add dissipation to all terms
    this->m_deriv.add_dissipation(matter_rhs, current_cell, this->m_sigma);

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(matter_rhs);
}

// Function to add in EM Tensor matter terms to CCZ4 rhs
template <class matter_t, class gauge_t, class deriv_t>
template <class data_t>
void MatterCCZ4RHS<matter_t, gauge_t, deriv_t>::add_emtensor_rhs(
    Vars<data_t> &matter_rhs, const Vars<data_t> &matter_vars,
    const Vars<Tensor<1, data_t>> &d1) const
{
    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(matter_vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // Calculate elements of the decomposed stress energy tensor
    const auto emtensor =
        my_matter.compute_emtensor(matter_vars, d1, h_UU, chris.ULL);

    // Update RHS for K and Theta depending on formulation
    if (this->m_formulation == CCZ4RHS<>::USE_BSSN)
    {
        matter_rhs.K += 4.0 * M_PI * m_G_Newton * matter_vars.lapse *
                        (emtensor.S + emtensor.rho);
        matter_rhs.Theta += 0.0;
    }
    else
    {
        matter_rhs.K += 4.0 * M_PI * m_G_Newton * matter_vars.lapse *
                        (emtensor.S - 3 * emtensor.rho);
        matter_rhs.Theta +=
            -8.0 * M_PI * m_G_Newton * matter_vars.lapse * emtensor.rho;
    }

    // Update RHS for other variables
    Tensor<2, data_t> Sij_TF = emtensor.Sij;
    make_trace_free(Sij_TF, matter_vars.h, h_UU);

    FOR(i, j)
    {
        matter_rhs.A[i][j] += -8.0 * M_PI * m_G_Newton * matter_vars.chi *
                              matter_vars.lapse * Sij_TF[i][j];
    }

    FOR(i)
    {
        data_t matter_term_Gamma = 0.0;
        FOR(j)
        {
            matter_term_Gamma += -16.0 * M_PI * m_G_Newton * matter_vars.lapse *
                                 h_UU[i][j] * emtensor.Si[j];
        }

        matter_rhs.Gamma[i] += matter_term_Gamma;
        matter_rhs.B[i] += matter_term_Gamma;
    }
}

#endif /* MATTERCCZ4RHS_IMPL_HPP_ */
