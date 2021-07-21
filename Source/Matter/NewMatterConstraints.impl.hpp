/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(NEWMATTERCONSTRAINTS_HPP_)
#error "This file should only be included through NewMatterConstraints.hpp"
#endif

#ifndef NEWMATTERCONSTRAINTS_IMPL_HPP_
#define NEWMATTERCONSTRAINTS_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class matter_t>
MatterConstraints<matter_t>::MatterConstraints(
    const matter_t a_matter, double dx, double G_Newton, int a_c_Ham,
    const Interval &a_c_Moms, int a_c_Ham_abs_terms /* defaulted*/,
    const Interval &a_c_Moms_abs_terms /*defaulted*/)
    : Constraints(dx, a_c_Ham, a_c_Moms, a_c_Ham_abs_terms, a_c_Moms_abs_terms,
                  0.0 /*No cosmological constant*/),
      my_matter(a_matter), m_G_Newton(G_Newton)
{
}

template <class matter_t>
template <class data_t>
void MatterConstraints<matter_t>::compute(Cell<data_t> current_cell) const
{
    // Load local vars and calculate derivs
    const auto vars = current_cell.template load_vars<BSSNMatterVars>();
    const auto d1 = m_deriv.template diff1<BSSNMatterVars>(current_cell);
    const auto d2 = m_deriv.template diff2<BSSNMatterVars>(current_cell);

    // Inverse metric and Christoffel symbol
    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    // Get the non matter terms for the constraints
    Vars<data_t> out = constraint_equations(vars, d1, d2, h_UU, chris);

    // Energy Momentum Tensor
    const auto emtensor = my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    // Hamiltonian constraint
    if (m_c_Ham >= 0 || m_c_Ham_abs_terms >= 0)
    {
        out.Ham += -16.0 * M_PI * m_G_Newton * emtensor.rho;
        out.Ham_abs_terms += 16.0 * M_PI * m_G_Newton * abs(emtensor.rho);
    }

    // Momentum constraints
    if (m_c_Moms.size() > 0 || m_c_Moms_abs_terms.size() > 0)
    {
        FOR(i)
        {
            out.Mom[i] += -8.0 * M_PI * m_G_Newton * emtensor.Si[i];
            out.Mom_abs_terms[i] +=
                8.0 * M_PI * m_G_Newton * abs(emtensor.Si[i]);
        }
    }
    // Write the constraints into the output FArrayBox
    store_vars(out, current_cell);
}

#endif /* NEWMATTERCONSTRAINTS_IMPL_HPP_ */
