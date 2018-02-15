/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CONSTRAINTS_HPP_)
#error "This file should only be included through Constraints.hpp"
#endif

#ifndef CONSTRAINTS_IMPL_HPP_
#define CONSTRAINTS_IMPL_HPP_

#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "VarsTools.hpp"

inline Constraints::Constraints(double dx,
                                double cosmological_constant /*defaulted*/)
    : m_deriv(dx), m_cosmological_constant(cosmological_constant)
{
}

template <class data_t>
void Constraints::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

    constraints_t<data_t> out = constraint_equations(vars, d1, d2);

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(out.Ham, c_Ham);
    current_cell.store_vars(out.Mom, GRInterval<c_Mom1, c_Mom3>());
}

template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
Constraints::constraints_t<data_t> Constraints::constraint_equations(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2) const
{
    constraints_t<data_t> out;

    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    auto ricci = CCZ4Geometry::compute_ricci(vars, d1, d2, h_UU, chris);

    auto A_UU = TensorAlgebra::raise_all(vars.A, h_UU);
    data_t tr_A2 = TensorAlgebra::compute_trace(vars.A, A_UU);

    out.Ham = ricci.scalar +
              (GR_SPACEDIM - 1.) * vars.K * vars.K / GR_SPACEDIM - tr_A2;
    out.Ham -= 2 * m_cosmological_constant;

    Tensor<2, data_t> covd_A[CH_SPACEDIM];
    FOR3(i, j, k)
    {
        covd_A[i][j][k] = d1.A[j][k][i];
        FOR1(l)
        {
            covd_A[i][j][k] += -chris.ULL[l][i][j] * vars.A[l][k] -
                               chris.ULL[l][i][k] * vars.A[l][j];
        }
    }

    FOR1(i) { out.Mom[i] = -(GR_SPACEDIM - 1.) * d1.K[i] / GR_SPACEDIM; }
    FOR3(i, j, k)
    {
        out.Mom[i] += h_UU[j][k] *
                      (covd_A[k][j][i] - GR_SPACEDIM * vars.A[i][j] *
                                             d1.chi[k] / (2 * chi_regularised));
    }

    return out;
}

#endif /* CONSTRAINTS_IMPL_HPP_ */
