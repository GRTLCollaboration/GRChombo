/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(NEWCONSTRAINTS_HPP_)
#error "This file should only be included through NewConstraints.hpp"
#endif

#ifndef NEWCONSTRAINTS_IMPL_HPP_
#define NEWCONSTRAINTS_IMPL_HPP_

#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "VarsTools.hpp"

inline Constraints::Constraints(
    double dx, int a_c_Ham, const Interval &a_c_Moms,
    int a_c_Ham_abs_terms /*defaulted*/,
    const Interval &a_c_Moms_abs_terms /*defaulted*/,
    double cosmological_constant /*defaulted*/)
    : m_deriv(dx), m_c_Ham(a_c_Ham), m_c_Moms(a_c_Moms),
      m_c_Ham_abs_terms(a_c_Ham_abs_terms),
      m_c_Moms_abs_terms(a_c_Moms_abs_terms),
      m_cosmological_constant(cosmological_constant)
{
}

template <class data_t>
void Constraints::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<MetricVars>();
    const auto d1 = m_deriv.template diff1<MetricVars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    Vars<data_t> out = constraint_equations(vars, d1, d2, h_UU, chris);

    store_vars(out, current_cell);
}

template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
Constraints::Vars<data_t> Constraints::constraint_equations(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const Tensor<2, data_t> &h_UU,
    const chris_t<data_t> &chris) const
{
    Vars<data_t> out;

    if (m_c_Ham >= 0 || m_c_Ham_abs_terms >= 0)
    {
        auto ricci = CCZ4Geometry::compute_ricci(vars, d1, d2, h_UU, chris);

        auto A_UU = TensorAlgebra::raise_all(vars.A, h_UU);
        data_t tr_A2 = TensorAlgebra::compute_trace(vars.A, A_UU);

        out.Ham = ricci.scalar +
                  (GR_SPACEDIM - 1.) * vars.K * vars.K / GR_SPACEDIM - tr_A2;
        out.Ham -= 2 * m_cosmological_constant;

        out.Ham_abs_terms =
            abs(ricci.scalar) + abs(tr_A2) +
            abs((GR_SPACEDIM - 1.) * vars.K * vars.K / GR_SPACEDIM);
        out.Ham_abs_terms += 2.0 * abs(m_cosmological_constant);
    }

    if (m_c_Moms.size() > 0 || m_c_Moms_abs_terms.size() > 0)
    {
        Tensor<2, data_t> covd_A[CH_SPACEDIM];
        FOR(i, j, k)
        {
            covd_A[i][j][k] = d1.A[j][k][i];
            FOR(l)
            {
                covd_A[i][j][k] += -chris.ULL[l][i][j] * vars.A[l][k] -
                                   chris.ULL[l][i][k] * vars.A[l][j];
            }
        }
        FOR(i)
        {
            out.Mom[i] = -(GR_SPACEDIM - 1.) * d1.K[i] / GR_SPACEDIM;
            out.Mom_abs_terms[i] = abs(out.Mom[i]);
        }
        Tensor<1, data_t> covd_A_term = 0.0;
        Tensor<1, data_t> d1_chi_term = 0.0;
        const data_t chi_regularised = simd_max(1e-6, vars.chi);
        FOR(i, j, k)
        {
            covd_A_term[i] += h_UU[j][k] * covd_A[k][j][i];
            d1_chi_term[i] += -GR_SPACEDIM * h_UU[j][k] * vars.A[i][j] *
                              d1.chi[k] / (2 * chi_regularised);
        }
        FOR(i)
        {
            out.Mom[i] += covd_A_term[i] + d1_chi_term[i];
            out.Mom_abs_terms[i] += abs(covd_A_term[i]) + abs(d1_chi_term[i]);
        }
    }
    return out;
}

template <class data_t>
void Constraints::store_vars(Vars<data_t> &out,
                             Cell<data_t> &current_cell) const
{
    if (m_c_Ham >= 0)
        current_cell.store_vars(out.Ham, m_c_Ham);
    if (m_c_Ham_abs_terms >= 0)
        current_cell.store_vars(out.Ham_abs_terms, m_c_Ham_abs_terms);
    if (m_c_Moms.size() == GR_SPACEDIM)
    {
        FOR(i)
        {
            int ivar = m_c_Moms.begin() + i;
            current_cell.store_vars(out.Mom[i], ivar);
        }
    }
    else if (m_c_Moms.size() == 1)
    {
        data_t Mom_sq = 0.0;
        FOR(i) { Mom_sq += out.Mom[i] * out.Mom[i]; }
        data_t Mom = sqrt(Mom_sq);
        current_cell.store_vars(Mom, m_c_Moms.begin());
    }
    if (m_c_Moms_abs_terms.size() == GR_SPACEDIM)
    {
        FOR(i)
        {
            int ivar = m_c_Moms_abs_terms.begin() + i;
            current_cell.store_vars(out.Mom_abs_terms[i], ivar);
        }
    }
    else if (m_c_Moms_abs_terms.size() == 1)
    {
        data_t Mom_abs_terms_sq = 0.0;
        FOR(i)
        {
            Mom_abs_terms_sq += out.Mom_abs_terms[i] * out.Mom_abs_terms[i];
        }
        data_t Mom_abs_terms = sqrt(Mom_abs_terms_sq);
        current_cell.store_vars(Mom_abs_terms, m_c_Moms_abs_terms.begin());
    }
}

#endif /* NEWCONSTRAINTS_IMPL_HPP_ */
