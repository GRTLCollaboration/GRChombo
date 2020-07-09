/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(FIXEDBGPROCACONSTRAINTTEST_HPP_)
#error
"This file should only be included through FixedBGProcaConstraintTest.hpp"
#endif

#ifndef FIXEDBGPROCACONSTRAINTTEST_IMPL_HPP_
#define FIXEDBGPROCACONSTRAINTTEST_IMPL_HPP_

    template <class potential_t, class background_t>
    FixedBGProcaConstraintTest<potential_t, background_t>::
        FixedBGProcaConstraintTest(background_t a_background, double dx,
                                   double a_vector_mass,
                                   double a_vector_damping,
                                   const potential_t potential)
    : m_background(a_background), m_deriv(dx), m_vector_mass(a_vector_mass),
      m_vector_damping(a_vector_damping), m_potential(potential)
{
}

template <class potential_t, class background_t>
template <class data_t>
void FixedBGProcaConstraintTest<potential_t, background_t>::compute(
    Cell<data_t> current_cell) const
{
    // get the metric vars
    MetricVars<data_t> metric_vars;
    m_background.compute_metric_background(metric_vars, current_cell);
    const auto vars = current_cell.template load_vars<MatterVars>();
    const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);

    data_t gauss_constraint = constraint_equations(vars, metric_vars, d1);

    // write the gauss constraint onto the grid at this gridpoint
    current_cell.store_vars(gauss_constraint, c_gauss);
}

template <class potential_t, class background_t>
template <class data_t, template <typename> class vars_t>
data_t
FixedBGProcaConstraintTest<potential_t, background_t>::constraint_equations(
    const vars_t<data_t> &vars, const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1, data_t>> &d1) const
{
    // calculate full spatial christoffel symbols
    using namespace TensorAlgebra;
    const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
    const auto chris_phys = compute_christoffel(metric_vars.d1_gamma, gamma_UU);

    // the output constraint
    data_t gauss_constraint;

    // compute potential
    data_t dVdA = 0;
    data_t dphidt = 0;
    m_potential.compute_potential(dVdA, dphidt, vars, d1, metric_vars);

    // this is the second part of eqn 27
    // ie dVdA = mu^2 ( 1 + 4 c4 A^k A_k - 12 c4 phi^2)
    gauss_constraint = dVdA * vars.phi;

    // Now add D_i E^i
    FOR1(i)
    {
        gauss_constraint += d1.Evec[i][i];
        FOR1(j) { gauss_constraint += chris_phys.ULL[i][i][j] * vars.Evec[j]; }
    }

    // this should be zero for constraints satisfied
    return gauss_constraint;
}

#endif /* FIXEDBGPROCACONSTRAINTTEST_IMPL_HPP_ */
