/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SCALARBUBBLE_HPP_)
#error "This file should only be included through ScalarBubble.hpp"
#endif

#ifndef SCALARBUBBLE_IMPL_HPP_
#define SCALARBUBBLE_IMPL_HPP_

inline ScalarBubble::ScalarBubble(params_t a_params, double a_dx)
    : m_dx(a_dx), m_params(a_params)
{
}

// Compute the value of the initial vars on the grid
template <class data_t>
void ScalarBubble::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4RHS<ScalarField<>>::Vars<data_t> vars;
    VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx, m_params.centerSF);

    // set the field vars
    vars.phi = compute_phi(coords);
    vars.Pi = 0;

    // start with unit lapse and flat metric (must be relaxed for chi)
    vars.lapse = 1;
    vars.chi = 1;

    // conformal metric is flat
    FOR(i) vars.h[i][i] = 1.;

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

// Compute the value of phi at the current point
template <class data_t>
data_t ScalarBubble::compute_phi(Coordinates<data_t> coords) const
{
    data_t rr = coords.get_radius();
    data_t rr2 = rr * rr;
    data_t out_phi = m_params.amplitudeSF * rr2 *
                     exp(-pow((rr - m_params.r_zero) / m_params.widthSF, 2.0));

    return out_phi;
}

#endif /* SCALARBUBBLE_IMPL_HPP_ */
