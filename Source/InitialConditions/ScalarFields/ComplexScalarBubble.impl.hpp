/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(COMPLEXSCALARBUBBLE_HPP_)
#error "This file should only be included through ComplexScalarBubble.hpp"
#endif

#ifndef COMPLEXSCALARBUBBLE_IMPL_HPP_
#define COMPLEXSCALARBUBBLE_IMPL_HPP_

inline ComplexScalarBubble::ComplexScalarBubble(params_t a_params, double a_dx)
    : m_dx(a_dx), m_params(a_params)
{
}

// Compute the value of the initial vars on the grid
template <class data_t>
void ComplexScalarBubble::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    VarsTools::assign(vars, 0.); // Set all vars to zero
    Coordinates<data_t> coords(current_cell, m_dx, m_params.centerSF);

    // set the field vars - so real and im parts out of phase
    vars.phi_Re = compute_gaussian(coords);
    vars.phi_Im = 0;
    vars.Pi_Re = 0;
    vars.Pi_Im = compute_gaussian(coords);

    // start with unit lapse and flat metric (must be relaxed for chi)
    vars.lapse = 1;
    vars.chi = 1;

    // conformal metric is flat
    FOR1(i) vars.h[i][i] = 1.;

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

// Compute the value of the gaussian at the current point
template <class data_t>
data_t ComplexScalarBubble::compute_gaussian(Coordinates<data_t> coords) const
{
    data_t rr = coords.get_radius();
    data_t rr2 = rr * rr;
    data_t out = m_params.amplitudeSF * rr2 *
                 exp(-pow(rr - m_params.r_zero / m_params.widthSF, 2.0));

    return out;
}

#endif /* COMPLEXSCALARBUBBLE_IMPL_HPP_ */
