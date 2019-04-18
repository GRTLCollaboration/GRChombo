/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SCALARGAUSS_HPP_)
#error "This file should only be included through ScalarGauss.hpp"
#endif

#ifndef SCALARGAUSS_IMPL_HPP_
#define SCALARGAUSS_IMPL_HPP_

inline ScalarGauss::ScalarGauss(params_t a_params, double a_dx)
    : m_dx(a_dx), m_params(a_params)
{
}

// Compute the value of the initial vars on the grid
template <class data_t>
void ScalarGauss::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<ScalarField<>>::Vars<data_t> vars;
    VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx, m_params.centerSF);

    // set the field vars
    vars.phi = compute_phi(coords);
    vars.Pi = 0;

    // start with unit lapse and flat metric (must be relaxed for chi)
    vars.lapse = 1;
    vars.chi = 1;

    // conformal metric is flat
    FOR1(i) vars.h[i][i] = 1.;

    // set the K var
    vars.K = compute_K(coords);

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

// Compute the value of phi at the current point
template <class data_t>
data_t ScalarGauss::compute_phi(Coordinates<data_t> coords) const
{
    data_t rr = coords.get_radius();
    // data_t rr2 = rr * rr;
    data_t out_phi =  m_params.bkgSF +
		     m_params.amplitudeSF *
                     exp(-pow( (rr - m_params.r_zero) / m_params.widthSF, 2.0));


    return out_phi;
}


// Compute the value of phi at the current point
template <class data_t>
data_t ScalarGauss::compute_K(Coordinates<data_t> coords) const
{
    // data_t rr = coords.get_radius();
    // data_t rr2 = rr * rr;
    data_t out_K =  m_params.kfactor;


    return out_K;
}

#endif /* SCALARGAUSS_IMPL_HPP_ */
