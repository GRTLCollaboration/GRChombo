/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SCALARSIN3D_HPP_)
#error "This file should only be included through ScalarSin3d.hpp"
#endif

#ifndef SCALARSIN3D_IMPL_HPP_
#define SCALARSIN3D_IMPL_HPP_

inline ScalarSin3d::ScalarSin3d(params_t a_params, double a_dx)
    : m_dx(a_dx), m_params(a_params)
{
}

// Compute the value of the initial vars on the grid
template <class data_t>
void ScalarSin3d::compute(Cell<data_t> current_cell) const
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

    // set the K value
    vars.K = compute_K(coords);

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

// Compute the value of phi at the current point
template <class data_t>
data_t ScalarSin3d::compute_phi(Coordinates<data_t> coords) const
{
    // data_t rr = coords.get_radius();
    // data_t rr2 = rr * rr;
    data_t out_phi =  m_params.bkgSF +
                      + m_params.amplitudeSF * sin(2*M_PI * coords[0] * 1. / m_params.widthSF)
                      + m_params.amplitudeSF * sin(2*M_PI * coords[1] * 1. / m_params.widthSF)
                      + m_params.amplitudeSF * sin(2*M_PI * coords[2] * 1. / m_params.widthSF);


    return out_phi;
}


// Compute the value of phi at the current point
template <class data_t>
data_t ScalarSin3d::compute_K(Coordinates<data_t> coords) const
{
    //data_t rr = coords.get_radius();
    // data_t rr2 = rr * rr;
    data_t out_K =  m_params.kfactor; // + 0 * rr;


    return out_K;
}

#endif /* SCALARSIN3D_IMPL_HPP_ */
