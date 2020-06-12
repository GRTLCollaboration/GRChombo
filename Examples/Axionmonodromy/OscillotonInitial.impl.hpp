/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(OSCILLOTONINITIAL_HPP_)
#error "This file should only be included through OscillotonInitial.hpp"
#endif

#ifndef OSCILLOTONINITIAL_IMPL_HPP_
#define OSCILLOTONINITIAL_IMPL_HPP_

// Compute the value of the initial vars on the grid
void OscillotonInitial::compute(Cell<double> current_cell) const
{
    // Make a vars object
    auto vars =
        current_cell.template load_vars<MatterCCZ4<ScalarField<>>::Vars>();
    // Define Coordinates
    Coordinates<double> coords(current_cell, m_dx, m_center);

    // Interpolate data from read in values
    double rr = coords.get_radius();
    int indxL = static_cast<int>(floor(rr / m_spacing));
    int indxH = static_cast<int>(ceil(rr / m_spacing));

    vars.lapse = m_lapse_values[indxL] +
                  (rr / m_spacing - indxL) *
                      (m_lapse_values[indxH] - m_lapse_values[indxL]);
    vars.Pi =
        m_Pi_values[indxL] +
        (rr / m_spacing - indxL) * (m_Pi_values[indxH] - m_Pi_values[indxL]);

    double psi =
        m_psi_values[indxL] +
        (rr / m_spacing - indxL) * (m_psi_values[indxH] - m_psi_values[indxL]);

    // assign value of chi
    vars.Pi = - m_sign_of_Pi * vars.Pi;
    vars.chi = pow(psi, -4.0);
    FOR1(i) {vars.h[i][i]=1.0;}

    current_cell.store_vars(vars);
}

#endif /* OSCILLOTONINITIAL_IMPL_HPP_ */
