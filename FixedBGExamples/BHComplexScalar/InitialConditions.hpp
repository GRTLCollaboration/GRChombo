/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALCONDITIONS_HPP_
#define INITIALCONDITIONS_HPP_

#include "ADMFixedBGVars.hpp"
#include "BoostedBHFixedBG.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates the initial conditions
class InitialConditions
{
  protected:
    const double m_dx;
    const double m_amplitude_re, m_amplitude_im;
    const double m_mu;
    const double m_velocity;
    const std::array<double, CH_SPACEDIM> m_center;
    const BoostedBHFixedBG::params_t m_bg_params;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  public:
    //! The constructor for the class
    InitialConditions(const double a_amplitude_re, const double a_amplitude_im,
                      const double a_mu, const double a_velocity,
                      const std::array<double, CH_SPACEDIM> a_center,
                      const BoostedBHFixedBG::params_t a_bg_params,
                      const double a_dx)
        : m_dx(a_dx), m_amplitude_re(a_amplitude_re),
          m_amplitude_im(a_amplitude_im), m_center(a_center), m_mu(a_mu),
          m_velocity(a_velocity), m_bg_params(a_bg_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // get the metric vars
        BoostedBHFixedBG boosted_bh(m_bg_params, m_dx);
        MetricVars<data_t> metric_vars;
        boosted_bh.compute_metric_background(metric_vars, current_cell);
        const data_t det_gamma =
            TensorAlgebra::compute_determinant_sym(metric_vars.gamma);

        const double t = 0.0;
        const data_t x = coords.x;
        const double k = m_mu * m_velocity;
        const double omega = sqrt(k*k + m_mu * m_mu);
        data_t phi_Re = m_amplitude_re * cos(omega * t + k * x);
        data_t phi_Im = m_amplitude_re * sin(omega * t + k * x);
        data_t Pi_Re = - m_amplitude_im * omega * sin(omega * t + k * x);
        data_t Pi_Im = m_amplitude_im * omega * cos(omega * t + k * x);

        // Store the initial values of the variables
        current_cell.store_vars(phi_Re, c_phi_Re);
        current_cell.store_vars(phi_Im, c_phi_Im);
        current_cell.store_vars(Pi_Re, c_Pi_Re);
        current_cell.store_vars(Pi_Im, c_Pi_Im);
    }
};

#endif /* INITIALCONDITIONS_HPP_ */
