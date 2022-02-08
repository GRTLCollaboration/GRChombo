/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SCALARCONSTANT_HPP_
#define SCALARCONSTANT_HPP_

#include "ADMFixedBGVars.hpp"
#include "BoostedBHFixedBG.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FixedBGComplexScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates a constant scalar field given params for initial
//! matter config
class ScalarConstant
{
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  protected:
    const double m_dx;
    const double m_amplitude;
    const double m_mu;
    const std::array<double, CH_SPACEDIM> m_center;
    const BoostedBHFixedBG::params_t m_bg_params;

    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  public:
    //! The constructor for the class
    ScalarConstant(const double a_amplitude, const double a_mu,
                   const std::array<double, CH_SPACEDIM> a_center,
                   const BoostedBHFixedBG::params_t a_bg_params,
                   const double a_dx)
        : m_amplitude(a_amplitude), m_mu(a_mu), m_center(a_center),
          m_bg_params(a_bg_params), m_dx(a_dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

        FixedBGComplexScalarField<>::Vars<data_t> vars;

        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // get the metric vars
        BoostedBHFixedBG boosted_bh(m_bg_params, m_dx);
        MetricVars<data_t> metric_vars;
        boosted_bh.compute_metric_background(metric_vars, current_cell);

        // set the field vars
        vars.phi_Re = m_amplitude;
        vars.phi_Im = 0;
        vars.Pi_Re = 0;
        vars.Pi_Im = m_amplitude * m_mu; // metric_vars.lapse;

        current_cell.store_vars(vars);
    }
};

#endif /* SCALARCONSTANT_HPP_ */
