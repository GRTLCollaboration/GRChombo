/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef XSQUARED_HPP_
#define XSQUARED_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FixedBGProcaFieldTest.hpp"
#include "KerrSchildFixedBG.hpp"
#include "Potential.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates the initial conditions
class XSquared
{
  protected:
    const Potential::params_t m_potential_params;
    const KerrSchildFixedBG::params_t m_bg_params;
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    // The evolution vars
    template <class data_t>
    using Vars = FixedBGProcaFieldTest<Potential>::template Vars<data_t>;

  public:
    //! The constructor for the class
    XSquared(const Potential::params_t a_potential_params,
             const KerrSchildFixedBG::params_t a_bg_params,
             const std::array<double, CH_SPACEDIM> a_center, const double a_dx)
        : m_dx(a_dx), m_center(a_center),
          m_potential_params(a_potential_params), m_bg_params(a_bg_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // get the metric vars
        KerrSchildFixedBG kerr_bh(m_bg_params, m_dx);
        MetricVars<data_t> metric_vars;
        kerr_bh.compute_metric_background(metric_vars, current_cell);
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
        // The matter vars
        const auto vars = current_cell.template load_vars<Vars>();

        // Xsquared = X^/mu X_/mu
        data_t Xsquared;
        Xsquared = -vars.phi * vars.phi;
        FOR2(i, j) { Xsquared += gamma_UU[i][j] * vars.Avec[j] * vars.Avec[i]; }

        // Store the initial values of the variables
        current_cell.store_vars(Xsquared, c_Xsquared);
    }
};

#endif /* XSQUARED_HPP_ */
