/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PROCACONSTRAINT_HPP_
#define PROCACONSTRAINT_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FixedBGProcaField.hpp"
#include "KerrSchildFixedBG.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates the initial conditions
class ProcaConstraint
{
  protected:
    const double m_dx;
    const double m_mu;
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const std::array<double, CH_SPACEDIM> m_center;
    const KerrSchildFixedBG::params_t m_bg_params;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    // The evolution vars
    template <class data_t> using Vars = FixedBGProcaField::Vars<data_t>;

  public:
    //! The constructor for the class
    ProcaConstraint(const std::array<double, CH_SPACEDIM> a_center,
                    const KerrSchildFixedBG::params_t a_bg_params,
                    const double a_mu, const double a_dx)
        : m_center(a_center), m_bg_params(a_bg_params), m_mu(a_mu), m_dx(a_dx),
          m_deriv(a_dx)
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

        // calculate full spatial christoffel symbols
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);

        // The matter vars
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);

        // The constraint sets phi from the derivs of Evec
        data_t phi = 0.0;
        FOR1(i)
        {
            phi += -d1.Evec[i][i];
            FOR1(j) { phi += -chris_phys.ULL[i][i][j] * vars.Evec[j]; }
        }
        phi = phi / m_mu / m_mu;

        // Store the initial values of the variables
        current_cell.store_vars(phi, c_phi);
    }
};

#endif /* PROCACONSTRAINTS_HPP_ */
