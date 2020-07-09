/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALCONDITIONS_HPP_
#define INITIALCONDITIONS_HPP_

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
class InitialConditions
{
  protected:
    const double m_dx;
    const double m_amplitude;
    const double m_mu;
    const std::array<double, CH_SPACEDIM> m_center;
    const KerrSchildFixedBG::params_t m_bg_params;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    // The evolution vars
    template <class data_t>
    using Vars = FixedBGProcaFieldTest<Potential>::template Vars<data_t>;

  public:
    //! The constructor for the class
    InitialConditions(const double a_amplitude, const double a_mu,
                      const std::array<double, CH_SPACEDIM> a_center,
                      const KerrSchildFixedBG::params_t a_bg_params,
                      const double a_dx)
        : m_dx(a_dx), m_amplitude(a_amplitude), m_center(a_center), m_mu(a_mu),
          m_bg_params(a_bg_params)
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
        const data_t det_gamma =
            TensorAlgebra::compute_determinant_sym(metric_vars.gamma);

        // useful coordinate quantities
        const double M = m_bg_params.mass;
        const double a = m_bg_params.spin;
        const double a2 = a * a;
        const double z = coords.z;
        const data_t rho = coords.get_radius();
        const data_t rho2 = rho * rho;

        // the Kerr Schild radius r
        const data_t r2 = 0.5 * (rho2 - a2) +
                          sqrt(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z);
        const data_t r = sqrt(r2);

        // the peak of the bound state for the n=l=0, j=m=1 state
        // See arXiv : 1704.05081 eqns 9 and 10
        double alpha = M * m_mu;
        double r_0 = 1.0 / (m_mu * alpha);

        // set the field variable to approx profile
        Vars<data_t> vars;
        VarsTools::assign(vars, 0.);
        data_t Avec = m_amplitude * (exp(-r / r_0)) / det_gamma;
        double coswt = 1.0;
        double sinwt = 0.0;

        // set the vector values
        vars.Avec[0] = -Avec * coswt;
        vars.Avec[1] = -Avec * sinwt;
        vars.Avec[2] = 0.0;
        /*
                // Evec approx E_i = - d_t A_i / lapse
                Tensor<1, data_t> Evec;
                Evec[0] = -m_mu * Avec * sinwt / metric_vars.lapse;
                Evec[1] = m_mu * Avec * coswt / metric_vars.lapse;
                Evec[2] = 0.0;

                // remember to raise the index for E^i
                const auto gamma_UU =
                    TensorAlgebra::compute_inverse_sym(metric_vars.gamma);
                FOR2(i, j) { vars.Evec[i] += Evec[j] * gamma_UU[i][j]; }
        */
        // Store the initial values of the variables
        current_cell.store_vars(vars);
    }
};

#endif /* INITIALCONDITIONS_HPP_ */
