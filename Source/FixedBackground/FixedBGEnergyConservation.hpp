/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGENERGYCONSERVATION_HPP_
#define FIXEDBGENERGYCONSERVATION_HPP_

#include "ADMFixedBGVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "CoordinateTransformations.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the momentum flux S_i with type matter_t and writes it to the
//! grid
template <class matter_t, class background_t> class FixedBGEnergyConservation
{
    // Use the variable definition in the matter class
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const matter_t m_matter;                        //!< The matter object
    const double m_dx;                              //!< The grid spacing
    const background_t m_background;                //!< The metric background
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center

  public:
    FixedBGEnergyConservation(matter_t a_matter, background_t a_background,
                              double a_dx,
                              std::array<double, CH_SPACEDIM> a_center)
        : m_matter(a_matter), m_deriv(a_dx), m_dx(a_dx),
          m_background(a_background), m_center(a_center)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and derivs
        const auto vars = current_cell.template load_vars<MatterVars>();
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);

        // get the metric vars from the background
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, coords);

        using namespace TensorAlgebra;
        using namespace CoordinateTransformations;
        //	const auto gamma = metric_vars.gamma;
        const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);
        const emtensor_t<data_t> emtensor = m_matter.compute_emtensor(
            vars, metric_vars, d1, gamma_UU, chris_phys.ULL);
        const data_t det_gamma = compute_determinant_sym(metric_vars.gamma);
        const data_t R = coords.get_radius();

        // the unit vector in the radial direction
        Tensor<1, data_t> si_L;
        si_L[0] = coords.x / R;
        si_L[1] = coords.y / R;
        si_L[2] = coords.z / R;

        data_t rhoEnergy = emtensor.rho * metric_vars.lapse;
        data_t fluxEnergy = 0.0;
        data_t sourceEnergy = 0.0;

        FOR1(i) { rhoEnergy += -emtensor.Si[i] * metric_vars.shift[i]; }
        rhoEnergy *= sqrt(det_gamma);

        FOR1(i)
        {
            fluxEnergy += -si_L[i] * metric_vars.shift[i] * emtensor.rho;

            FOR1(j)
            {
                fluxEnergy += metric_vars.lapse * si_L[i] * emtensor.Si[j] *
                              gamma_UU[i][j];
            }
        }
        fluxEnergy *= det_gamma;

        current_cell.store_vars(rhoEnergy, c_rhoEnergy);
        current_cell.store_vars(fluxEnergy, c_fluxEnergy);
    }
};

#endif /* FIXEDBGENERGYCONSERVATION_HPP_ */
