/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGMOMENTUMFLUX_HPP_
#define FIXEDBGMOMENTUMFLUX_HPP_

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

//! Calculates the density rho with type matter_t and writes it to the grid
template <class matter_t, class background_t> class FixedBGMomentumFlux
{
    // Use the variable definition in the matter class
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const matter_t m_matter;         //!< The matter object
    const double m_dx;               //!< The matter object
    const background_t m_background; //!< The matter object
    const std::array<double, CH_SPACEDIM> m_center;

  public:
    FixedBGMomentumFlux(matter_t a_matter, background_t a_background,
                        double a_dx, std::array<double, CH_SPACEDIM> a_center)
        : m_matter(a_matter), m_deriv(a_dx), m_dx(a_dx),
          m_background(a_background), m_center(a_center)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and calc 1st
        // derivs
        const auto vars = current_cell.template load_vars<MatterVars>();
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // get the metric vars
        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, current_cell);

        // first gather some useful geometric quantities
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);
        const emtensor_t<data_t> emtensor = m_matter.compute_emtensor(
            vars, metric_vars, d1, gamma_UU, chris_phys.ULL);
        const data_t det_gamma = compute_determinant_sym(metric_vars.gamma);

        // The covector in the radial direction N_i in cartesian coords
        // (normal to surfaces of constant r)
        Tensor<1, data_t> Ni_L;
        data_t R = coords.get_radius();
        Ni_L[0] = coords.x / R;
        Ni_L[1] = coords.y / R;
        Ni_L[2] = coords.z / R;
        // normalise this
        data_t mod_N2 = 0.0;
        FOR2(i, j) { mod_N2 += gamma_UU[i][j] * Ni_L[i] * Ni_L[j]; }
        FOR1(i) { Ni_L[i] = Ni_L[i] / sqrt(mod_N2); }

        // the area element of the sphere
        data_t rho2 =
            simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
        data_t r2sintheta = sqrt(rho2) * R;
        data_t det_Sigma = CoordinateTransformations::get_det_spherical_area(
            metric_vars.gamma, coords.x, coords.y, coords.z);

        // The integrand for the x-momentum flux out of a radial
        // shell at the current position
        data_t Mdot = 0;

        FOR1(i)
        {
            Mdot += -metric_vars.shift[i] * Ni_L[i] * emtensor.Si[0];
            FOR1(j)
            {
                Mdot += metric_vars.lapse * gamma_UU[i][j] *
                        emtensor.Sij[0][j] * Ni_L[i];
            }
        }

        // the r2sintheta is taken care of by the integration of the flux
        // so just need the dA relating to the metric
        Mdot *= sqrt(det_Sigma) / r2sintheta;

        // Now (minus) the x Momentum density with volume factor
        data_t xMom = -emtensor.Si[0] * sqrt(det_gamma);

        // How big is the source of i mom?
        Tensor<1, data_t> source;
        FOR1(i)
        {
            source[i] = -emtensor.rho * metric_vars.d1_lapse[i];

            FOR1(j)
            {
                source[i] += emtensor.Si[j] * metric_vars.d1_shift[j][i];
                FOR2(k, l)
                {
                    source[i] += metric_vars.lapse * gamma_UU[k][l] *
                                 emtensor.Sij[k][j] * chris_phys.ULL[j][l][i];
                }
            }
            source[i] = source[i] * sqrt(det_gamma);
        }

        // assign values of Stress integrand in the output box
        current_cell.store_vars(source[0], c_Source);
        current_cell.store_vars(Mdot, c_Stress);
        current_cell.store_vars(xMom, c_xMom);
    }
};

#endif /* FIXEDBGMOMENTUMFLUX_HPP_ */
