/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGENERGYANDANGULARMOMFLUX_HPP_
#define FIXEDBGENERGYANDANGULARMOMFLUX_HPP_

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
template <class matter_t, class background_t>
class FixedBGEnergyAndAngularMomFlux
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
    const double m_z_over_x;

  public:
    FixedBGEnergyAndAngularMomFlux(
        const matter_t a_matter, const background_t a_background,
        const double a_dx, const std::array<double, CH_SPACEDIM> a_center,
        const double a_z_over_x = 1.0)
        : m_matter(a_matter), m_deriv(a_dx), m_dx(a_dx),
          m_background(a_background), m_center(a_center), m_z_over_x(a_z_over_x)
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

        // The unit covector normal to the surfaces in cartesian coords
        // relative to the ellipsoid with axis ratio z_over_x
        Tensor<1, data_t> Ni_L;
        data_t R = coords.get_radius();
        Ni_L[0] = coords.x / R;
        Ni_L[1] = coords.y / R;
        Ni_L[2] = coords.z / R / m_z_over_x / m_z_over_x;
        // normalise this to 1 using full metric g_munu
        data_t mod_N2 = 0.0;
        FOR2(i, j)
        {
            mod_N2 += gamma_UU[i][j] * Ni_L[i] * Ni_L[j];
//                    -  metric_vars.shift[i] * metric_vars.shift[j] * Ni_L[i] *
//                          Ni_L[j] / metric_vars.lapse / metric_vars.lapse;

        }
        FOR1(i) { Ni_L[i] = Ni_L[i] / sqrt(mod_N2); }

        // the area element of the sphere - should be spheroid, but approx
        // ok this just corrects for the determinant of Sigma
        data_t rho2 =
            simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
        data_t r2sintheta = sqrt(rho2) * R;
        data_t det_Sigma = CoordinateTransformations::get_det_spherical_area(
            metric_vars.gamma, coords.x, coords.y, coords.z);

        // dxdphi to convert tensor components to spherical polar
        Tensor<1, data_t> dxdphi;
        dxdphi[0] = -coords.y;
        dxdphi[1] = coords.x;
        dxdphi[2] = 0;

        // The integrand for the energy flux out of a radial
        // shell at the current position
        data_t Edot = 0.0;
        FOR1(i)
        {
            Edot += metric_vars.lapse * Ni_L[i] * emtensor.rho *
                    metric_vars.shift[i];

            FOR1(j)
            {
                Edot +=
                    -Ni_L[i] * emtensor.Si[j] *
                    (metric_vars.shift[i] * metric_vars.shift[j] +
                     metric_vars.lapse * metric_vars.lapse * gamma_UU[i][j]);
                FOR1(k)
                {
                    Edot += Ni_L[i] * metric_vars.lapse * gamma_UU[i][j] *
                            metric_vars.shift[k] * emtensor.Sij[j][k];
                }
            }
        }

        // This factor of det_Sigma takes care of the surface element
        // The r2sintheta part is counted in the coordinate integration
        // so remove it here
        Edot *= sqrt(det_Sigma) / r2sintheta;

        // The integrand for the angular momentum flux out of a radial
        // shell at the current position
        data_t Jdot = 0;
        FOR2(i, j)
        {
            Jdot +=
                -Ni_L[i] * metric_vars.shift[i] * emtensor.Si[j] * dxdphi[j];
            FOR1(k)
            {
                Jdot += metric_vars.lapse * emtensor.Sij[i][j] * dxdphi[j] *
                        gamma_UU[i][k] * Ni_L[k];
            }
        }
        Jdot *= sqrt(det_Sigma) / r2sintheta;

        // assign values of fluxes in output box
        current_cell.store_vars(Edot, c_Edot);
        current_cell.store_vars(Jdot, c_Jdot);
    }
};

#endif /* FIXEDBGENERGYANDANGULARMOMFLUX_HPP_ */
