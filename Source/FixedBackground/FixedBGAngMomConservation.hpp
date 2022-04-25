/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGANGMOMCONSERVATION_HPP_
#define FIXEDBGANGMOMCONSERVATION_HPP_

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
template <class matter_t, class background_t> class FixedBGAngMomConservation
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
    FixedBGAngMomConservation(matter_t a_matter, background_t a_background,
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
        Tensor<2, data_t> spherical_gamma = cartesian_to_spherical_LL(
            metric_vars.gamma, coords.x, coords.y, coords.z);
        data_t dArea = area_element_sphere(spherical_gamma);
        const data_t R = coords.get_radius();
        data_t rho2 =
            simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
        data_t r2sintheta = sqrt(rho2) * R;

        // the unit vector in the radial direction
        Tensor<1, data_t> si_L;
        si_L[0] = coords.x / R;
        si_L[1] = coords.y / R;
        si_L[2] = coords.z / R;

        // Normalise
        data_t si_norm = 0.0;
        FOR2(i, j) { si_norm += gamma_UU[i][j] * si_L[i] * si_L[j]; }

        FOR1(i) { si_L[i] = si_L[i] / sqrt(si_norm); }

        Tensor<1, data_t> dxdphi;
        dxdphi[0] = -coords.y;
        dxdphi[1] = coords.x;
        dxdphi[2] = 0;

        data_t rhoAngMom = 0.0;
        data_t fluxAngMom = 0.0;
        data_t sourceAngMom = -emtensor.rho * metric_vars.d1_lapse[0];

        FOR1(i) { rhoAngMom += -emtensor.Si[i] * dxdphi[i]; }
        rhoAngMom *= sqrt(det_gamma);

        FOR2(i, j)
        {
            fluxAngMom +=
                -si_L[i] * metric_vars.shift[i] * emtensor.Si[j] * dxdphi[j];
            FOR1(k)
            {
                fluxAngMom += metric_vars.lapse * emtensor.Sij[i][j] *
                              dxdphi[j] * gamma_UU[i][k] * si_L[k];
            }
        }

        // dArea is the integration surface element; Divide by r2sintheta,
        // as that's accounted for in the SprericalExtraction
        fluxAngMom *= dArea / r2sintheta;

        Tensor<1, data_t> sourceLM; // TEMP LIN MOM
        FOR1(i)
        {
            sourceLM[i] = -emtensor.rho * metric_vars.d1_lapse[i];

            FOR1(j)
            {
                sourceLM[i] += emtensor.Si[j] * metric_vars.d1_shift[j][i];
                FOR2(k, l)
                {
                    sourceLM[i] += metric_vars.lapse * gamma_UU[k][l] *
                                   emtensor.Sij[k][j] * chris_phys.ULL[j][l][i];
                }
            }
            sourceAngMom += sourceLM[i] * dxdphi[i];
        }

        sourceAngMom = sourceAngMom * sqrt(det_gamma);

        current_cell.store_vars(rhoAngMom, c_rhoAngMom);
        current_cell.store_vars(fluxAngMom, c_fluxAngMom);
        current_cell.store_vars(sourceAngMom, c_sourceAngMom);
    }
};

#endif /* FIXEDBGANGMOMCONSERVATION_HPP_ */
