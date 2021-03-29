/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGCONSERVEDQUANTITIES_HPP_
#define FIXEDBGCONSERVEDQUANTITIES_HPP_

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
template <class matter_t, class background_t> class FixedBGConservedQuantities
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
    FixedBGConservedQuantities(matter_t a_matter, background_t a_background,
                               double a_dx,
                               std::array<double, CH_SPACEDIM> a_center)
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
        Tensor<1, data_t> Ni_L, si_L;
        data_t R = coords.get_radius();
        si_L[0] = coords.x / R;
        si_L[1] = coords.y / R;
        si_L[2] = coords.z / R;
        // normalise this
        data_t mod_N2 = 0.0;
        FOR2(i, j) { mod_N2 += gamma_UU[i][j] * si_L[i] * si_L[j]; }
        FOR1(i) { Ni_L[i] = si_L[i] / sqrt(mod_N2); }

        // the area element of the sphere
        // DONT NEED THIS IF DONT NORMALISE N
        data_t rxy2 = coords.x * coords.x + coords.y * coords.y;
        data_t rxy = sqrt(rxy2);
        rxy = simd_max(rxy, 1e-6);
        data_t r2sintheta = rxy * R;
        using namespace CoordinateTransformations;
        Tensor<2, data_t> spherical_gamma = cartesian_to_spherical_LL(
            metric_vars.gamma, coords.x, coords.y, coords.z);
        data_t sqrt_det_Sigma = area_element_sphere(spherical_gamma);

        // KV METHODS
        // ------------------------------

        // MOMENTUM QUANTITIES
        // KV METHOD
        // x Momentum density with volume factor
        data_t rhoM = emtensor.Si[0] * sqrt(det_gamma);

        // The integrand for the x-momentum flux
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
        Mdot *= sqrt_det_Sigma / r2sintheta;

        // How big is the sourceM of i mom?
        Tensor<1, data_t> sourceM;
        Tensor<1, data_t> BHMom; // TEMP CHECK ON BH MOM
        FOR1(i)
        {
            sourceM[i] = -emtensor.rho * metric_vars.d1_lapse[i];
            BHMom[i] = 0.0;

            FOR1(j)
            {
                sourceM[i] += emtensor.Si[j] * metric_vars.d1_shift[j][i];
                FOR1(k)
                {
                    BHMom[i] += -1.0 / 8.0 / M_PI * Ni_L[k] * gamma_UU[j][k] *
                                delta(i, j) * metric_vars.K;

                    FOR1(l)
                    {
                        sourceM[i] += metric_vars.lapse * gamma_UU[k][l] *
                                      emtensor.Sij[k][j] *
                                      chris_phys.ULL[j][l][i];
                        BHMom[i] += 1.0 / 8.0 / M_PI * Ni_L[j] *
                                    gamma_UU[k][i] * gamma_UU[l][j] *
                                    metric_vars.K_tensor[k][l];
                    }
                }
            }
            sourceM[i] = sourceM[i] * sqrt(det_gamma);
            BHMom[i] *= sqrt_det_Sigma / r2sintheta;
        }
        /*
                // ENERGY QUANTITIES
                // KV METHOD
                // Energy density
                data_t rhoE = emtensor.rho * metric_vars.lapse;
                FOR1(i) { rhoE += - emtensor.Si[i] * metric_vars.shift[i]; }
                rhoE *= sqrt(det_gamma);

                // Energy flux
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
                             metric_vars.lapse * metric_vars.lapse *
           gamma_UU[i][j]);
                        FOR1(k)
                        {
                            Edot += Ni_L[i] * metric_vars.lapse * gamma_UU[i][j]
           * metric_vars.shift[k] * emtensor.Sij[j][k];
                        }
                    }
                }
                Edot *= sqrt_det_Sigma / r2sintheta;

                // calculate the E source (should be zero?)
                data_t dlapsedt = 0.0;
                Tensor<1, data_t> dbetadt;
                FOR1(i) { dbetadt[i] = 0.0; }
                data_t sourceE = -emtensor.rho * dlapsedt;
                FOR1(i)
                {
                    sourceE += emtensor.Si[i] * dbetadt[i];
                    FOR1(j)
                    {
                        sourceE += emtensor.Si[i] *
                                   (-metric_vars.shift[i] / metric_vars.lapse *
                                        metric_vars.shift[j] *
           metric_vars.d1_lapse[j] + metric_vars.shift[j] *
           metric_vars.d1_shift[i][j]); FOR1(k)
                        {
                            sourceE +=
                                emtensor.Si[i] * metric_vars.shift[j] *
                                    (-metric_vars.lapse * gamma_UU[i][k] *
                                         metric_vars.K_tensor[k][j] +
                                     metric_vars.shift[i] * metric_vars.shift[k]
           * metric_vars.K_tensor[j][k] / metric_vars.lapse +
                                     metric_vars.shift[k] *
           chris_phys.ULL[i][j][k]) + emtensor.Sij[k][j] * gamma_UU[i][k] *
                                    metric_vars.lapse *
           metric_vars.d1_shift[j][i];
                            FOR1(l)
                            {
                                sourceE +=
                                    -metric_vars.lapse * metric_vars.lapse *
                                        gamma_UU[i][k] * gamma_UU[j][l] *
                                        emtensor.Sij[k][l] *
           metric_vars.K_tensor[i][j] + metric_vars.lapse * gamma_UU[i][k] *
                                        emtensor.Sij[k][j] *
           metric_vars.shift[l] * chris_phys.ULL[j][l][i];
                            }
                        }
                    }
                }
                sourceE *= sqrt(det_gamma);
        */
        // ANGULAR MOM J QUANTITIES
        // KV METHOD
        // Change basis
        Tensor<1, data_t> dxdphi;
        dxdphi[0] = -coords.y;
        dxdphi[1] = coords.x;
        dxdphi[2] = 0;

        // AM density
        data_t rhoJ = 0;
        FOR1(i) { rhoJ += -emtensor.Si[i] * dxdphi[i]; }
        rhoJ *= sqrt(det_gamma);

        // AM flux
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
        Jdot *= sqrt_det_Sigma / r2sintheta;

        // ADD J SOURCE??

        // LL METHODS
        // ------------------------------
        // ENERGY QUANTITIES
        // LL METHOD
        // Energy density
        data_t rhoE = emtensor.rho * det_gamma;

        // Energy flux
        data_t Edot = 0.0;
        FOR1(i)
        {
            Edot += -si_L[i] * metric_vars.shift[i] * emtensor.rho;

            FOR1(j)
            {
                Edot += metric_vars.lapse * si_L[i] * emtensor.Si[j] *
                        gamma_UU[i][j];
            }
        }
        Edot *= det_gamma;

        // calculate the E source (should be zero?)
        data_t dlapsedt = 0.0;
        data_t sourceE = emtensor.rho * dlapsedt / metric_vars.lapse;
        FOR1(i)
        {
            sourceE += -metric_vars.shift[i] * emtensor.rho *
                       metric_vars.d1_lapse[i] / metric_vars.lapse;
            FOR1(j)
            {
                sourceE += 2.0 * gamma_UU[i][j] * emtensor.Si[i] *
                           metric_vars.d1_lapse[j];
                FOR2(k, l)
                {
                    sourceE += -metric_vars.lapse * emtensor.Sij[i][j] *
                               gamma_UU[i][k] * gamma_UU[j][l] *
                               metric_vars.K_tensor[k][l];
                }
            }
        }
        sourceE *= det_gamma;

        // assign values of Momentum quantities in the output box
        current_cell.store_vars(rhoM, c_rhoM);
        current_cell.store_vars(Mdot, c_Mdot);
        current_cell.store_vars(sourceM[0], c_SourceM);

        // assign values of Energy quantities in the output box
        current_cell.store_vars(rhoE, c_rhoE);
        current_cell.store_vars(Edot, c_Edot);
        current_cell.store_vars(sourceE, c_SourceE);

        // assign values of Angular Mom quantities in output box
        current_cell.store_vars(rhoJ, c_rhoJ);
        current_cell.store_vars(Jdot, c_Jdot);

        // check BH Mom for fixed BG
        current_cell.store_vars(BHMom[0], c_BHMom);
    }
};

#endif /* FIXEDBGCONSERVEDQUANTITIES_HPP_ */
