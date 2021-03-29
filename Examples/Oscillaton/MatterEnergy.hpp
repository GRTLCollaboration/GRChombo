/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MATTERENERGY_HPP_
#define MATTERENERGY_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "CoordinateTransformations.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "MatterCCZ4.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the MatterEnergy rho with type matter_t and writes it to the grid
template <class matter_t> class MatterEnergy
{
    // Use the variable definition in CCZ4
    template <class data_t>
    using Vars = typename MatterCCZ4<matter_t>::template Vars<data_t>;

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const matter_t my_matter; //!< The matter object
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_dx;

  public:
    MatterEnergy(matter_t a_matter, double a_dx,
                 std::array<double, CH_SPACEDIM> a_center)
        : my_matter(a_matter), m_dx(a_dx), m_deriv(a_dx), m_center(a_center)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, derivs
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);

        // set up the metric quantities needed
        using namespace TensorAlgebra;
        const auto h_UU = compute_inverse_sym(vars.h);
        const auto chris = compute_christoffel(d1.h, h_UU);
        const emtensor_t<data_t> emtensor =
            my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);
        const data_t det_gamma = pow(vars.chi, -3.0);
        Tensor<2, data_t> vars_gamma, vars_Kij;
        FOR2(i, j)
        {
            vars_gamma[i][j] = vars.h[i][j] / vars.chi;
            vars_Kij[i][j] =
                1.0 / vars.chi * (vars.A[i][j] + 1.0 / 3.0 * vars.h[i][j]);
        }
        const auto gamma_UU = compute_inverse_sym(vars_gamma);
        const Tensor<3, data_t> chris_phys =
            compute_phys_chris(d1.chi, vars.chi, vars.h, h_UU, chris.ULL);
        // coordinate quantities
        // The unit covector normal to the spherical surface
        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        Tensor<1, data_t> si_L, Ni_L;
        data_t R = coords.get_radius();
        si_L[0] = coords.x / R;
        si_L[1] = coords.y / R;
        si_L[2] = coords.z / R;
        // normalise this to 1 using spatial metric
        data_t mod_N2 = 0.0;
        FOR2(i, j) { mod_N2 += gamma_UU[i][j] * si_L[i] * si_L[j]; }
        FOR1(i) { Ni_L[i] = si_L[i] / sqrt(mod_N2); }

        // the area element of the sphere - check matches up
        data_t rxy2 =
            simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
        data_t r2sintheta = sqrt(rxy2) * R;
        using namespace CoordinateTransformations;
        Tensor<2, data_t> spherical_gamma =
            cartesian_to_spherical_LL(vars_gamma, coords.x, coords.y, coords.z);
        data_t sqrt_det_Sigma = area_element_sphere(spherical_gamma);

        // calculate according to Landau Lifshitz method
        data_t g = vars.lapse * vars.lapse * det_gamma;
        data_t rho1 = sqrt_det_Sigma / r2sintheta * Ni_L[0]; // CHECK THE SAME
        data_t flux1 = sqrt(det_gamma) * si_L[0];            // AS ME
        data_t source1 = 0.0; // ADD EXPRESSIONS HERE

        // calculate according to Killing vector method
        data_t rho2 = -emtensor.rho * vars.lapse;
        FOR1(i) { rho2 += -emtensor.Si[i] * vars.shift[i]; }
        rho2 *= sqrt(det_gamma);

        data_t flux2 = 0.0;
        FOR1(i)
        {
            flux2 += vars.lapse * Ni_L[i] * emtensor.rho * vars.shift[i];
            FOR1(j)
            {
                flux2 += -Ni_L[i] * emtensor.Si[j] *
                         (vars.shift[i] * vars.shift[j] +
                          vars.lapse * vars.lapse * gamma_UU[i][j]);
                FOR1(k)
                {
                    flux2 += Ni_L[i] * vars.lapse * gamma_UU[i][j] *
                             vars.shift[k] * emtensor.Sij[j][k];
                }
            }
        }
        flux2 *= sqrt_det_Sigma / r2sintheta;

        data_t dlapsedt = -2.0 * vars.lapse * vars.K;
        Tensor<1, data_t> dbetadt;
        FOR1(i) { dbetadt[i] = vars.B[i]; }
        data_t source2 = -emtensor.rho * dlapsedt;
        FOR1(i)
        {
            source2 += emtensor.Si[i] * dbetadt[i];
            FOR1(j)
            {
                source2 += emtensor.Si[i] * (-vars.shift[i] / vars.lapse *
                                                 vars.shift[j] * d1.lapse[j] +
                                             vars.shift[j] * d1.shift[i][j]);
                FOR1(k)
                {
                    source2 +=
                        emtensor.Si[i] * vars.shift[j] *
                            (-vars.lapse * gamma_UU[i][k] * vars_Kij[k][j] +
                             vars.shift[i] * vars.shift[k] * vars_Kij[j][k] /
                                 vars.lapse +
                             vars.shift[k] * chris_phys[i][j][k]) +
                        emtensor.Sij[k][j] * gamma_UU[i][k] * vars.lapse *
                            d1.shift[j][i];
                    FOR1(l)
                    {
                        source2 += -vars.lapse * vars.lapse * gamma_UU[i][k] *
                                       gamma_UU[j][l] * emtensor.Sij[k][l] *
                                       vars_Kij[i][j] +
                                   vars.lapse * gamma_UU[i][k] *
                                       emtensor.Sij[k][j] * vars.shift[l] *
                                       chris_phys[j][l][i];
                    }
                }
            }
        }
        source2 *= sqrt(det_gamma);

        // assign values of MatterEnergy in output box
        current_cell.store_vars(rho1, c_rho1);
        current_cell.store_vars(flux1, c_flux1);
        current_cell.store_vars(source1, c_source1);
        current_cell.store_vars(rho2, c_rho2);
        current_cell.store_vars(flux2, c_flux2);
        current_cell.store_vars(source2, c_source2);
    }
};

#endif /* MATTERENERGY_HPP_ */
