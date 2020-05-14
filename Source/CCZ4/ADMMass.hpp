/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ADMMASS_HPP_
#define ADMMASS_HPP_

#include "ADMConformalVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the ADM massa
class ADMMass
{
    // Use the variable definition in ADMConformalVars - only require the key
    // vars
    template <class data_t> using Vars = ADMConformalVars::VarsNoGauge<data_t>;

    template <class data_t>
    using Diff1Vars = ADMConformalVars::Diff2VarsNoGauge<data_t>;

  public:
    enum DIR
    {
        NONE,
        X,
        Y,
        Z
    };

    ADMMass(const std::array<double, CH_SPACEDIM> &a_center, double a_dx,
            DIR spin_direction = Z, double a_G_Newton = 1.0)
        : m_deriv(a_dx), m_center(a_center), m_G_Newton(a_G_Newton),
          m_dir(spin_direction)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and calc 1st
        // derivs
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Diff1Vars>(current_cell);

        using namespace TensorAlgebra;
        const auto h_UU = compute_inverse_sym(vars.h);

        // Surface element for integration
        Coordinates<data_t> coords(current_cell, m_deriv.m_dx, m_center);
        data_t r = coords.get_radius();
        Tensor<1, data_t> x = {coords.x, coords.y, coords.z};
        // This is multiplied by r^2 as SphericalExtraction assumes it is
        // normalised as such.
        Tensor<1, data_t> dS_U = x;

        data_t dS_norm = 0.;
        FOR2(i, j) { dS_norm += vars.h[i][j] / vars.chi * dS_U[i] * dS_U[j]; }
        dS_norm = sqrt(dS_norm);
        FOR1(i) { dS_U[i] /= dS_norm; }

        Tensor<1, data_t> dS_L;
        FOR1(i)
        {
            // dS_L[i] = dS_U[i];
            dS_L[i] = 0.;
            FOR1(j) { dS_L[i] += vars.h[i][j] / vars.chi * dS_U[j]; }
        }

        data_t Madm = 0.0;
        FOR4(i, j, k, l)
        {
            Madm += dS_L[i] / (16. * M_PI * m_G_Newton) * pow(vars.chi, -1.5) *
                    h_UU[j][k] * h_UU[i][l] *
                    (vars.chi * (d1.h[l][k][j] - d1.h[j][k][l]) -
                     (vars.h[l][k] * d1.chi[j] - vars.h[j][k] * d1.chi[l]));
        }

        // assign values of ADMMass in output box
        current_cell.store_vars(Madm, c_Madm);

        if (m_dir == NONE)
            return;

        // user should be able to run this just for the ADM mass
        auto int_Jadm = ChomboParameters::variable_name_to_enum("Jadm");
        if (int_Jadm < 0)
            MayDay::Error("Please include 'c_Jadm' in UserVariables with name "
                          "'Jadm'");

        // spin about z axis
        data_t Jadm = 0.0;

        // note this is the levi civita symbol,
        // not tensor (eps_tensor = eps_symbol * chi^-1.5)
        const Tensor<3, double> epsilon = TensorAlgebra::epsilon();

        FOR3(i, j, k)
        {
            Jadm += -dS_L[i] / (8. * M_PI * m_G_Newton) *
                    epsilon[m_dir - 1][j][k] * x[j] * vars.K *
                    TensorAlgebra::delta(i, k);

            FOR2(l, m)
            {
                Jadm += dS_L[i] / (8. * M_PI * m_G_Newton) *
                        epsilon[m_dir - 1][j][k] * x[j] * h_UU[i][l] *
                        h_UU[k][m] * vars.chi *
                        (vars.A[l][m] + vars.K * vars.h[l][m] / 3.);
            }
        }
        current_cell.store_vars(Jadm, int_Jadm);
    }

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const std::array<double, CH_SPACEDIM> &m_center;
    const double m_G_Newton; //!< Newton's constant

    DIR m_dir;
};

#endif /* ADMMASS_HPP_ */
