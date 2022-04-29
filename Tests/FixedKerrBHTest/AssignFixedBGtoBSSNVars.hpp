/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ASSIGNFIXEDBGTOBSSNVARS_HPP_
#define ASSIGNFIXEDBGTOBSSNVARS_HPP_

#include "ADMFixedBGVars.hpp"
#include "BoxLoops.hpp"
#include "CCZ4.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

template <class background_t> class AssignFixedBGtoBSSNVars
{
  public:
    AssignFixedBGtoBSSNVars(const background_t a_background, const double a_dx,
                            const std::array<double, CH_SPACEDIM> a_center)
        : m_dx(a_dx), m_center(a_center), m_background(a_background)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // compute the ADM vars
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        ADMFixedBGVars::Vars<data_t> adm_vars;
        m_background.compute_metric_background(adm_vars, coords);

        // assign values to the BSSN vars
        using namespace TensorAlgebra;
        CCZ4::Vars<data_t> bssn_vars;
        bssn_vars.K = adm_vars.K;
        bssn_vars.lapse = adm_vars.lapse;
        FOR1(i) { bssn_vars.shift[i] = adm_vars.shift[i]; }
        bssn_vars.chi = compute_determinant_sym(adm_vars.gamma);
        bssn_vars.chi = pow(bssn_vars.chi, -1.0 / 3.0);
        FOR2(i, j)
        {
            bssn_vars.h[i][j] = adm_vars.gamma[i][j] * bssn_vars.chi;
            bssn_vars.A[i][j] = adm_vars.K_tensor[i][j] * bssn_vars.chi -
                                1.0 / 3.0 * bssn_vars.K * bssn_vars.h[i][j];
        }

        // set the field to something arbitrary
        data_t r = coords.get_radius();
        current_cell.store_vars(0.1 * sin(coords.x) * exp(-r) * coords.z,
                                c_phi);
        current_cell.store_vars(0.1 * cos(coords.y) * exp(-r) * coords.z, c_Pi);
        current_cell.store_vars(bssn_vars);
    }

  protected:
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;
    const background_t m_background;
};

#endif /* ASSIGNFIXEDBGTOBSSNVARS_HPP_ */
