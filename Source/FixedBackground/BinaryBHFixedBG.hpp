/*
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BINARYBHFIXEDBG_HPP_
#define BINARYBHFIXEDBG_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "InitialDataTools.hpp"
#include "IsotropicKerrFixedBG.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes the initial conditions per arXiv 1401.1548
class BinaryBHFixedBG
{
  public:
    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;

  protected:
    const IsotropicKerrFixedBG m_bh_1;
    const IsotropicKerrFixedBG m_bh_2;
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;

  public:
    BinaryBHFixedBG(IsotropicKerrFixedBG::params_t a_params_1,
                    IsotropicKerrFixedBG::params_t a_params_2, double a_dx,
                    std::array<double, CH_SPACEDIM> a_center)
        : m_bh_1(a_params_1, a_dx), m_bh_2(a_params_2, a_dx), m_dx(a_dx),
          m_center(a_center)
    {
    }

    /// This just calculates chi which helps with regridding, debug etc
    /// it is only done once on setup as the BG is fixed
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // get position and set vars
        Vars<data_t> metric_vars;
        compute_metric_background(metric_vars, current_cell);

        // calculate and save chi
        data_t chi = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        chi = pow(chi, -1.0 / 3.0);

        current_cell.store_vars(chi, c_chi);
    }

    /// Schwarzschild boosted solution per gr-qc 9805023
    /// Use superposition of two solutions
    template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Cell<data_t> &current_cell) const
    {
        const Coordinates<data_t> coords1(current_cell, m_dx,
                                          m_bh_1.m_params.center);
        const Coordinates<data_t> coords2(current_cell, m_dx,
                                          m_bh_2.m_params.center);
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        Vars<data_t> vars_1;
        Vars<data_t> vars_2;

        m_bh_1.compute_metric_background(vars_1, current_cell);
        m_bh_2.compute_metric_background(vars_2, current_cell);

        // combine the ADM vars, using analytic continuation of lapse
        vars.lapse =
            vars_1.lapse * vars_1.lapse + vars_2.lapse * vars_2.lapse - 1.0;
        data_t sign_of_lapse = vars.lapse / abs(vars.lapse);
        vars.lapse = sqrt(abs(vars.lapse)) * sign_of_lapse;

        using namespace TensorAlgebra;
        FOR2(i, j)
        {
            vars.gamma[i][j] =
                vars_1.gamma[i][j] + vars_2.gamma[i][j] - delta(i, j);
        }

        FOR1(i) { vars.shift[i] = vars_1.shift[i] + vars_2.shift[i]; }

        // Calculate partial derivative of spatial metric
        FOR3(i, j, k)
        {
            vars.d1_gamma[i][j][k] =
                vars_1.d1_gamma[i][j][k] + vars_2.d1_gamma[i][j][k];
        }

        // calculate derivs of lapse and shift
        FOR1(i)
        {
            vars.d1_lapse[i] = 0.5 / vars.lapse *
                               (vars_1.d1_lapse[i] * 2.0 * vars_1.lapse +
                                vars_2.d1_lapse[i] * 2.0 * vars_2.lapse);
        }
        FOR2(i, j)
        {
            vars.d1_shift[i][j] = vars_1.d1_shift[i][j] + vars_2.d1_shift[i][j];
        }

        // calculate the extrinsic curvature, using the fact that
        // 2 * lapse * K_ij = D_i \beta_j + D_j \beta_i - dgamma_ij dt
        // and dgamma_ij dt = 0 in chosen fixed gauge
        const auto gamma_UU = compute_inverse_sym(vars.gamma);
        const auto chris_phys = compute_christoffel(vars.d1_gamma, gamma_UU);
        FOR2(i, j)
        {
            vars.K_tensor[i][j] = 0.0;
            vars.K_tensor[i][j] +=
                vars_1.K_tensor[i][j] + vars_2.K_tensor[i][j];
        }
        vars.K = compute_trace(gamma_UU, vars.K_tensor);
    }

    // used to decide when to excise - ie when within the horizon of the BH
    // note that this is not templated over data_t
    double excise(const Cell<double> &current_cell) const
    {
        // compare this to horizon in kerr schild coords
        const double r_bh1 = m_bh_1.excise(current_cell);
        const double r_bh2 = m_bh_2.excise(current_cell);

        // return the minimum value, as we care only if
        // either is less than 1.0
        return min(r_bh1, r_bh2);
    }
};

#endif /* BINARYBHFIXEDBG_HPP_ */
