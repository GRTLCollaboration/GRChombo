/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef METRICBG_HPP_
#define METRICBG_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "InitialDataTools.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes the metric BG to evolve on
class MetricBG
{
  public:
    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;
    const double m_dx;
    const double m_L;
    const std::array<double, CH_SPACEDIM> m_center;

    MetricBG(const double a_L, const double a_dx,
             const std::array<double, CH_SPACEDIM> a_center)
        : m_dx(a_dx), m_center(a_center), m_L(a_L)
    {
    }

    /// This just calculates chi which helps with regridding, debug etc
    /// it is only done once on setup as the BG is fixed
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // get position and set vars
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        Vars<data_t> metric_vars;
        compute_metric_background(metric_vars, current_cell);

        // calculate and save chi
        data_t chi = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        chi = pow(chi, -1.0 / 3.0);
        current_cell.store_vars(chi, c_chi);
    }

    // Kerr Schild solution
    template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Cell<data_t> &current_cell) const
    {
        using namespace TensorAlgebra;
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        const Tensor<1, data_t> x = {coords.x, coords.y, coords.z};
        const data_t r = coords.get_radius();

        // populate ADM vars
        vars.lapse = 1.0;
        FOR2(i, j) { vars.gamma[i][j] = delta(i, j); }
        const auto gamma_UU = compute_inverse_sym(vars.gamma);
        FOR1(i) { vars.shift[i] = 0.0; }

        // Calculate partial derivative of spatial metric
        FOR3(i, j, k) { vars.d1_gamma[i][j][k] = 0.0; }

        // calculate derivs of lapse and shift
        FOR1(i) { vars.d1_lapse[i] = 0.0; }

        // use the fact that shift^i = lapse^2 * shift_i
        FOR2(i, j) { vars.d1_shift[i][j] = 0.0; }

        // calculate the extrinsic curvature, using the fact that
        // 2 * lapse * K_ij = D_i \beta_j + D_j \beta_i - dgamma_ij dt
        // and dgamma_ij dt = 0 in chosen fixed gauge
        const auto chris_phys = compute_christoffel(vars.d1_gamma, gamma_UU);
        FOR2(i, j)
        {
            vars.K_tensor[i][j] = 0.0;
            FOR1(k)
            {
                vars.K_tensor[i][j] +=
                    vars.gamma[k][j] * vars.d1_shift[k][i] +
                    vars.gamma[k][i] * vars.d1_shift[k][j] +
                    (vars.d1_gamma[k][i][j] + vars.d1_gamma[k][j][i]) *
                        vars.shift[k];
                FOR1(m)
                {
                    vars.K_tensor[i][j] += -2.0 * chris_phys.ULL[k][i][j] *
                                           vars.gamma[k][m] * vars.shift[m];
                }
            }
            vars.K_tensor[i][j] *= 0.5 / vars.lapse;
        }
        vars.K = compute_trace(gamma_UU, vars.K_tensor);
    }
};

#endif /* METRICBG_HPP_ */
