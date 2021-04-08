/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOOSTEDKERRSCHILDBHFIXEDBG_HPP_
#define BOOSTEDKERRSCHILDBHFIXEDBG_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes the initial conditions per arXiv 1401.1548
class BoostedKerrSchildBHFixedBG
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        double mass = 1.0;                      //!<< The mass of the BH
        std::array<double, CH_SPACEDIM> center; //!< The center of the BH
        double velocity = 0.0; //!< The boost velocity in the x direction
    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;

    const params_t m_params;
    const double m_dx;

    BoostedKerrSchildBHFixedBG(params_t a_params, double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
        // check this boost param is sensible
        if ((m_params.velocity > 1.0) || (m_params.velocity < -1.0))
        {
            MayDay::Error("The boost velocity parameter must be in the range "
                          "-1.0 < velocity < 1.0");
        }
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
    /// NB we use x' = x - vt to prevent the movement of the BH
    /// on the grid
    /// x_p = gamma_boost * x'
    /// but we are in the rest frame of the SF (the BH is boosted)
    /// NB2 eqn 11 is missing a squared under the sqrt
    template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Cell<data_t> &current_cell) const
    {
        // where am i?
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

        // black hole params - mass M and boost v
        // "boost" is the gamma factor for the boost
        const double M = m_params.mass;
        const double v = m_params.velocity;
        const double boost = pow(1 - v * v, -0.5);

        // work out where we are on the grid including effect of boost
        // on x direction (length contraction)
        const data_t x = coords.x;
        const data_t x_p = coords.x * boost;
        const double y = coords.y;
        const double z = coords.z;

        // the Kerr Schild radius (boosted)
        const data_t r2 = x_p * x_p + y * y + z * z;
        const data_t r = sqrt(r2);

        // find the H and el quantities (el decomposed into space and time)
        const data_t H = M / r;
        const Tensor<1, data_t> el = {boost * (x_p / r - v), y / r, z / r};
        const data_t el_t = boost * (1.0 - v * x_p / r);

        // Calculate the gradients in el and H
        Tensor<1, data_t> dHdx;
        Tensor<1, data_t> dltdx;
        Tensor<2, data_t> dldx;
        get_KS_derivs(dHdx, dldx, dltdx, H, coords);

        // populate ADM vars
        vars.lapse = pow(1.0 + 2.0 * H * el_t * el_t, -0.5);
        FOR2(i, j)
        {
            vars.gamma[i][j] =
                TensorAlgebra::delta(i, j) + 2.0 * H * el[i] * el[j];
        }
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(vars.gamma);
        FOR1(i)
        {
            vars.shift[i] = 0;
            FOR1(j)
            {
                vars.shift[i] += gamma_UU[i][j] * 2.0 * H * el[j] * el_t;
            }
        }
        // this adjustment gives the shift which achieves x' = x - vt
        vars.shift[0] += v;

        // Calculate partial derivative of spatial metric
        FOR3(i, j, k)
        {
            vars.d1_gamma[i][j][k] =
                2.0 * (el[i] * el[j] * dHdx[k] + H * el[i] * dldx[j][k] +
                       H * el[j] * dldx[i][k]);
        }

        // calculate derivs of lapse and shift
        FOR1(i)
        {
            vars.d1_lapse[i] = -pow(vars.lapse, 3.0) * el_t *
                               (el_t * dHdx[i] + 2.0 * H * dltdx[i]);
        }

        // use the fact that shift^i = lapse^2 * shift_i + v^i
        // and v^i is a constant vector
        FOR2(i, j)
        {
            vars.d1_shift[i][j] =
                2.0 * el_t * dHdx[j] * pow(vars.lapse, 2.0) * el[i] +
                4.0 * el_t * H * vars.lapse * vars.d1_lapse[j] * el[i] +
                2.0 * el_t * H * pow(vars.lapse, 2.0) * dldx[i][j] +
                2.0 * dltdx[j] * H * pow(vars.lapse, 2.0) * el[i];
        }

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

  protected:
    /// Work out the gradients of the quantities H and el appearing in the Kerr
    /// Schild solution
    template <class data_t>
    void get_KS_derivs(Tensor<1, data_t> &dHdx, Tensor<2, data_t> &dldx,
                       Tensor<1, data_t> &dltdx, const data_t &H,
                       const Coordinates<data_t> &coords) const
    {
        // black hole params - mass M and boost v
        const double M = m_params.mass;
        const double v = m_params.velocity;
        const double boost = pow(1 - v * v, -0.5);

        // work out where we are on the grid including effect of boost
        // on x direction (length contraction)
        Tensor<1, data_t> x;
        x[0] = coords.x;
        x[1] = coords.y;
        x[2] = coords.z;
        const data_t x_p = coords.x * boost;

        // the coordinate radius (boosted), subject to a floor
        const data_t r2 = x_p * x_p + x[1] * x[1] + x[2] * x[2];
        const data_t r = sqrt(r2);

        using namespace TensorAlgebra;

        // derivatives of r wrt actual grid coords
        Tensor<1, data_t> drdx;
        FOR1(i) { drdx[i] = x[i] / r; }
        drdx[0] *= boost * boost;

        FOR1(i) { dHdx[i] = -H / r * drdx[i]; }

        // note to use convention as in rest of tensors the last index is the
        // derivative index so these are d_j l_i
        FOR2(i, j)
        {
            dldx[i][j] = -x[i] / r2 * drdx[j] + delta(i, j) / r;
            if (i == 0)
            {
                dldx[i][j] *= boost * boost;
            }
        }

        // then dltdi
        FOR1(i)
        {
            dltdx[i] =
                v * boost * boost * (-delta(0, i) / r + x[0] / r2 * drdx[i]);
        }
    }

  public:
    // used to decide when to excise - ie when within the horizon of the BH
    // note that this is not templated over data_t
    double excise(const Cell<double> &current_cell) const
    {
        // black hole params - mass M and boost v
        // "boost" is the gamma factor for the boost
        const Coordinates<double> coords(current_cell, m_dx, m_params.center);
        const double M = m_params.mass;
        const double boost =
            pow(1.0 - m_params.velocity * m_params.velocity, -0.5);

        // work out where we are on the grid including effect of boost
        // on x direction (length contraction)
        const double x_p = coords.x * boost;
        const double y = coords.y;
        const double z = coords.z;

        // the coordinate radius (boosted)
        const double r2 = x_p * x_p + y * y + z * z;

        // compare this to horizon in kerr schild coords
        const double r_horizon = 2.0 * M;

        return sqrt(r2) / r_horizon;
    }
};

#endif /* BOOSTEDKERRSCHILDBHFIXEDBG_HPP_ */
