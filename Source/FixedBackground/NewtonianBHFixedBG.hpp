/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef NEWTONIANBHFIXEDBG_HPP_
#define NEWTONIANBHFIXEDBG_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes the background metric for a newtonian potential
class NewtonianBHFixedBG
{
  public:
    //! Struct for the params of the object
    struct params_t
    {
        double mass = 1.0;                      //!<< The mass of the star
        std::array<double, CH_SPACEDIM> center; //!< The center of the star
    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;

    const params_t m_params;
    const double m_dx;
    const double m_time;
    const double m_omega;
    const double m_separation;

    NewtonianBHFixedBG(const params_t a_params, const double a_dx,
                       const double a_time = 0.0,
                       const double a_separation = 0.0,
                       const double a_omega = 0.0)
        : m_params(a_params), m_dx(a_dx), m_time(a_time), m_omega(a_omega),
          m_separation(a_separation)
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

    /// Newtonian solution ds^2 =
    ///       -(1 + 2 Phi) dt^2 + (1 - 2 Phi) \delta_ij dx^i dx^j
    /// Phi = - M / r
    template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Cell<data_t> &current_cell) const
    {
        // where am i?
        std::array<double, CH_SPACEDIM> center_now;
        if (m_time != 0.0)
        {
            double sign = m_time / abs(m_time);
            center_now[0] =
                m_params.center[0] -
                sign * m_separation * (1 - cos(m_omega * abs(m_time)));
            center_now[1] = m_params.center[1] +
                            sign * m_separation * sin(m_omega * abs(m_time));
            center_now[2] = m_params.center[2];
        }
        else
        {
            center_now = m_params.center;
        }
        const Coordinates<data_t> coords(current_cell, m_dx, center_now);

        // object mass M
        const double M = m_params.mass;

        // work out where we are on the grid
        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;

        // the coordinate radius
        const data_t r = coords.get_radius();

        // find the potential
        data_t Phi = -M / r;

        // Calculate the gradients in el and H
        Tensor<1, data_t> dPhidx;
        get_potential_derivs(dPhidx, Phi, coords);

        // populate ADM vars
        using namespace TensorAlgebra;
        const data_t sign = (r - 2.0 * M) / abs(r - 2.0 * M);
        vars.lapse = sign * sqrt(abs(1.0 + 2.0 * Phi));

        FOR2(i, j) { vars.gamma[i][j] = (1.0 - 2.0 * Phi) * delta(i, j); }
        FOR1(i) { vars.shift[i] = 0; }

        // Calculate partial derivative of spatial metric
        FOR3(i, j, k)
        {
            vars.d1_gamma[i][j][k] = -2.0 * dPhidx[k] * delta(i, j);
        }

        // calculate derivs of lapse and shift
        FOR1(i) { vars.d1_lapse[i] = dPhidx[i] / vars.lapse; }

        FOR2(i, j) { vars.d1_shift[i][j] = 0.0; }

        // calculate the extrinsic curvature
        FOR2(i, j) { vars.K_tensor[i][j] = 0.0; }
        vars.K = 0.0;
    }

  protected:
    /// Work out the gradients of the quantity Phi
    template <class data_t>
    void get_potential_derivs(Tensor<1, data_t> &dPhidx, const data_t &Phi,
                              const Coordinates<data_t> &coords) const
    {
        // black hole params - mass M
        const double M = m_params.mass;

        // work out where we are on the grid
        Tensor<1, data_t> x;
        x[0] = coords.x;
        x[1] = coords.y;
        x[2] = coords.z;

        // the coordinate radius, subject to a floor
        const data_t r = coords.get_radius();
        const auto inside_horizon = simd_compare_lt(r, 2 * M);

        // derivatives of Phi wrt grid coords
        data_t dPhidr = -Phi / r;
        FOR1(i) { dPhidx[i] = dPhidr * x[i] / r; }
    }

  public:
    // used to decide when to excise
    double excise(const Cell<double> &current_cell) const
    {
        // black hole params - where am i now?
        double sign = m_time / abs(m_time);
        const double M = m_params.mass;
        std::array<double, CH_SPACEDIM> center_now;
        center_now[0] = m_params.center[0] -
                        sign * m_separation * (1 - cos(m_omega * m_time));
        center_now[1] = m_params.center[1] +
                        sign * m_separation * sin(m_omega * abs(m_time));
        center_now[2] = m_params.center[2];
        const Coordinates<double> coords(current_cell, m_dx, center_now);

        // work out where we are on the grid
        const double x = coords.x;
        const double y = coords.y;
        const double z = coords.z;

        // the coordinate radius (boosted)
        const double r = coords.get_radius();

        // compare this to horizon coords
        const double r_horizon = 2.0 * M;

        return r / r_horizon;
    }
};

#endif /* NEWTONIANBHFIXEDBG_HPP_ */
