/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ISOTROPICKERRFIXEDBG_HPP_
#define ISOTROPICKERRFIXEDBG_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes the initial conditions per arXiv 1401.1548
//! For a highly spinning BH in quasi isotropic coords
class IsotropicKerrFixedBG
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        double mass = 1.0;                      //!<< The mass of the BH
        std::array<double, CH_SPACEDIM> center; //!< The center of the BH
        double spin = 0.0; //!< The spin 'a' in the z direction
    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;
    const params_t m_params;

  protected:
    const double m_dx;

  public:
    IsotropicKerrFixedBG(params_t a_params, double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
        // check this spin param is sensible
        if ((m_params.spin > m_params.mass) || (m_params.spin < -m_params.mass))
        {
            MayDay::Error(
                "The dimensionless spin parameter must be in the range "
                "-1.0 < spin < 1.0");
        }
    }

    /// This just calculates chi which helps with regridding, debug etc
    /// it is only done once on setup as the BG is fixed
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // get position and set vars
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        Vars<data_t> metric_vars;
        compute_metric_background(metric_vars, current_cell);

        // calculate and save chi
        data_t chi = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        chi = pow(chi, -1.0 / 3.0);

        current_cell.store_vars(chi, c_chi);
    }

    /// Refer to Witek et al 1401.1548 for reference for
    /// Quasi Isotropic Kerr and variable conventions used here
    template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Cell<data_t> &current_cell) const
    {
        // where am i?
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

        // black hole params - mass M and spin a
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;

        // work out where we are on the grid
        // R is the quasi isotropic radial coord
        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t R = coords.get_radius();
        const data_t R2 = R * R;
        // the radius in xy plane, subject to a floor
        const data_t rho2 = simd_max(x * x + y * y, 1e-8);
        const data_t rho = sqrt(rho2);

        // useful position quantities
        const data_t cos_theta = z / R;
        const data_t sin_theta = rho / R;
        const data_t cos_theta2 = cos_theta * cos_theta;
        const data_t sin_theta2 = sin_theta * sin_theta;
        const data_t cos_phi = x / rho;
        const data_t sin_phi = y / rho;
        const double r_plus = M + sqrt(M * M - a * a);
        const double r_minus = M - sqrt(M * M - a * a);

        // the Boyer Lindquist coord (eqn (46) 1401.1548)
        const data_t r_BL = R + 0.5 * r_plus + 0.0625 * r_plus * r_plus / R;
        const data_t r_BL2 = r_BL * r_BL;

        // find the quantities in eqns (45)
        const data_t Sigma = r_BL2 + a2 * cos_theta2;
        const data_t Delta = r_BL2 - 2.0 * M * r_BL + a2;
        // In the paper this is script 'A', not to be confused with A_ij
        const data_t AA = (r_BL2 + a2) * (r_BL2 + a2) - Delta * a2 * sin_theta2;

        // now the components per eqn (47) using psi04 = Sigma / R2;
        const data_t gamma_RR = Sigma * (1.0 + 0.25 * r_plus / R) *
                                (1.0 + 0.25 * r_plus / R) /
                                (R * (r_BL - r_minus));
        // gamma_pp = psi04 * AA * R2 / Sigma^2
        const data_t gamma_pp = AA / Sigma * sin_theta2;

        // Need to convert spherical metric to cartesian
        Tensor<2, data_t> spherical_g;
        FOR2(i, j) { spherical_g[i][j] = 0.0; }
        spherical_g[0][0] = gamma_RR;
        spherical_g[1][1] = Sigma; // gamma_tt = psi04 * R2 = Sigma
        spherical_g[2][2] = gamma_pp;

        // jacobian derivatives drdx, dthetadx etc
        using namespace TensorAlgebra;
        const Tensor<1, data_t> x_i = {x, y, z};
        const Tensor<1, data_t> drhodx = {x / rho, y / rho, 0.0};
        Tensor<1, data_t> dRdx = {x / R, y / R, z / R};
        Tensor<2, data_t> jac;
        Tensor<3, data_t> djacdx;
        jac[0][0] = x / R;
        jac[1][0] = cos_phi * z / R2;
        jac[2][0] = -y / rho2;
        jac[0][1] = y / R;
        jac[1][1] = sin_phi * z / R2;
        jac[2][1] = x / rho2;
        jac[0][2] = z / R;
        jac[1][2] = -rho / R2;
        jac[2][2] = 0.0;
        FOR2(i, j)
        {
            djacdx[0][i][j] = (delta(i, j) - x_i[i] * x_i[j] / R2) / R;
        }
        FOR1(i)
        {
            djacdx[1][0][i] =
                1.0 / R2 / rho * (delta(i, 0) * z + delta(i, 2) * x) +
                x * z / R2 / rho * (-drhodx[i] / rho - 2.0 * dRdx[i] / R);
            djacdx[1][1][i] =
                1.0 / R2 / rho * (delta(i, 1) * z + delta(i, 2) * y) +
                y * z / R2 / rho * (-drhodx[i] / rho - 2.0 * dRdx[i] / R);
            djacdx[1][2][i] = -drhodx[i] / R2 + 2.0 * rho * dRdx[i] / R2 / R;
            djacdx[2][0][i] = (-delta(i, 1) + 2.0 * y * drhodx[i] / rho) / rho2;
            djacdx[2][1][i] = (+delta(i, 0) - 2.0 * x * drhodx[i] / rho) / rho2;
            djacdx[2][2][i] = 0.0;
        }

        // derivs wrt R
        const data_t drBLdR = 1.0 - 0.0625 * r_plus * r_plus / R2;
        const data_t dDeltadR = drBLdR * (2.0 * r_BL - 2.0 * M);
        const data_t dSigmadR = 2.0 * r_BL * drBLdR;
        const data_t dAAdR =
            4.0 * drBLdR * r_BL * (r_BL2 + a2) - dDeltadR * a2 * sin_theta2;
        const data_t dgammaRRdR =
            gamma_RR * (dSigmadR / Sigma - 1.0 / R - drBLdR / (r_BL - r_minus) -
                        0.5 * r_plus / R2 / (1.0 + 0.25 * r_plus / R));
        const data_t dgammappdR = gamma_pp * (dAAdR / AA - dSigmadR / Sigma);

        // derivs wrt theta
        const data_t dSigmadtheta = -2.0 * a2 * cos_theta * sin_theta;
        const data_t dAAdtheta = -2.0 * Delta * a2 * sin_theta * cos_theta;
        const data_t dgammaRRdtheta = gamma_RR * dSigmadtheta / Sigma;
        const data_t dgammappdtheta =
            gamma_pp * (dAAdtheta / AA - dSigmadtheta / Sigma) +
            2.0 * cos_theta * sin_theta * AA / Sigma;

        // Calculate the gradients needed (wrt x, y, z)
        Tensor<1, data_t> drBLdx;
        Tensor<1, data_t> dSigmadx;
        Tensor<1, data_t> dDeltadx;
        Tensor<1, data_t> dAAdx;
        Tensor<1, data_t> dgammaRRdx;
        Tensor<1, data_t> dgammappdx;
        FOR1(i)
        {
            drBLdx[i] = jac[0][i] * drBLdR;
            dDeltadx[i] = jac[0][i] * dDeltadR;
            dSigmadx[i] = jac[0][i] * dSigmadR + jac[1][i] * dSigmadtheta;
            dAAdx[i] = jac[0][i] * dAAdR + jac[1][i] * dAAdtheta;
            dgammaRRdx[i] = jac[0][i] * dgammaRRdR + jac[1][i] * dgammaRRdtheta;
            dgammappdx[i] = jac[0][i] * dgammappdR + jac[1][i] * dgammappdtheta;
        }

        // populate ADM vars - lapse and shift
        // use analytic continuation for lapse within horizon
        data_t sign_lapse = (R - 0.25 * r_plus) / abs(R - 0.25 * r_plus);
        vars.lapse = sign_lapse * sqrt(Delta * Sigma / AA);

        // now the shift
        const data_t beta_phi = -2.0 * M * a * r_BL / AA;
        FOR1(i) { vars.shift[i] = 0.0; }
        vars.shift[0] = -y * beta_phi;
        vars.shift[1] = x * beta_phi;

        // spatial metric - convert spherical to cartesian
        FOR2(i, j)
        {
            vars.gamma[i][j] = 0.0;
            FOR2(k, m)
            {
                vars.gamma[i][j] += spherical_g[k][m] * jac[k][i] * jac[m][j];
            }
        }
        const auto gamma_UU = compute_inverse_sym(vars.gamma);

        // Calculate partial derivative of spatial metric
        FOR3(i, j, k)
        {
            vars.d1_gamma[i][j][k] = djacdx[0][j][k] * jac[0][i] * gamma_RR +
                                     djacdx[0][i][k] * jac[0][j] * gamma_RR +
                                     djacdx[1][j][k] * jac[1][i] * Sigma +
                                     djacdx[1][i][k] * jac[1][j] * Sigma +
                                     djacdx[2][j][k] * jac[2][i] * gamma_pp +
                                     djacdx[2][i][k] * jac[2][j] * gamma_pp +
                                     jac[0][i] * jac[0][j] * dgammaRRdx[k] +
                                     jac[1][i] * jac[1][j] * dSigmadx[k] +
                                     jac[2][i] * jac[2][j] * dgammappdx[k];
        }

        // calculate derivs of lapse and shift
        // use analytic continuation of lapse within the horizon
        // (taken care of automatically by use of vars.lapse)
        FOR1(i)
        {
            vars.d1_lapse[i] =
                0.5 * vars.lapse *
                (dDeltadx[i] / Delta + dSigmadx[i] / Sigma - dAAdx[i] / AA);
        }

        FOR2(i, j) { vars.d1_shift[i][j] = 0.0; }
        FOR1(i)
        {
            vars.d1_shift[0][i] =
                vars.shift[0] *
                (drBLdx[i] / r_BL - dAAdx[i] / AA + delta(i, 1) / y);
            vars.d1_shift[1][i] =
                vars.shift[1] *
                (drBLdx[i] / r_BL - dAAdx[i] / AA + delta(i, 0) / x);
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
        vars.K = compute_trace(vars.K_tensor, gamma_UU);
    }

  public:
    // used to decide when to excise - ie when within the horizon of the BH
    // note that this is not templated over data_t
    double excise(const Cell<double> &current_cell) const
    {
        // where am i?
        const Coordinates<double> coords(current_cell, m_dx, m_params.center);

        // black hole params - mass M and boost v
        // "boost" is the gamma factor for the boost
        const double M = m_params.mass;
        const double a = m_params.spin;

        // the quasi isotropic Kerr radius
        const double R = coords.get_radius();

        // compare this to horizon in quasi isotropic Kerr coords
        // which is r+/4
        const double r_horizon = 0.25 * (M + sqrt(M * M - a * a));

        return R / r_horizon;
    }
};

#endif /* ISOTROPICKERRFIXEDBG_HPP_ */
