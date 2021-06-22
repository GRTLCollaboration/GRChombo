/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "ADMFixedBGVars.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"

class Potential
{
  public:
    struct params_t
    {
        double mass;
        double self_interaction;
    };

    // class params
    params_t m_params;

    // add alias for metric vars
    template <class data_t>
    using MetricVars = typename ADMFixedBGVars::template Vars<data_t>;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the proca field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &dVdA, data_t &dphidt,
                           const vars_t<data_t> &vars,
                           const vars_t<Tensor<1, data_t>> &d1,
                           const MetricVars<data_t> &metric_vars) const
    {
        // calculate full spatial christoffel symbols and gamma^ij
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);

        // for ease of reading
        double c4 = m_params.self_interaction;

        // Here we are defining often used terms
        // DA[i][j] = D_i A_j
        Tensor<2, data_t> DA;
        FOR2(i, j)
        {
            DA[i][j] = d1.Avec[j][i];
            FOR1(k) { DA[i][j] += -chris_phys.ULL[k][i][j] * vars.Avec[k]; }
        }

        // DAScalar = D_i A^i
        data_t DA_scalar;
        DA_scalar = 0;
        FOR2(i, j) { DA_scalar += DA[i][j] * gamma_UU[i][j]; }

        // Xsquared = X^/mu X_/mu
        data_t Xsquared;
        Xsquared = -vars.phi * vars.phi;
        FOR2(i, j) { Xsquared += gamma_UU[i][j] * vars.Avec[j] * vars.Avec[i]; }

        // C = 1 + 4 c4 A^k A_k - 12 c4 phi^2
        data_t C = 1.0 - 12.0 * c4 * vars.phi * vars.phi;
        FOR2(i, j)
        {
            C += 4.0 * c4 * gamma_UU[i][j] * vars.Avec[j] * vars.Avec[i];
        }

        // dVdA = mu^2 ( 1 + 4 c4 (A^k A_k - phi^2))
        dVdA = pow(m_params.mass, 2.0) * (1.0 + 4.0 * c4 * Xsquared);

        // dphidt - for now the whole thing is here since it depends mainly
        // on the form of the potential - except the advection term which is in
        // the ProcaField code
        dphidt = 0;
        FOR2(i, j)
        {
            dphidt += -gamma_UU[i][j] * vars.Avec[i] * metric_vars.d1_lapse[j];
        }
        // QUESTION: Should this be lapse * Z / C  or lapse * Z??
        dphidt += -metric_vars.lapse * vars.Z / C;
        FOR4(i, j, k, l)
        {
            dphidt += -8.0 * c4 * metric_vars.lapse / C * gamma_UU[i][k] *
                      gamma_UU[j][l] * vars.Avec[i] * vars.Avec[j] * DA[k][l];
        }
        dphidt += metric_vars.lapse / C * (1.0 + 4.0 * c4 * Xsquared) *
                  (metric_vars.K * vars.phi - DA_scalar);
        FOR1(i)
        {
            dphidt += 8.0 * c4 * vars.phi * metric_vars.lapse / C *
                      (vars.Evec[i] * vars.Avec[i]);
        }
        FOR4(i, j, k, l)
        {
            dphidt += 8.0 * c4 * vars.phi * metric_vars.lapse / C *
                      (-metric_vars.K_tensor[i][j] * vars.Avec[k] *
                       vars.Avec[l] * gamma_UU[i][k] * gamma_UU[j][l]);
        }
        FOR2(i, j)
        {
            dphidt += 8.0 * c4 * vars.phi * metric_vars.lapse / C *
                      (2.0 * vars.Avec[i] * d1.phi[j] * gamma_UU[i][j]);
        }
    }
};

#endif /* POTENTIAL_HPP_ */
