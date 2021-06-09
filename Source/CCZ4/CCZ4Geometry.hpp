/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This file calculates CCZ4 geometric quantities (or a similar 3+1 split).
#ifndef CCZ4GEOMETRY_HPP_
#define CCZ4GEOMETRY_HPP_

#include "DimensionDefinitions.hpp"
#include "TensorAlgebra.hpp"

//! A structure for the decomposed elements of the Energy Momentum Tensor in
//! 3+1D
template <class data_t> struct emtensor_t
{
    Tensor<2, data_t> Sij; //!< S_ij = T_ij
    Tensor<1, data_t> Si;  //!< S_i = T_ia_n^a
    data_t S;              //!< S = S^i_i
    data_t rho;            //!< rho = T_ab n^a n^b
};

template <class data_t> struct ricci_t
{
    Tensor<2, data_t> LL; // Ricci with two indices down
    data_t scalar;        // Ricci scalar
};

class CCZ4Geometry
{
  protected:
    template <class data_t>
    ALWAYS_INLINE static data_t
    compute_z_terms(const int i, const int j,
                    const Tensor<1, data_t> &Z_over_chi,
                    const Tensor<2, data_t> &h, const Tensor<1, data_t> &d1_chi)
    {
        data_t out = 0.;
        FOR(k)
        {
            out += Z_over_chi[k] * (h[i][k] * d1_chi[j] + h[j][k] * d1_chi[i] -
                                    h[i][j] * d1_chi[k]);
        }
        return out;
    }

  public:
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    static ricci_t<data_t>
    compute_ricci_Z(const vars_t<data_t> &vars,
                    const vars_t<Tensor<1, data_t>> &d1,
                    const diff2_vars_t<Tensor<2, data_t>> &d2,
                    const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris,
                    const Tensor<1, data_t> &Z_over_chi)
    {
        ricci_t<data_t> out;

        Tensor<2, data_t> covdtilde2chi;
        FOR(k, l)
        {
            covdtilde2chi[k][l] = d2.chi[k][l];
            FOR(m) { covdtilde2chi[k][l] -= chris.ULL[m][k][l] * d1.chi[m]; }
        }

        Tensor<3, data_t> chris_LLU = {0.};
        data_t boxtildechi = 0.;
        data_t dchi_dot_dchi = 0.;
        FOR(i, j)
        {
            boxtildechi += covdtilde2chi[i][j] * h_UU[i][j];
            dchi_dot_dchi += d1.chi[i] * d1.chi[j] * h_UU[i][j];
            FOR(k, l) { chris_LLU[i][j][k] += h_UU[k][l] * chris.LLL[i][j][l]; }
        }

        FOR(i, j)
        {
            data_t ricci_hat = 0;
            FOR(k)
            {
                // We call this ricci_hat rather than ricci_tilde as we have
                // replaced what should be \tilde{Gamma} with \hat{Gamma} in
                // order to avoid adding terms that cancel later on
                ricci_hat += 0.5 * (vars.h[k][i] * d1.Gamma[k][j] +
                                    vars.h[k][j] * d1.Gamma[k][i]);
                ricci_hat += 0.5 * vars.Gamma[k] * d1.h[i][j][k];
                FOR(l)
                {
                    ricci_hat += -0.5 * h_UU[k][l] * d2.h[i][j][k][l] +
                                 (chris.ULL[k][l][i] * chris_LLU[j][k][l] +
                                  chris.ULL[k][l][j] * chris_LLU[i][k][l] +
                                  chris.ULL[k][i][l] * chris_LLU[k][j][l]);
                }
            }

            data_t ricci_chi =
                0.5 * ((GR_SPACEDIM - 2) * covdtilde2chi[i][j] +
                       vars.h[i][j] * boxtildechi -
                       ((GR_SPACEDIM - 2) * d1.chi[i] * d1.chi[j] +
                        GR_SPACEDIM * vars.h[i][j] * dchi_dot_dchi) /
                           (2 * vars.chi));

            data_t z_terms = compute_z_terms(i, j, Z_over_chi, vars.h, d1.chi);

            out.LL[i][j] =
                (ricci_chi + vars.chi * ricci_hat + z_terms) / vars.chi;
        }

        out.scalar = vars.chi * TensorAlgebra::compute_trace(out.LL, h_UU);

        return out;
    }

    template <class data_t>
    static Tensor<2, data_t>
    compute_d1_chris_contracted(const Tensor<2, data_t> &h_UU,
                                const Tensor<2, Tensor<1, data_t>> &d1_h,
                                const Tensor<2, Tensor<2, data_t>> &d2_h)
    {
        Tensor<2, data_t> d1_chris_contracted = 0.0;
        FOR(i, j)
        {
            FOR(m, n, p)
            {
                data_t d1_terms = 0.0;
                FOR(q, r)
                {
                    d1_terms += -h_UU[q][r] * (d1_h[n][q][j] * d1_h[m][p][r] +
                                               d1_h[m][n][j] * d1_h[p][q][r]);
                }
                d1_chris_contracted[i][j] +=
                    h_UU[i][m] * h_UU[n][p] * (d2_h[m][n][j][p] + d1_terms);
            }
        }
        return d1_chris_contracted;
    }

    // This function allows adding arbitrary multiples of D_{(i}Z_{j)}
    // to the Ricci scalar rather than the default of 2 in compute_ricci_Z
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    static ricci_t<data_t>
    compute_ricci_Z_general(const vars_t<data_t> &vars,
                            const vars_t<Tensor<1, data_t>> &d1,
                            const diff2_vars_t<Tensor<2, data_t>> &d2,
                            const Tensor<2, data_t> &h_UU,
                            const chris_t<data_t> &chris, const double dZ_coeff)
    {
        // get contributions from conformal metric and factor with zero Z vector
        Tensor<1, data_t> Z0 = 0.;
        auto ricci = compute_ricci_Z(vars, d1, d2, h_UU, chris, Z0);

        // need to add term to correct for d1.Gamma (includes Z contribution)
        // and Gamma in ricci_hat
        auto d1_chris_contracted =
            compute_d1_chris_contracted(h_UU, d1.h, d2.h);
        Tensor<1, data_t> Z_over_chi;
        FOR(i) { Z_over_chi[i] = 0.5 * (vars.Gamma[i] - chris.contracted[i]); }
        FOR(i, j)
        {
            FOR(m)
            {
                // This corrects for the \hat{Gamma}s in ricci_hat
                ricci.LL[i][j] +=
                    (1. - 0.5 * dZ_coeff) * 0.5 *
                    (vars.h[m][i] *
                         (d1_chris_contracted[m][j] - d1.Gamma[m][j]) +
                     vars.h[m][j] *
                         (d1_chris_contracted[m][i] - d1.Gamma[m][i]) +
                     (chris.contracted[m] - vars.Gamma[m]) * d1.h[i][j][m]);
            }
            data_t z_terms = compute_z_terms(i, j, Z_over_chi, vars.h, d1.chi);
            ricci.LL[i][j] += 0.5 * dZ_coeff * z_terms / vars.chi;
        }
        ricci.scalar = vars.chi * TensorAlgebra::compute_trace(ricci.LL, h_UU);
        return ricci;
    }

    // This function returns the pure Ricci scalar with no contribution from the
    // Z vector - used e.g. in the constraint calculations.
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    static ricci_t<data_t>
    compute_ricci(const vars_t<data_t> &vars,
                  const vars_t<Tensor<1, data_t>> &d1,
                  const diff2_vars_t<Tensor<2, data_t>> &d2,
                  const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris)
    {
        return compute_ricci_Z_general(vars, d1, d2, h_UU, chris, 0.0);
    }
};

#endif /* CCZ4GEOMETRY_HPP_ */
