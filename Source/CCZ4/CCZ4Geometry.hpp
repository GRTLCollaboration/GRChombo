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

        data_t boxtildechi = 0;

        Tensor<2, data_t> covdtilde2chi;
        FOR2(k, l)
        {
            covdtilde2chi[k][l] = d2.chi[k][l];
            FOR1(m) { covdtilde2chi[k][l] -= chris.ULL[m][k][l] * d1.chi[m]; }
        }

        FOR2(k, l) { boxtildechi += h_UU[k][l] * covdtilde2chi[k][l]; }

        data_t dchi_dot_dchi = 0;
        {
            FOR2(m, n) { dchi_dot_dchi += h_UU[m][n] * d1.chi[m] * d1.chi[n]; }
        }

        FOR2(i, j)
        {
            data_t ricci_tilde = 0;
            FOR1(k)
            {
                // Trick: For CCZ4, we can add Z terms to ricci by changing
                // Gamma to chrisvec This way of writing it allows the user to
                // pass Z/chi = {0};
                ricci_tilde += 0.5 * (vars.h[k][i] * d1.Gamma[k][j] +
                                      vars.h[k][j] * d1.Gamma[k][i]);
                ricci_tilde += 0.5 * (vars.Gamma[k] - 2 * Z_over_chi[k]) *
                               (chris.LLL[i][j][k] + chris.LLL[j][i][k]);
                FOR1(l)
                {
                    ricci_tilde -= 0.5 * h_UU[k][l] * d2.h[i][j][k][l];
                    FOR1(m)
                    {
                        ricci_tilde +=
                            h_UU[l][m] *
                            (chris.ULL[k][l][i] * chris.LLL[j][k][m] +
                             chris.ULL[k][l][j] * chris.LLL[i][k][m] +
                             chris.ULL[k][i][m] * chris.LLL[k][l][j]);
                    }
                }
            }

            data_t ricci_chi =
                0.5 * ((GR_SPACEDIM - 2) * covdtilde2chi[i][j] +
                       vars.h[i][j] * boxtildechi -
                       ((GR_SPACEDIM - 2) * d1.chi[i] * d1.chi[j] +
                        GR_SPACEDIM * vars.h[i][j] * dchi_dot_dchi) /
                           (2 * vars.chi));

            data_t z_terms = 0;
            FOR1(k)
            {
                z_terms +=
                    Z_over_chi[k] *
                    (vars.h[i][k] * d1.chi[j] + vars.h[j][k] * d1.chi[i] -
                     vars.h[i][j] * d1.chi[k] + d1.h[i][j][k] * vars.chi);
            }

            out.LL[i][j] =
                (ricci_chi + vars.chi * ricci_tilde + z_terms) / vars.chi;
        }

        out.scalar = vars.chi * TensorAlgebra::compute_trace(out.LL, h_UU);

        return out;
    }

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    static ricci_t<data_t>
    compute_ricci(const vars_t<data_t> &vars,
                  const vars_t<Tensor<1, data_t>> &d1,
                  const diff2_vars_t<Tensor<2, data_t>> &d2,
                  const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris)
    {
        Tensor<1, data_t> Z0 = 0.;
        return compute_ricci_Z(vars, d1, d2, h_UU, chris, Z0);
    }
};

#endif /* CCZ4GEOMETRY_HPP_ */
