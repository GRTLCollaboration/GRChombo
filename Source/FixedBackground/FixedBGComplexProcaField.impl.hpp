/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(FIXEDBGCOMPLEXPROCAFIELD_HPP_)
#error "This file should only be included through FixedBGComplexProcaField.hpp"
#endif

#ifndef FIXEDBGCOMPLEXPROCAFIELD_IMPL_HPP_
#define FIXEDBGCOMPLEXPROCAFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> FixedBGComplexProcaField::compute_emtensor(
    const vars_t<data_t> &vars, const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<2, data_t> &gamma_UU,
    const Tensor<3, data_t> &chris_phys_ULL) const
{
    emtensor_t<data_t> out;

    // Some useful quantities
    const double msquared = pow(m_vector_mass, 2.0);

    // Real parts first then Im
    // D_i A_j
    Tensor<2, data_t> DA_Re;
    Tensor<2, data_t> DA_Im;
    FOR2(i, j)
    {
        DA_Re[i][j] = d1.Avec_Re[j][i];
        DA_Im[i][j] = d1.Avec_Im[j][i];
        FOR1(k)
        {
            DA_Re[i][j] += -chris_phys_ULL[k][i][j] * vars.Avec_Re[k];
            DA_Im[i][j] += -chris_phys_ULL[k][i][j] * vars.Avec_Im[k];
        }
    }

    // D_i A_j - D_j A_i (NB exterior derivative, so christoffel symbols cancel)
    Tensor<2, data_t> diff_DA_Re;
    Tensor<2, data_t> diff_DA_Im;
    FOR2(i, j)
    {
        diff_DA_Re[i][j] = d1.Avec_Re[j][i] - d1.Avec_Re[i][j];
        diff_DA_Im[i][j] = d1.Avec_Im[j][i] - d1.Avec_Im[i][j];
    }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i, j)
    {
        out.Sij[i][j] =
            // Re part
            msquared * (vars.Avec_Re[i] * vars.Avec_Re[j] +
                        0.5 * metric_vars.gamma[i][j] * vars.Avec0_Re *
                            vars.Avec0_Re) +
            // Im part
            msquared *
                (vars.Avec_Im[i] * vars.Avec_Im[j] +
                 0.5 * metric_vars.gamma[i][j] * vars.Avec0_Im * vars.Avec0_Im);

        FOR2(k, l)
        {
            out.Sij[i][j] +=
                // Re part
                gamma_UU[k][l] * (diff_DA_Re[i][k] * diff_DA_Re[j][l]) -
                metric_vars.gamma[i][k] * metric_vars.gamma[j][l] *
                    vars.Evec_Re[k] * vars.Evec_Re[l] +
                0.5 * metric_vars.gamma[k][l] * metric_vars.gamma[i][j] *
                    vars.Evec_Re[k] * vars.Evec_Re[l] -
                0.5 * gamma_UU[k][l] * metric_vars.gamma[i][j] * msquared *
                    vars.Avec_Re[k] * vars.Avec_Re[l] +
                // Im part
                gamma_UU[k][l] * (diff_DA_Im[i][k] * diff_DA_Im[j][l]) -
                metric_vars.gamma[i][k] * metric_vars.gamma[j][l] *
                    vars.Evec_Im[k] * vars.Evec_Im[l] +
                0.5 * metric_vars.gamma[k][l] * metric_vars.gamma[i][j] *
                    vars.Evec_Im[k] * vars.Evec_Im[l] -
                0.5 * gamma_UU[k][l] * metric_vars.gamma[i][j] * msquared *
                    vars.Avec_Im[k] * vars.Avec_Im[l];

            FOR2(m, n)
            {
                out.Sij[i][j] += -0.5 * metric_vars.gamma[i][j] *
                                 gamma_UU[k][m] * gamma_UU[l][n] *
                                 (DA_Re[k][l] * diff_DA_Re[m][n] +
                                  DA_Im[k][l] * diff_DA_Im[m][n]);
            }
        }
    }

    // S = Tr_S_ij
    out.S = 0.0;
    FOR2(i, j) { out.S += out.Sij[i][j] * gamma_UU[i][j]; }

    // S_i (note lower index) = n^a T_a0
    FOR1(i)
    {
        out.Si[i] = msquared * vars.Avec0_Re * vars.Avec_Re[i] +
                    msquared * vars.Avec0_Im * vars.Avec_Im[i];

        FOR1(j)
        {
            out.Si[i] += vars.Evec_Re[j] * diff_DA_Re[i][j] +
                         vars.Evec_Im[j] * diff_DA_Im[i][j];
        }
    }

    // rho = n^a n^b T_ab
    out.rho = 0.5 * msquared *
              (vars.Avec0_Re * vars.Avec0_Re + vars.Avec0_Im * vars.Avec0_Im);
    FOR2(i, j)
    {
        out.rho += 0.5 * metric_vars.gamma[i][j] *
                       (vars.Evec_Re[i] * vars.Evec_Re[j] +
                        vars.Evec_Im[i] * vars.Evec_Im[j]) +
                   0.5 * gamma_UU[i][j] * msquared *
                       (vars.Avec_Re[i] * vars.Avec_Re[j] +
                        vars.Avec_Im[i] * vars.Avec_Im[j]);

        FOR2(k, l)
        {
            out.rho += 0.5 * gamma_UU[i][k] * gamma_UU[j][l] *
                       (DA_Re[k][l] * diff_DA_Re[i][j] +
                        DA_Im[k][l] * diff_DA_Im[i][j]);
        }
    }

    return out;
}

// Adds VF evolution to the RHS

template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FixedBGComplexProcaField::matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // calculate full spatial christoffel symbols
    using namespace TensorAlgebra;
    const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
    const auto chris_phys = compute_christoffel(metric_vars.d1_gamma, gamma_UU);

    // evolution equations for vector fields A_0, A_i (note indices down) and
    // the conjugate momentum E^i (index up)

    // Start with real parts
    total_rhs.Avec0_Re = metric_vars.lapse * metric_vars.K * vars.Avec0_Re +
                         advec.Avec0_Re - metric_vars.lapse * vars.Zvec_Re;

    FOR2(i, j)
    {
        total_rhs.Avec0_Re +=
            -gamma_UU[i][j] * (vars.Avec_Re[i] * metric_vars.d1_lapse[j] +
                               metric_vars.lapse * d1.Avec_Re[i][j]);

        FOR1(k)
        {
            total_rhs.Avec0_Re += gamma_UU[i][j] * metric_vars.lapse *
                                  chris_phys.ULL[k][i][j] * vars.Avec_Re[k];
        }
    }

    FOR1(i)
    {
        total_rhs.Avec_Re[i] = -metric_vars.lapse * d1.Avec0_Re[i] -
                               vars.Avec0_Re * metric_vars.d1_lapse[i] +
                               advec.Avec_Re[i];
        FOR1(j)
        {
            total_rhs.Avec_Re[i] +=
                -metric_vars.lapse * metric_vars.gamma[i][j] * vars.Evec_Re[j] +
                vars.Avec_Re[j] * metric_vars.d1_shift[j][i];
        }
    }

    // variable for term (D_i A_j - D_j A_i)
    // NB Christoffel symbols cancel and take care with
    // indices - the second index is the derivative index
    Tensor<2, data_t> diff_DA_Re;
    FOR2(i, j) { diff_DA_Re[i][j] = d1.Avec_Re[j][i] - d1.Avec_Re[i][j]; }

    // NB This is for E^i with indices up
    FOR1(i)
    {
        total_rhs.Evec_Re[i] =
            metric_vars.lapse * metric_vars.K * vars.Evec_Re[i] +
            advec.Evec_Re[i];
        FOR1(j)
        {
            total_rhs.Evec_Re[i] +=
                gamma_UU[i][j] * (metric_vars.lapse * d1.Zvec_Re[j] +
                                  metric_vars.lapse * pow(m_vector_mass, 2.0) *
                                      vars.Avec_Re[j]) -
                vars.Evec_Re[j] * metric_vars.d1_shift[i][j];
        }

        FOR3(j, k, l)
        {
            total_rhs.Evec_Re[i] +=
                gamma_UU[k][j] * gamma_UU[i][l] *
                (metric_vars.d1_lapse[k] * diff_DA_Re[l][j] +
                 metric_vars.lapse *
                     (d2.Avec_Re[j][l][k] - d2.Avec_Re[l][j][k]));

            FOR1(m)
            {
                total_rhs.Evec_Re[i] +=
                    -gamma_UU[k][j] * gamma_UU[i][l] * metric_vars.lapse *
                    (chris_phys.ULL[m][k][l] * diff_DA_Re[m][j] +
                     chris_phys.ULL[m][k][j] * diff_DA_Re[l][m]);
            }
        }
    }

    // evolution equations for vector field Avec and its conjugate momentum Evec
    total_rhs.Zvec_Re =
        metric_vars.lapse * (pow(m_vector_mass, 2.0) * vars.Avec0_Re -
                             m_vector_damping * vars.Zvec_Re) +
        advec.Zvec_Re;

    FOR1(i)
    {
        total_rhs.Zvec_Re += metric_vars.lapse * d1.Evec_Re[i][i];
        FOR1(j)
        {
            total_rhs.Zvec_Re +=
                metric_vars.lapse * chris_phys.ULL[i][i][j] * vars.Evec_Re[j];
        }
    }

    // Now for the imaginary parts
    total_rhs.Avec0_Im = metric_vars.lapse * metric_vars.K * vars.Avec0_Im +
                         advec.Avec0_Im - metric_vars.lapse * vars.Zvec_Im;

    FOR2(i, j)
    {
        total_rhs.Avec0_Im +=
            -gamma_UU[i][j] * (vars.Avec_Im[i] * metric_vars.d1_lapse[j] +
                               metric_vars.lapse * d1.Avec_Im[i][j]);

        FOR1(k)
        {
            total_rhs.Avec0_Im += gamma_UU[i][j] * metric_vars.lapse *
                                  chris_phys.ULL[k][i][j] * vars.Avec_Im[k];
        }
    }

    FOR1(i)
    {
        total_rhs.Avec_Im[i] = -metric_vars.lapse * d1.Avec0_Im[i] -
                               vars.Avec0_Im * metric_vars.d1_lapse[i] +
                               advec.Avec_Im[i];
        FOR1(j)
        {
            total_rhs.Avec_Im[i] +=
                -metric_vars.lapse * metric_vars.gamma[i][j] * vars.Evec_Im[j] +
                vars.Avec_Im[j] * metric_vars.d1_shift[j][i];
        }
    }

    // variable for term (D_i A_j - D_j A_i)
    // NB Christoffel symbols cancel and take care with
    // indices - the second index is the derivative index
    Tensor<2, data_t> diff_DA_Im;
    FOR2(i, j) { diff_DA_Im[i][j] = d1.Avec_Im[j][i] - d1.Avec_Im[i][j]; }

    // NB This is for E^i with indices up
    FOR1(i)
    {
        total_rhs.Evec_Im[i] =
            metric_vars.lapse * metric_vars.K * vars.Evec_Im[i] +
            advec.Evec_Im[i];
        FOR1(j)
        {
            total_rhs.Evec_Im[i] +=
                gamma_UU[i][j] * (metric_vars.lapse * d1.Zvec_Im[j] +
                                  metric_vars.lapse * pow(m_vector_mass, 2.0) *
                                      vars.Avec_Im[j]) -
                vars.Evec_Im[j] * metric_vars.d1_shift[i][j];
        }

        FOR3(j, k, l)
        {
            total_rhs.Evec_Im[i] +=
                gamma_UU[k][j] * gamma_UU[i][l] *
                (metric_vars.d1_lapse[k] * diff_DA_Im[l][j] +
                 metric_vars.lapse *
                     (d2.Avec_Im[j][l][k] - d2.Avec_Im[l][j][k]));

            FOR1(m)
            {
                total_rhs.Evec_Im[i] +=
                    -gamma_UU[k][j] * gamma_UU[i][l] * metric_vars.lapse *
                    (chris_phys.ULL[m][k][l] * diff_DA_Im[m][j] +
                     chris_phys.ULL[m][k][j] * diff_DA_Im[l][m]);
            }
        }
    }

    // evolution equations for vector field Avec and its conjugate momentum Evec
    total_rhs.Zvec_Im =
        metric_vars.lapse * (pow(m_vector_mass, 2.0) * vars.Avec0_Im -
                             m_vector_damping * vars.Zvec_Im) +
        advec.Zvec_Im;

    FOR1(i)
    {
        total_rhs.Zvec_Im += metric_vars.lapse * d1.Evec_Im[i][i];
        FOR1(j)
        {
            total_rhs.Zvec_Im +=
                metric_vars.lapse * chris_phys.ULL[i][i][j] * vars.Evec_Im[j];
        }
    }
}

#endif /* FIXEDBGCOMPLEXPROCAFIELD_IMPL_HPP_ */
