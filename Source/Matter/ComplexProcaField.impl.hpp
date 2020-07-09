/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(COMPLEXPROCAFIELD_HPP_)
#error "This file should only be included through ComplexProcaField.hpp"
#endif

#ifndef COMPLEXPROCAFIELD_IMPL_HPP_
#define COMPLEXPROCAFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> ComplexProcaField::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;

    // Some useful quantities
    const double msquared = pow(m_vector_mass, 2.0);
    const data_t one_over_chi2 = pow(vars.chi, -2.0);

    // calculate full spatial christoffel symbols
    const Tensor<3, data_t> chris_phys = TensorAlgebra::compute_phys_chris(
        d1.chi, vars.chi, vars.h, h_UU, chris_ULL);

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
            DA_Re[i][j] += -chris_phys[k][i][j] * vars.Avec_Re[k];
            DA_Im[i][j] += -chris_phys[k][i][j] * vars.Avec_Im[k];
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
                        0.5 * vars.h[i][j] / vars.chi * vars.Avec0_Re *
                            vars.Avec0_Re) +
            // Im part
            msquared *
                (vars.Avec_Im[i] * vars.Avec_Im[j] +
                 0.5 * vars.h[i][j] / vars.chi * vars.Avec0_Im * vars.Avec0_Im);

        FOR2(k, l)
        {
            out.Sij[i][j] +=
                // Re part
                vars.chi * h_UU[k][l] * (diff_DA_Re[i][k] * diff_DA_Re[j][l]) -
                vars.h[i][k] * vars.h[j][l] * one_over_chi2 * vars.Evec_Re[k] *
                    vars.Evec_Re[l] +
                0.5 * vars.h[k][l] * vars.h[i][j] * one_over_chi2 *
                    vars.Evec_Re[k] * vars.Evec_Re[l] -
                0.5 * h_UU[k][l] * vars.h[i][j] * msquared * vars.Avec_Re[k] *
                    vars.Avec_Re[l] +
                // Im part
                vars.chi * h_UU[k][l] * (diff_DA_Im[i][k] * diff_DA_Im[j][l]) -
                vars.h[i][k] * vars.h[j][l] * one_over_chi2 * vars.Evec_Im[k] *
                    vars.Evec_Im[l] +
                0.5 * vars.h[k][l] * vars.h[i][j] * one_over_chi2 *
                    vars.Evec_Im[k] * vars.Evec_Im[l] -
                0.5 * h_UU[k][l] * vars.h[i][j] * msquared * vars.Avec_Im[k] *
                    vars.Avec_Im[l];

            FOR2(m, n)
            {
                out.Sij[i][j] += -0.5 * vars.chi * vars.h[i][j] * h_UU[k][m] *
                                 h_UU[l][n] *
                                 (DA_Re[k][l] * diff_DA_Re[m][n] +
                                  DA_Im[k][l] * diff_DA_Im[m][n]);
            }
        }
    }

    // S = Tr_S_ij
    out.S = 0.0;
    FOR2(i, j) { out.S += out.Sij[i][j] * vars.chi * h_UU[i][j]; }

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
        out.rho += 0.5 * vars.h[i][j] / vars.chi *
                       (vars.Evec_Re[i] * vars.Evec_Re[j] +
                        vars.Evec_Im[i] * vars.Evec_Im[j]) +
                   0.5 * vars.chi * h_UU[i][j] * msquared *
                       (vars.Avec_Re[i] * vars.Avec_Re[j] +
                        vars.Avec_Im[i] * vars.Avec_Im[j]);

        FOR2(k, l)
        {
            out.rho += 0.5 * vars.chi * vars.chi * h_UU[i][k] * h_UU[j][l] *
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
void ComplexProcaField::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // calculate full spatial christoffel symbols
    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);
    const Tensor<3, data_t> chris_phys =
        compute_phys_chris(d1.chi, vars.chi, vars.h, h_UU, chris.ULL);

    // evolution equations for vector fields A_0, A_i (note indices down) and
    // the conjugate momentum E^i (index up)

    // Start with real parts
    total_rhs.Avec0_Re = vars.lapse * vars.K * vars.Avec0_Re + advec.Avec0_Re -
                         vars.lapse * vars.Zvec_Re;

    FOR2(i, j)
    {
        total_rhs.Avec0_Re +=
            -vars.chi * h_UU[i][j] *
            (vars.Avec_Re[i] * d1.lapse[j] + vars.lapse * d1.Avec_Re[i][j]);

        FOR1(k)
        {
            total_rhs.Avec0_Re += vars.chi * h_UU[i][j] * vars.lapse *
                                  chris_phys[k][i][j] * vars.Avec_Re[k];
        }
    }

    FOR1(i)
    {
        total_rhs.Avec_Re[i] = -vars.lapse * d1.Avec0_Re[i] -
                               vars.Avec0_Re * d1.lapse[i] + advec.Avec_Re[i];
        FOR1(j)
        {
            total_rhs.Avec_Re[i] +=
                -vars.lapse * vars.h[i][j] * vars.Evec_Re[j] / vars.chi +
                vars.Avec_Re[j] * d1.shift[j][i];
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
            vars.lapse * vars.K * vars.Evec_Re[i] + advec.Evec_Re[i];
        FOR1(j)
        {
            total_rhs.Evec_Re[i] +=
                vars.chi * h_UU[i][j] *
                    (vars.lapse * d1.Zvec_Re[j] +
                     vars.lapse * pow(m_vector_mass, 2.0) * vars.Avec_Re[j]) -
                vars.Evec_Re[j] * d1.shift[i][j];
        }

        FOR3(j, k, l)
        {
            total_rhs.Evec_Re[i] +=
                vars.chi * vars.chi * h_UU[k][j] * h_UU[i][l] *
                (d1.lapse[k] * diff_DA_Re[l][j] +
                 vars.lapse * (d2.Avec_Re[j][l][k] - d2.Avec_Re[l][j][k]));

            FOR1(m)
            {
                total_rhs.Evec_Re[i] +=
                    -vars.chi * vars.chi * h_UU[k][j] * h_UU[i][l] *
                    vars.lapse *
                    (chris_phys[m][k][l] * diff_DA_Re[m][j] +
                     chris_phys[m][k][j] * diff_DA_Re[l][m]);
            }
        }
    }

    // evolution equations for vector field Avec and its conjugate momentum Evec
    total_rhs.Zvec_Re = vars.lapse * (pow(m_vector_mass, 2.0) * vars.Avec0_Re -
                                      m_vector_damping * vars.Zvec_Re) +
                        advec.Zvec_Re;

    FOR1(i)
    {
        total_rhs.Zvec_Re += vars.lapse * d1.Evec_Re[i][i];
        FOR1(j)
        {
            total_rhs.Zvec_Re +=
                vars.lapse * chris_phys[i][i][j] * vars.Evec_Re[j];
        }
    }

    // Now for the imaginary parts
    total_rhs.Avec0_Im = vars.lapse * vars.K * vars.Avec0_Im + advec.Avec0_Im -
                         vars.lapse * vars.Zvec_Im;

    FOR2(i, j)
    {
        total_rhs.Avec0_Im +=
            -vars.chi * h_UU[i][j] *
            (vars.Avec_Im[i] * d1.lapse[j] + vars.lapse * d1.Avec_Im[i][j]);

        FOR1(k)
        {
            total_rhs.Avec0_Im += vars.chi * h_UU[i][j] * vars.lapse *
                                  chris_phys[k][i][j] * vars.Avec_Im[k];
        }
    }

    FOR1(i)
    {
        total_rhs.Avec_Im[i] = -vars.lapse * d1.Avec0_Im[i] -
                               vars.Avec0_Im * d1.lapse[i] + advec.Avec_Im[i];
        FOR1(j)
        {
            total_rhs.Avec_Im[i] +=
                -vars.lapse * vars.h[i][j] * vars.Evec_Im[j] / vars.chi +
                vars.Avec_Im[j] * d1.shift[j][i];
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
            vars.lapse * vars.K * vars.Evec_Im[i] + advec.Evec_Im[i];
        FOR1(j)
        {
            total_rhs.Evec_Im[i] +=
                vars.chi * h_UU[i][j] *
                    (vars.lapse * d1.Zvec_Im[j] +
                     vars.lapse * pow(m_vector_mass, 2.0) * vars.Avec_Im[j]) -
                vars.Evec_Im[j] * d1.shift[i][j];
        }

        FOR3(j, k, l)
        {
            total_rhs.Evec_Im[i] +=
                vars.chi * vars.chi * h_UU[k][j] * h_UU[i][l] *
                (d1.lapse[k] * diff_DA_Im[l][j] +
                 vars.lapse * (d2.Avec_Im[j][l][k] - d2.Avec_Im[l][j][k]));

            FOR1(m)
            {
                total_rhs.Evec_Im[i] +=
                    -vars.chi * vars.chi * h_UU[k][j] * h_UU[i][l] *
                    vars.lapse *
                    (chris_phys[m][k][l] * diff_DA_Im[m][j] +
                     chris_phys[m][k][j] * diff_DA_Im[l][m]);
            }
        }
    }

    // evolution equations for vector field Avec and its conjugate momentum Evec
    total_rhs.Zvec_Im = vars.lapse * (pow(m_vector_mass, 2.0) * vars.Avec0_Im -
                                      m_vector_damping * vars.Zvec_Im) +
                        advec.Zvec_Im;

    FOR1(i)
    {
        total_rhs.Zvec_Im += vars.lapse * d1.Evec_Im[i][i];
        FOR1(j)
        {
            total_rhs.Zvec_Im +=
                vars.lapse * chris_phys[i][i][j] * vars.Evec_Im[j];
        }
    }
}

#endif /* COMPLEXPROCAFIELD_IMPL_HPP_ */
