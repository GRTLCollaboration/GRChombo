/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(COMPLEXSCALARFIELD_HPP_)
#error "This file should only be included through ComplexScalarField.hpp"
#endif

#ifndef COMPLEXSCALARFIELD_IMPL_HPP_
#define COMPLEXSCALARFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> ComplexScalarField<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;

    // call the function which computes the em tensor excluding the potential
    emtensor_excl_potential(out, vars, d1, h_UU, chris_ULL);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi_Re = 0.0;
    data_t dVdphi_Im = 0.0;

    // compute potential and add constributions to EM Tensor
    my_potential.compute_potential(V_of_phi, dVdphi_Re, dVdphi_Im, vars);

    out.rho += V_of_phi;
    out.S += -3.0 * V_of_phi;
    FOR(i, j) { out.Sij[i][j] += -vars.h[i][j] * V_of_phi / vars.chi; }

    return out;
}

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void ComplexScalarField<potential_t>::emtensor_excl_potential(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<2, data_t> &h_UU,
    const Tensor<3, data_t> &chris_ULL)
{
    // Useful quantity Vt
    data_t Vt = -vars.Pi_Re * vars.Pi_Re - vars.Pi_Im * vars.Pi_Im;
    FOR2(i, j)
    {
        Vt += vars.chi * h_UU[i][j] *
              (d1.phi_Re[i] * d1.phi_Re[j] + d1.phi_Im[i] * d1.phi_Im[j]);
    }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR(i, j)
    {
        out.Sij[i][j] = -0.5 * vars.h[i][j] * Vt / vars.chi +
                        d1.phi_Re[i] * d1.phi_Re[j] +
                        d1.phi_Im[i] * d1.phi_Im[j];
    }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    // S_i (note lower index) = - n^a T_ai
    FOR(i)
    {
        out.Si[i] = -d1.phi_Re[i] * vars.Pi_Re - d1.phi_Im[i] * vars.Pi_Im;
    }

    // rho = n^a n^b T_ab
    out.rho = vars.Pi_Re * vars.Pi_Re + vars.Pi_Im * vars.Pi_Im + 0.5 * Vt;
}

// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void ComplexScalarField<potential_t>::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // first get the non potential part of the rhs
    // this may seem a bit long winded, but it makes the function
    // work for more multiple fields

    // call the function for the rhs excluding the potential
    matter_rhs_excl_potential(total_rhs, vars, d1, d2, advec);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi_Re = 0.0;
    data_t dVdphi_Im = 0.0;
    my_potential.compute_potential(V_of_phi, dVdphi_Re, dVdphi_Im, vars);

    // adjust RHS for the potential term
    total_rhs.Pi_Re += -vars.lapse * dVdphi_Re;
    total_rhs.Pi_Im += -vars.lapse * dVdphi_Im;
}

// the RHS excluding the potential terms
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void ComplexScalarField<potential_t>::matter_rhs_excl_potential(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec)
{
    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs.phi_Re = vars.lapse * vars.Pi_Re + advec.phi_Re;
    rhs.Pi_Re = vars.lapse * vars.K * vars.Pi_Re + advec.Pi_Re;
    rhs.phi_Im = vars.lapse * vars.Pi_Im + advec.phi_Im;
    rhs.Pi_Im = vars.lapse * vars.K * vars.Pi_Im + advec.Pi_Im;

    FOR(i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs.Pi_Re +=
            h_UU[i][j] * (-0.5 * d1.chi[j] * vars.lapse * d1.phi_Re[i] +
                          vars.chi * vars.lapse * d2.phi_Re[i][j] +
                          vars.chi * d1.lapse[i] * d1.phi_Re[j]);
        rhs.Pi_Im +=
            h_UU[i][j] * (-0.5 * d1.chi[j] * vars.lapse * d1.phi_Im[i] +
                          vars.chi * vars.lapse * d2.phi_Im[i][j] +
                          vars.chi * d1.lapse[i] * d1.phi_Im[j]);
        FOR(k)
        {
            rhs.Pi_Re += -vars.chi * vars.lapse * h_UU[i][j] *
                         chris.ULL[k][i][j] * d1.phi_Re[k];
            rhs.Pi_Im += -vars.chi * vars.lapse * h_UU[i][j] *
                         chris.ULL[k][i][j] * d1.phi_Im[k];
        }
    }
}

#endif /* COMPLEXSCALARFIELD_IMPL_HPP_ */
