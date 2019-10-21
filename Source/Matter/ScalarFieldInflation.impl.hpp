/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SCALARFIELD_HPP_)
#error "This file should only be included through ScalarField.hpp"
#endif

#ifndef SCALARFIELD_IMPL_HPP_
#define SCALARFIELD_IMPL_HPP_
#include "DimensionDefinitions.hpp"

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> ScalarField<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;

    // Copy the field vars into SFObject
    SFObject<data_t> vars_sf;
    vars_sf.phi = vars.phi;
    vars_sf.Pi = vars.Pi;

    // call the function which computes the em tensor excluding the potential
    emtensor_excl_potential(out, vars, vars_sf, d1.phi, h_UU, chris_ULL);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi = 0.0;

    // compute potential and add constributions to EM Tensor
    my_potential.compute_potential(V_of_phi, dVdphi, vars);

    out.rho += V_of_phi;
    out.S += -3.0 * V_of_phi;
    FOR2(i, j) { out.Sij[i][j] += -vars.h[i][j] * V_of_phi / vars.chi; }

    return out;
}

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void ScalarField<potential_t>::emtensor_excl_potential(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars,
    const SFObject<data_t> &vars_sf, const Tensor<1, data_t> &d1_phi,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL)
{
    // Useful quantity Vt
    data_t Vt = -vars_sf.Pi * vars_sf.Pi;
    FOR2(i, j) { Vt += vars.chi * h_UU[i][j] * d1_phi[i] * d1_phi[j]; }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i, j)
    {
        out.Sij[i][j] =
            -0.5 * vars.h[i][j] * Vt / vars.chi + d1_phi[i] * d1_phi[j];
    }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    // S_i (note lower index) = - n^a T_ai
    FOR1(i) { out.Si[i] = -d1_phi[i] * vars_sf.Pi; }

    // rho = n^a n^b T_ab
    out.rho = vars_sf.Pi * vars_sf.Pi + 0.5 * Vt;
}

// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void ScalarField<potential_t>::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // first get the non potential part of the rhs
    // this may seem a bit long winded, but it makes the function
    // work for more multiple fields

    SFObject<data_t> rhs_sf;
    // advection terms
    SFObject<data_t> advec_sf;
    advec_sf.phi = advec.phi;
    advec_sf.Pi = advec.Pi;
    // the vars
    SFObject<data_t> vars_sf;
    vars_sf.phi = vars.phi;
    vars_sf.Pi = vars.Pi;

    // call the function for the rhs excluding the potential
    matter_rhs_excl_potential(rhs_sf, vars, vars_sf, d1, d1.phi, d2.phi,
                              advec_sf);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi = 0.0;

    // compute potential
    my_potential.compute_potential(V_of_phi, dVdphi, vars);

    // adjust RHS for the potential term
    total_rhs.phi = rhs_sf.phi;
    total_rhs.Pi = rhs_sf.Pi - vars.lapse * dVdphi;
}

// the RHS excluding the potential terms
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void ScalarField<potential_t>::matter_rhs_excl_potential(
    SFObject<data_t> &rhs_sf, const vars_t<data_t> &vars,
    const SFObject<data_t> &vars_sf, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<1, data_t> &d1_phi, const Tensor<2, data_t> &d2_phi,
    const SFObject<data_t> &advec_sf)
{
    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs_sf.phi = vars.lapse * vars_sf.Pi + advec_sf.phi;
    rhs_sf.Pi = vars.lapse * vars.K * vars_sf.Pi + advec_sf.Pi;

    FOR2(i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs_sf.Pi += h_UU[i][j] * (-0.5 * d1.chi[j] * vars.lapse * d1_phi[i] +
                                   vars.chi * vars.lapse * d2_phi[i][j] +
                                   vars.chi * d1.lapse[i] * d1_phi[j]);
        FOR1(k)
        {
            rhs_sf.Pi += -vars.chi * vars.lapse * h_UU[i][j] *
                         chris.ULL[k][i][j] * d1_phi[k];
        }
    }
}




template <class potential_t>
// template <class data_t>
template <class data_t, template <typename> class vars_t>
void ScalarField<potential_t>::compute(
  Cell<data_t> current_cell) const
{

    // Load local vars and calculate derivs
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Vars>(current_cell);

    auto up_vars = current_cell.template load_vars<Vars>();

    // Inverse metric and Christoffel symbol
    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    // Energy Momentum Tensor
    const auto emtensor = my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    data_t psi_dot = 0.0;
    data_t psi_dotdot = 0.0;
    // data_t t_rho = emtensor.rho;
    // data_t t_S = emtensor.S;
    up_vars.rho =  emtensor.rho;
    up_vars.S =  emtensor.S;


      FOR1(i)
      {
        V_i[i] = vars.Z[i] / (vars.E + vars.D + pressure_guess);
      }
      FOR2(i,j)
      {
        V2 +=  V_i[i] * h_UU[i][j] * V_i[j];
      }

      up_vars.W = 1.0 / sqrt(1.0 - V2);
      up_vars.density = vars.D / up_vars.W;
      up_vars.energy = (vars.E + vars.D * ( 1 - up_vars.W)
                       + pressure_guess * (1 - up_vars.W * up_vars.W))
                       / vars.D / up_vars.W;


      my_eos.compute_eos(pressure, enthalpy, vars);
      residual =  (pressure - pressure_guess);
      pressure += 0.5 * residual;


      criterion = simd_compare_gt(
              abs(residual), threshold_residual );

      condition = criterion;
      cont += 1;

    up_vars.pressure = pressure;
    up_vars.enthalpy = enthalpy;


    // FOR1(i) { up_vars.u[i] = V_i[i] * up_vars.W; }
    //
    // up_vars.u0 = up_vars.W / vars.newlapse;
    //
    // FOR1(i) { up_vars.V[i] = up_vars.u[i] / geo_vars.lapse / up_vars.u0
    //                           + geo_vars.shift[i] / geo_vars.lapse;  }


    // Overwrite new values for fluid variables
    current_cell.store_vars(up_vars.rho, c_rho);
    current_cell.store_vars(up_vars.S, c_S);
    current_cell.store_vars(up_vars.psi_dot, c_psid);
    current_cell.store_vars(up_vars.enthalpy, c_psidd);
    // current_cell.store_vars(up_vars.V, GRInterval<c_V1, c_V3>());
    // current_cell.store_vars(up_vars.u, GRInterval<c_u1, c_u3>());
}


#endif /* SCALARFIELD_IMPL_HPP_ */
