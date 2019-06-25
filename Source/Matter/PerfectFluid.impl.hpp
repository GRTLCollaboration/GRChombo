/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(PERFECTFLUID_HPP_)
#error "This file should only be included through PerfectFluid.hpp"
#endif

#ifndef PERFECTFLUID_IMPL_HPP_
#define PERFECTFLUID_IMPL_HPP_
#include "DimensionDefinitions.hpp"

// Calculate the stress energy tensor elements
template <class eos_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> PerfectFluid<eos_t>::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;
    GeoVars<data_t> geo_vars;

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i, j)
    {
        out.Sij[i][j] =
          vars.density * vars.u[i] * vars.u[j] +
          vars.pressure * vars.h[i][j] / vars.chi;
    }

    // S_i (note lower index) = - n^a T_ai
    FOR1(i) { out.Si[i] =  geo_vars.lapse *                                               //FIXME:  if one uses lapse it doesn't work, I donno why
                        vars.density * vars.enthalpy * vars.u[i] * vars.u0; }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);


    // rho = n^a n^b T_ab
    out.rho =  geo_vars.lapse * geo_vars.lapse *
             (vars.density * vars.enthalpy * vars.u0 * vars.u0 -
              vars.pressure);

    return out;
}

// Adds in the RHS for the matter vars
template <class eos_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void PerfectFluid<eos_t>::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // the rhs vars
    FluidObject<data_t> rhs_fl;
    rhs_fl.D = 0;
    rhs_fl.E = 0;

    // advection terms
    FluidObject<data_t> advec_fl;

    // the vars
    Vars<data_t> vars_fl;
    vars_fl.W = vars.W;
    vars_fl.D = vars.D;
    vars_fl.E = vars.E;
    vars_fl.Z0 = vars.Z0;

    FOR1(i) {
      advec_fl.Z[i] = advec.Z[i];

      vars_fl.V[i] = vars.V[i];
      vars_fl.Z[i] = vars.Z[i];

      // rhs_fl.V[i] = 0;
      rhs_fl.Z[i] = 0;
    }

    // useful variable
    data_t dt_W = 0.0;                                                                  // FIXME: need to be defined

    {  // templated from (ScalarField)  matter_rhs_excl_potential                       // TODO: create indp function?
    /* ** starts braket */

    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // evolution equations for scalar field and (minus) its conjugate momentum
    FOR1(i)
    {
        rhs_fl.D += - d1.D[i] * vars_fl.V[i] - vars_fl.D * d1.V[i][i];
        rhs_fl.E += - d1.E[i] * vars_fl.V[i] - vars_fl.E * d1.V[i][i] +
                 ( - d1.W[i] * vars_fl.V[i] - vars_fl.W * d1.V[i][i] - dt_W)
                 * vars.pressure;

        rhs_fl.Z[i] +=  - vars_fl.Z0 * vars.lapse * d1.lapse[i] +
                    - vars.lapse * pow(vars.chi, 1.5) * d1.pressure[i]
                    +  advec_fl.Z[i];
    }

    FOR2(i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs_fl.Z[i] += - d1.V[j][j] * vars_fl.Z[i] - d1.Z[j][i] * vars_fl.V[j];


        FOR1(k)
        {

            FOR2(m,n)
            {
                rhs_fl.Z[i] +=  0.5 * h_UU[m][j] * h_UU[n][k]
                                    * vars_fl.Z[j] * vars_fl.Z[n]
                                    * d1.h[i][j][k] / vars_fl.Z0;
            }
        }
    }

    /* ** ends braket */
    }
}


template <class eos_t>
template <class data_t>
void PerfectFluid<eos_t>::compute(
  Cell<data_t> current_cell) const
{

    const auto vars = current_cell.template load_vars<Vars>();
    const auto geo_vars = current_cell.template load_vars<GeoVars>();
    auto up_vars = current_cell.template load_vars<Vars>();

    data_t pressure = 0.0;
    data_t enthalpy = 0.0;
    my_eos.compute_eos(pressure, enthalpy, vars);

    // Inverse metric
    const auto h_UU = TensorAlgebra::compute_inverse_sym(geo_vars.h);
    // data_t determinant = 1.0/vars.chi/vars.chi/vars.chi;

    //  useful vars
    data_t uiui = 0.0;                                                                  //FIXME   use data_t?
    data_t shift_ui = 0.0;
    data_t lapse2 = geo_vars.lapse * geo_vars.lapse;

    FOR2(i, j)
    {
      uiui += h_UU[i][j] * vars.u[i] * vars.u[j];
    }
    FOR1(i)
    {
      shift_ui += geo_vars.shift[i] * vars.u[i];     // See page 250 Shibata's book
    }
    shift_ui = shift_ui / lapse2;

    // fluit variables
    up_vars.density = vars.D / vars.W;
    up_vars.energy = vars.E / vars.D;
    up_vars.pressure = pressure;
    up_vars.enthalpy = enthalpy;

    FOR1(i) { up_vars.u[i] = vars.Z[i] / vars.D / up_vars.enthalpy;  }

    // from : -\lapse^2 * u_0^2 + (u_0 \shift^i u_i  + u_i h^{ij} u_l = -1
    up_vars.u0 = 0.5 * shift_ui +
                sqrt(0.25 * shift_ui * shift_ui + (uiui + 1) / lapse2);

    up_vars.W = geo_vars.lapse  * pow(geo_vars.chi, -3) / vars.u0;

    // Overwrite new values for fluid variables
    current_cell.store_vars(up_vars.density, c_density);
    current_cell.store_vars(up_vars.energy, c_energy);
    current_cell.store_vars(up_vars.pressure, c_pressure);
    current_cell.store_vars(up_vars.enthalpy, c_enthalpy);
    current_cell.store_vars(up_vars.u, GRInterval<c_u1, c_u3>());
    current_cell.store_vars(up_vars.u0, c_u0);
    current_cell.store_vars(up_vars.W, c_W);
}


#endif /* PERFECTFLUID_IMPL_HPP_ */
