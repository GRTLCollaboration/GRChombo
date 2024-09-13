/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(INFLATIONDIAGNOSTICS_HPP_)
#error "This file should only be included through InflationDiagnostics.hpp"
#endif

#ifndef INFLATIONDIAGNOSTICS_IMPL_HPP_
#define INFLATIONDIAGNOSTICS_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class matter_t>
template <class data_t>
void InflationDiagnostics<matter_t>::compute(Cell<data_t> current_cell) const
{
    // Load local vars and calculate derivs
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Vars>(current_cell);

    // Inverse metric and Christoffel symbol
    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    // Define quantities
    data_t rho;
    data_t sqrt_gamma;
    data_t S;
    data_t rho_scaled;
    data_t S_scaled;
    data_t K_scaled;
    data_t a;

    // Energy Momentum Tensor
    const auto emtensor = m_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    sqrt_gamma = pow(vars.chi, -3. / 2.);
    rho = emtensor.rho;
    S = emtensor.S;
    K_scaled = vars.K / pow(vars.chi, 3. / 2.);
    rho_scaled = emtensor.rho / pow(vars.chi, 3. / 2.);
    S_scaled = emtensor.S / pow(vars.chi, 3. / 2.);

    // Write the constraints into the output FArrayBox
    current_cell.store_vars(sqrt_gamma, c_sqrt_gamma);
    current_cell.store_vars(rho, c_rho);
    current_cell.store_vars(rho_scaled, c_rho_scaled);
    current_cell.store_vars(S_scaled, c_S_scaled);
    current_cell.store_vars(K_scaled, c_K_scaled);

    // store_vars(out, current_cell);
}

#endif /* INFLATIONDIAGNOSTICS_IMPL_HPP_ */
