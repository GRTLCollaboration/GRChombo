/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_Ham,
    c_Mom,
    c_Ham_abs_sum, // c_Ham_abs_sum
    c_Mom_abs_sum, // c_Mom_abs_sum
    c_rho,
    c_sqrt_gamma, // sqrt(gamma) = pow(chi,-3/2) volume factor of spatial metric
    c_rho_scaled,
    c_S_scaled,
    c_K_scaled,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "Mom",

    "Ham_abs",    "Mom_abs",  "rho",     "sqrt_gamma",
    "rho_scaled", "S_scaled", "K_scaled"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
