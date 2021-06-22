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
    c_Ham_abs_sum,

    c_Mom,
    c_Mom_abs_sum,

    c_rho1,
    c_rho2,
    c_flux1,
    c_flux2,
    c_source1,
    c_source2,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",  "Ham_abs_sum", "Mom",   "Mom_abs_sum", "rho1",
    "rho2", "flux1",       "flux2", "source1",     "source2"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
