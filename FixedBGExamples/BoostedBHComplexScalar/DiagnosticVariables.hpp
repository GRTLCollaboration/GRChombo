/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_chi,
    c_rho,
    c_Source,
    c_xMom,
    c_Mdot, // Momentum flux
    c_Edot, // Energy flux

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi", "rho", "Source", "xMom", "Mdot", "Edot"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
