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
    c_rhoJ,
    c_Edot,
    c_Jdot,
    c_xMom,
    c_Stress,
    c_Source,
    c_BHMom,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi", "rho", "rhoJ", "Edot", "Jdot", "xMom", "Stress", "Source", "BHMom"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
