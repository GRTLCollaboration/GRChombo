/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_chi,     // the conformal factor of the metric, not evolved
    c_rhoE,    // the energy density of the SF
    c_rhoJ,    // the energy density of the SF
    c_rhoM,    // the energy density of the SF
    c_SourceE, // the energy density of the SF
    c_SourceM, // the energy density of the SF
    c_Edot,    // the energy density of the SF
    c_Jdot,    // the energy density of the SF
    c_Mdot,    // the energy density of the SF
    c_BHMom,   // the energy density of the SF

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi",     "rhoE", "rhoJ", "rhoM", "SourceE",
    "SourceM", "Edot", "Jdot", "Mdot", "BHMom"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
