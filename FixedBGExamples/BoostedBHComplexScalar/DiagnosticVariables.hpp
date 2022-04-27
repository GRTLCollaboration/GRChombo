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
    c_rhoEnergy,
    c_fluxEnergy,
    c_rhoLinMom,
    c_fluxLinMom,
    c_sourceLinMom,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi",       "rhoEnergy",  "fluxEnergy",
    "rhoLinMom", "fluxLinMom", "sourceLinMom"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
