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
    c_rhoLinMom,
    c_rhoAngMom,
    c_rhoEnergy,
    c_fluxLinMom,
    c_fluxAngMom,
    c_fluxEnergy,
    c_sourceLinMom,
    c_sourceAngMom,
    c_sourceEnergy,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi",        "rhoLinMom",  "rhoAngMom",    "rhoEnergy",    "fluxLinMom",
    "fluxAngMom", "fluxEnergy", "sourceLinMom", "sourceAngMom", "sourceEnergy"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
