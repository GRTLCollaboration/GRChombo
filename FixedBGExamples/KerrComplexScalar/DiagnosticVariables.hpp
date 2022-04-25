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
    c_rhoAngMom,
    c_rhoEnergy,
    c_fluxAngMom,
    c_fluxEnergy,
    c_sourceAngMom,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi",        "rhoAngMom",  "rhoEnergy",
    "fluxAngMom", "fluxEnergy", "sourceAngMom"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
