/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_sf,
    c_sfd,
    c_a,
    c_H,
    c_sf2,

    c_pot,
    c_kin,

    c_Ham,
    c_Mom,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Phi",
    "Pi",
    "ScaleFactor",
    "HubbleFactor",
    "PhiSq",
    "PotED",
    "KinED",
    "Ham",
    "Mom"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
