/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_phi,

    c_Avec1,
    c_Avec2,
    c_Avec3,

    c_Evec1,
    c_Evec2,
    c_Evec3,

    c_Z,

    c_chi,

    c_rho,
    c_rhoJ,

    c_gauss,

    NUM_VARS
};

namespace UserVariables
{
static constexpr char const *variable_names[NUM_VARS] = {
    "phi",

    "Avec1", "Avec2", "Avec3",

    "Evec1", "Evec2", "Evec3",

    "Z",

    "chi",

    "rho",   "rhoJ",

    "gauss"};

} // namespace UserVariables

#endif /* USERVARIABLES_HPP */
