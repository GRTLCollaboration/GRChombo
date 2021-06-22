/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_chi,

    c_h11,
    c_h12,
    c_h13,
    c_h22,
    c_h23,
    c_h33,

    c_K,

    c_A11,
    c_A12,
    c_A13,
    c_A22,
    c_A23,
    c_A33,

    c_Theta,

    c_Gamma1,
    c_Gamma2,
    c_Gamma3,

    c_lapse,

    c_shift1,
    c_shift2,
    c_shift3,

    c_B1,
    c_B2,
    c_B3,

    c_rho, // Energy momentum comps
    c_S1,
    c_S2,
    c_S3,
    c_S11,
    c_S12,
    c_S13,
    c_S22,
    c_S23,
    c_S33,

    c_Ham,
    c_Weyl4_Re,
    c_Weyl4_Im,

    c_Mom1,
    c_Mom2,
    c_Mom3,

    NUM_VARS
};

namespace UserVariables
{
static constexpr char const *variable_names[NUM_VARS] = {
    "chi",

    "h11",      "h12",      "h13",    "h22", "h23", "h33",

    "K",

    "A11",      "A12",      "A13",    "A22", "A23", "A33",

    "Theta",

    "Gamma1",   "Gamma2",   "Gamma3",

    "lapse",

    "shift1",   "shift2",   "shift3",

    "B1",       "B2",       "B3",

    "rho",      "S1",       "S2",     "S3",

    "S11",      "S12",      "S13",    "S22", "S23", "S33",

    "Ham",

    "Weyl4_Re", "Weyl4_Im",

    "Mom1",     "Mom2",     "Mom3"};
}

#endif /* USERVARIABLES_HPP */
