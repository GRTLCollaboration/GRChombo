/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP_
#define USERVARIABLES_HPP_

enum
{
    c_chi,

    c_h,
    c_h11 = c_h,
    c_h12,
    c_h13,
    c_h22,
    c_h23,
    c_h33,

    c_K,

    c_A,
    c_A11 = c_A,
    c_A12,
    c_A13,
    c_A22,
    c_A23,
    c_A33,

    c_Theta,

    c_Gamma,
    c_Gamma1 = c_Gamma,
    c_Gamma2,
    c_Gamma3,

    c_lapse,

    c_shift,
    c_shift1 = c_shift,
    c_shift2,
    c_shift3,

    c_B,
    c_B1 = c_B,
    c_B2,
    c_B3,

    c_rho,
    c_S,
    c_S1,
    c_S2,
    c_S3,

    c_Avec0,

    c_Avec1,
    c_Avec2,
    c_Avec3,

    c_Evec1,
    c_Evec2,
    c_Evec3,

    c_Zvec,

    c_Avec0_Re,

    c_Avec1_Re,
    c_Avec2_Re,
    c_Avec3_Re,

    c_Evec1_Re,
    c_Evec2_Re,
    c_Evec3_Re,

    c_Zvec_Re,

    c_Avec0_Im,

    c_Avec1_Im,
    c_Avec2_Im,
    c_Avec3_Im,

    c_Evec1_Im,
    c_Evec2_Im,
    c_Evec3_Im,

    c_Zvec_Im,

    NUM_VARS
};

namespace UserVariables
{
static constexpr char const *variable_names[NUM_VARS] = {
    "chi",

    "h11",      "h12",      "h13",      "h22", "h23", "h33",

    "K",

    "A11",      "A12",      "A13",      "A22", "A23", "A33",

    "Theta",

    "Gamma1",   "Gamma2",   "Gamma3",

    "lapse",

    "shift1",   "shift2",   "shift3",

    "B1",       "B2",       "B3",

    "rho",      "S",

    "S1",       "S2",       "S3",

    "Avec0",

    "Avec1",    "Avec2",    "Avec3",

    "Evec1",    "Evec2",    "Evec3",

    "Zvec",

    "Avec0_Re",

    "Avec1_Re", "Avec2_Re", "Avec3_Re",

    "Evec1_Re", "Evec2_Re", "Evec3_Re",

    "Zvec_Re",

    "Avec0_Im",

    "Avec1_Im", "Avec2_Im", "Avec3_Im",

    "Evec1_Im", "Evec2_Im", "Evec3_Im",

    "Zvec_Im"};
}

#endif /* USERVARIABLES_HPP */
