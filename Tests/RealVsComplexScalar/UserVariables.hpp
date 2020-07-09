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

    c_phi,
    c_Pi,

    c_phi_Re,
    c_Pi_Re,

    c_phi_Im,
    c_Pi_Im,

    NUM_VARS
};

namespace UserVariables
{
static constexpr char const *variable_names[NUM_VARS] = {
    "chi",

    "h11",    "h12",    "h13",    "h22", "h23", "h33",

    "K",

    "A11",    "A12",    "A13",    "A22", "A23", "A33",

    "Theta",

    "Gamma1", "Gamma2", "Gamma3",

    "lapse",

    "shift1", "shift2", "shift3",

    "B1",     "B2",     "B3",

    "rho",    "S",

    "S1",     "S2",     "S3",

    "phi",    "Pi",

    "phi_Re", "Pi_Im",

    "phi_Im", "Pi_Im"};
}

#endif /* USERVARIABLES_HPP */
