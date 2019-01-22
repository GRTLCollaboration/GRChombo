/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

// assign enum to each variable
enum
{
    c_chi,
    c_phi,

    NUM_VARS
};

namespace UserVariables
{
static constexpr char const *variable_names[NUM_VARS] = {"chi", "phi"};
}

#endif /* USERVARIABLES_HPP */
