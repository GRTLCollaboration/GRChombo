/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "CCZ4UserVariables.hpp"
#include <array>
#include <string>

enum
{
    c_h = c_h11,
    c_A = c_A11,
    c_Gamma = c_Gamma1,
    c_shift = c_shift1,
    c_B = c_B1,

    NUM_VARS = NUM_CCZ4_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS> variable_names =
    ccz4_variable_names;
}

#endif /* USERVARIABLES_HPP */
