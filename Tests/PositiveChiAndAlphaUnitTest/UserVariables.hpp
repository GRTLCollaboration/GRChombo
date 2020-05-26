/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include <array>
#include <string>

enum
{
    c_chi,
    c_lapse,
    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS> variable_names = {
    "chi",
    "lapse",
};
}

#endif /* USERVARIABLES_HPP */
