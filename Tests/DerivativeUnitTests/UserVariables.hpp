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
    c_d1,
    c_d2,
    c_d2_mixed,
    c_diss,
    c_advec_up,
    c_advec_down,
    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS> variable_names = {
    "d1", "d2", "d2_mixed", "diss", "advec_up", "advec_down"};
}

#endif /* USERVARIABLES_HPP */
