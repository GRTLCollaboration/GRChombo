/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "DiagnosticVariables.hpp"
#include <array>
#include <string>

// assign an enum to each variable
enum
{
    c_phi, // matter field added
    c_Pi,  //(minus) conjugate momentum

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS> variable_names = {"phi", "Pi"};
}

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
