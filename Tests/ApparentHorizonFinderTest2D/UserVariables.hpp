/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "EmptyDiagnosticVariables.hpp"
#include <array>
#include <string>

// assign enum to each variable
enum
{
    c_V,

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS> variable_names = {"V"};
}

#include "UserVariables.inc.hpp"

#include "AHStringGeometry.hpp"
#define AHSurfaceGeometry AHStringGeometry

#include "AHTest2DFunction.hpp"
#define AHFunction AHTest2DFunction

#endif /* USERVARIABLES_HPP */
