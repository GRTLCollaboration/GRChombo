/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(USERVARIABLES_HPP)
#error "This file should only be included through UserVariables.hpp"
#endif

#ifndef USERVARIABLES_INC_HPP
#define USERVARIABLES_INC_HPP

#include "parstream.H"
#include <algorithm>
#include <array>
#include <string>

// This file must be included at the end of UserVariables.hpp

namespace UserVariables
{
/// Takes a string and returns the variable enum number if the string
/// matches one of those in UserVariables::variable_names, or returns -1
/// otherwise
static int variable_name_to_enum(const std::string &a_var_name)
{
    const auto var_name_it =
        std::find(variable_names.begin(), variable_names.end(), a_var_name);

    int var = std::distance(variable_names.begin(), var_name_it);
    if (var != NUM_VARS)
        return var;
    else
    {
        return -1;
    }
}

}; // namespace UserVariables

namespace DiagnosticVariables
{
/// Takes a string and returns the variable enum number if the string
/// matches one of those in UserVariables::variable_names, or returns -1
/// otherwise
static int variable_name_to_enum(const std::string &a_var_name)
{
    const auto var_name_it =
        std::find(variable_names.begin(), variable_names.end(), a_var_name);

    int var = std::distance(variable_names.begin(), var_name_it);
    if (var != NUM_DIAGNOSTIC_VARS)
        return var;
    else
    {
        return -1;
    }
}

}; // namespace DiagnosticVariables

#endif /* USERVARIABLES_INC_HPP_ */
