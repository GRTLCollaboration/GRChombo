#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4UserVariables.hpp"
#include "DiagnosticVariables.hpp"

// TODO: This file can be auto-generated from a list of variable names
// Also, we should probably scope this enum too...
//
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

} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
