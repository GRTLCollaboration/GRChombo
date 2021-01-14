#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4UserVariables.hpp"

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

    c_phi = NUM_CCZ4_VARS,
    c_Pi,

    c_Rho,

    c_chi2,

    c_Weyl4_Re,
    c_Weyl4_Im,

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    user_variable_names = {

        "phi",      "Pi",

        "rho",

        "chi2",

        "Weyl4_Re", "Weyl4_Im"};

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, user_variable_names);
} // namespace UserVariables

#endif /* USERVARIABLES_HPP */
