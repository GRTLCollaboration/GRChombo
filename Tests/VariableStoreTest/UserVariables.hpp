/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP_
#define USERVARIABLES_HPP_

enum
{
    c_var,
    c_sym_var,
    NUM_VARS
};

namespace UserVariables
{
static constexpr char const *variable_names[NUM_VARS] = {
    "var", "sym_var",
};
}

#endif /* USERVARIABLES_HPP_ */
