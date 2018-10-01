/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP_
#define USERVARIABLES_HPP_

// assigns number to each variable
enum
{
    c_phi,

    NUM_VARS
};

namespace UserVariables
{
static constexpr char const *variable_names[NUM_VARS] = {
    "phi"};
}

#endif /* USERVARIABLES_HPP_ */
