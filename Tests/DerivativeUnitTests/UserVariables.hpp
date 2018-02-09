/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP_
#define USERVARIABLES_HPP_

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
static constexpr char const *variable_names[NUM_VARS] = {
    "d1", "d2", "d2_mixed", "diss", "advec_up", "advec_down"};
}

#endif /* USERVARIABLES_HPP_ */
