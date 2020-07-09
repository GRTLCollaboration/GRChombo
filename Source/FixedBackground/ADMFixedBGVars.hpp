/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ADMFIXEDBGVARS_HPP_
#define ADMFIXEDBGVARS_HPP_

#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"

/// Namespace for ADM vars in conformally decomposed form
/** The structs in this namespace collect all the ADM variables. It's main use
 *is to make a local, nicely laid-out, copy of the ADM variables for the
 *current grid cell (Otherwise, this data would only exist on the grid in the
 *huge, flattened Chombo array). \sa {CCZ4Vars, BSSNVars}
 **/
namespace ADMFixedBGVars
{
/// Vars object for ADM vars used in FixedBG evolution
template <class data_t> struct Vars
{
    // ADM vars needed in matter only rhs (ok for Proca and SF)
    Tensor<2, data_t> gamma;
    Tensor<2, data_t> K_tensor;
    data_t K;
    data_t lapse;
    Tensor<1, data_t> shift;
    Tensor<2, Tensor<1, data_t>> d1_gamma;
    Tensor<1, data_t> d1_lapse;
    Tensor<2, data_t> d1_shift;
};

} // namespace ADMFixedBGVars

#endif /* ADMFIXEDBGVARS_HPP_ */
