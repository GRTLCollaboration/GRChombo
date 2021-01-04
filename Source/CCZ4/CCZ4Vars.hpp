/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CCZ4VARS_HPP_
#define CCZ4VARS_HPP_

#include "ADMConformalVars.hpp"
#include "BSSNVars.hpp"
#include "Tensor.hpp"
#include "VarsTools.hpp"

/// Namespace for CCZ4 vars
/** The structs in this namespace collect all the CCZ4 variables. It's main use
 *  is to make a local, nicely laid-out, copy of the CCZ4 variables for the
 *  current grid cell (Otherwise, this data would only exist on the grid in
 *  the huge, flattened Chombo array). \sa {CCZ4Vars, ADMConformalVars}
 **/
namespace CCZ4Vars
{
/// Vars object for CCZ4 vars, including gauge vars
template <class data_t>
struct VarsNoGauge : public BSSNVars::VarsNoGauge<data_t>
{
    data_t Theta; //!< CCZ4 quantity associated to Hamiltonian constraint

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        BSSNVars::VarsNoGauge<data_t>::enum_mapping(mapping_function);
        define_enum_mapping(mapping_function, c_Theta, Theta);
    }
};

/// Vars object for CCZ4 vars, including gauge vars
template <class data_t>
struct VarsWithGauge : public BSSNVars::VarsWithGauge<data_t>
{
    data_t Theta; //!< CCZ4 quantity associated to hamiltonian constraint

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        BSSNVars::VarsWithGauge<data_t>::enum_mapping(mapping_function);
        define_enum_mapping(mapping_function, c_Theta, Theta);
    }
};

/// Vars object for CCZ4 vars needing second derivs, excluding gauge vars
template <class data_t>
struct Diff2VarsNoGauge : public ADMConformalVars::Diff2VarsNoGauge<data_t>
{
};

/// Vars object for CCZ4 vars needing second derivs, including gauge vars
template <class data_t>
struct Diff2VarsWithGauge : public ADMConformalVars::Diff2VarsWithGauge<data_t>
{
};
} // namespace CCZ4Vars

#endif /* CCZ4VARS_HPP_ */
