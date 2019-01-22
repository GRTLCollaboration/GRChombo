/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BSSNVARS_HPP_
#define BSSNVARS_HPP_

#include "ADMConformalVars.hpp"
#include "Tensor.hpp"
#include "VarsTools.hpp"

/// Namespace for BSSN vars
/** The structs in this namespace collect all the BSSN variables. It's main use
 *  is to make a local, nicely laid-out, copy of the BSSN variables for the
 *  current grid cell (Otherwise, this data would only exist on the grid in
 *  the huge, flattened Chombo array). \sa {CCZ4Vars, ADMConformalVars}
 **/
namespace BSSNVars
{
/// Vars object for BSSN vars excluding gauge vars
template <class data_t>
struct VarsNoGauge : public ADMConformalVars::VarsNoGauge<data_t>
{
    Tensor<1, data_t> Gamma; //!< Conformal connection functions

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        ADMConformalVars::VarsNoGauge<data_t>::enum_mapping(mapping_function);
        define_enum_mapping(mapping_function, GRInterval<c_Gamma1, c_Gamma3>(),
                            Gamma); //!< The auxilliary variable Gamma^i
    }
};

/// Vars object for BSSN vars, including gauge vars
template <class data_t> struct VarsWithGauge : public VarsNoGauge<data_t>
{
    data_t lapse;
    Tensor<1, data_t> shift;
    Tensor<1, data_t> B; //!< \f$B^i = \partial_t \beta^i\f$, this is used
                         //! for second order shift conditions

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        VarsNoGauge<data_t>::enum_mapping(mapping_function);
        define_enum_mapping(mapping_function, c_lapse, lapse);
        define_enum_mapping(mapping_function, GRInterval<c_shift1, c_shift3>(),
                            shift);
        define_enum_mapping(mapping_function, GRInterval<c_B1, c_B3>(), B);
    }
};

/// Vars object for BSSN vars needing second derivs, excluding gauge vars
template <class data_t>
struct Diff2VarsNoGauge : public ADMConformalVars::Diff2VarsNoGauge<data_t>
{
};

/// Vars object for BSSN vars needing second derivs, including gauge vars
template <class data_t>
struct Diff2VarsWithGauge : public ADMConformalVars::Diff2VarsWithGauge<data_t>
{
};
} // namespace BSSNVars

#endif /* BSSNVARS_HPP_ */
