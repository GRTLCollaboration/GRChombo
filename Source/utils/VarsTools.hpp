/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef VARSTOOLS_HPP_
#define VARSTOOLS_HPP_

// Chombo includes
#include "parstream.H" //Gives us pout()

// Our includes
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp"

// Chombo namespace
#include "UsingNamespace.H"

namespace VarsTools
{
template <typename mapping_function_t, typename data_t>
void define_enum_mapping(mapping_function_t mapping_function, const int &ivar,
                         data_t &scalar)
{
    mapping_function(ivar, scalar);
}

template <typename mapping_function_t, typename data_t, int start_var,
          int end_var>
void define_enum_mapping(mapping_function_t mapping_function,
                         const GRInterval<start_var, end_var> interval,
                         Tensor<1, data_t, end_var - start_var + 1> &tensor)
{
    for (int ivar = 0; ivar < interval.size(); ++ivar)
        mapping_function(start_var + ivar, tensor[ivar]);
}

template <typename mapping_function_t, typename data_t, int start_var,
          int end_var>
void define_symmetric_enum_mapping(
    mapping_function_t mapping_function,
    const GRInterval<start_var, end_var> interval, Tensor<2, data_t> &tensor)
{
    static_assert(interval.size() ==
                      DEFAULT_TENSOR_DIM * (DEFAULT_TENSOR_DIM + 1) / 2,
                  "Interval has wrong size");
    int idx = 0;
    FOR(idir1)
    {
        for (int idir2 = idir1; idir2 < DEFAULT_TENSOR_DIM; ++idir2, ++idx)
        {
            mapping_function(start_var + idx, tensor[idir1][idir2]);
            if (idir1 != idir2)
                mapping_function(start_var + idx, tensor[idir2][idir1]);
        }
    }
}

//--> Begin: Helper for the assign function
template <class nested_template> struct strip_nested_template;

template <template <typename> class outermost_layer, class inner_part>
struct strip_nested_template<outermost_layer<inner_part>>
{
    using type = inner_part;
};
//<-- End: Helper for the assign function

/// Writes data directly into all variables
/**if this variables has multiple components (e.g. if it is an array of
 *derivatives) the data can be written directly into these components by
 *specifying an arbitrary number of icomps
 */
template <class vars_t, typename value_t>
ALWAYS_INLINE void assign(vars_t &vars, const value_t &value)
{
    // The template magic below is needed to make sure that we can write
    // assign(vars, 0.)  and 0. gets correctly cast from double to simd<double>
    // if necessary.
    using data_t = typename strip_nested_template<vars_t>::type;
    vars.enum_mapping([&value](const int &ivar, data_t &var)
                      { var = static_cast<data_t>(value); });
}

/// Prints all elements of the vars element with component names
/// (Very useful for debugging)
template <template <typename> class vars_t, typename data_t>
void print(const vars_t<data_t> &vars)
{
    vars.enum_mapping(
        [](const int &ivar, data_t &var) {
            pout() << UserVariables::variable_names[ivar] << ": " << var
                   << "\n";
        });
}
} // namespace VarsTools

#endif /* VARSTOOLS_HPP_ */
