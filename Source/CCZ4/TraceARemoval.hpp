/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This class enforces A to be trace-free
#ifndef TRACEAREMOVAL_HPP_
#define TRACEAREMOVAL_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"

class TraceARemoval
{
  public:
    template <class data_t> struct Vars
    {
        Tensor<2, data_t> h;
        Tensor<2, data_t> A;

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function);
    };

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto vars = current_cell.template load_vars<Vars>();

        const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
        TensorAlgebra::make_trace_free(vars.A, vars.h, h_UU);

        current_cell.store_vars(vars);
    }
};

template <class data_t>
template <typename mapping_function_t>
void TraceARemoval::Vars<data_t>::enum_mapping(
    mapping_function_t mapping_function)
{
    VarsTools::define_symmetric_enum_mapping(mapping_function,
                                             GRInterval<c_h11, c_h33>(), h);
    VarsTools::define_symmetric_enum_mapping(mapping_function,
                                             GRInterval<c_A11, c_A33>(), A);
}

#endif /* TRACEAREMOVAL_HPP_ */
