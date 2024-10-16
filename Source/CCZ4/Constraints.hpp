/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This compute class calculates Hamiltonian and Momentum constraints

#ifndef CONSTRAINTS_HPP_
#define CONSTRAINTS_HPP_

#include "BSSNVars.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "simd.hpp"

#include "CCZ4Geometry.hpp"

#include <array>

class [[deprecated(
    "Use new Constraints class in NewConstraints.hpp")]] Constraints
{
  public:
    /// CCZ4 variables
    template <class data_t> using MetricVars = BSSNVars::VarsNoGauge<data_t>;

    /// CCZ4 variables
    template <class data_t>
    using Diff2Vars = BSSNVars::Diff2VarsNoGauge<data_t>;

    /// Vars object for Constraints
    template <class data_t> struct Vars
    {
        data_t Ham;
        Tensor<1, data_t> Mom;

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools;
            define_enum_mapping(mapping_function, c_Ham, Ham);
            define_enum_mapping(mapping_function, GRInterval<c_Mom1, c_Mom3>(),
                                Mom);
        }
    };

    Constraints(double dx, double cosmological_constant = 0);

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const FourthOrderDerivatives m_deriv;
    double m_cosmological_constant;

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    Vars<data_t>
    constraint_equations(const vars_t<data_t> &vars,
                         const vars_t<Tensor<1, data_t>> &d1,
                         const diff2_vars_t<Tensor<2, data_t>> &d2) const;
};

#include "Constraints.impl.hpp"

#endif /* CONSTRAINTS_HPP_ */
