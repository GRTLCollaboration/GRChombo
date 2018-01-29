/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This compute class calculates Hamiltonian and Momentum constraints

#ifndef CONSTRAINTS_HPP_
#define CONSTRAINTS_HPP_

#include "Cell.hpp"
#include "FArrayBox.H"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "simd.hpp"

#include "CCZ4Geometry.hpp"

#include <array>

class Constraints
{
  public:
    template <class data_t> struct Vars
    {
        data_t chi;
        Tensor<2, data_t> h;
        data_t K;
        Tensor<2, data_t> A;
        Tensor<1, data_t> Gamma;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function);
    };

    template <class data_t> struct constraints_t
    {
        data_t Ham;
        Tensor<1, data_t> Mom;
    };

    Constraints(double dx, double cosmological_constant = 0);

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const FourthOrderDerivatives m_deriv;
    double m_cosmological_constant;

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    constraints_t<data_t>
    constraint_equations(const vars_t<data_t> &vars,
                         const vars_t<Tensor<1, data_t>> &d1,
                         const diff2_vars_t<Tensor<2, data_t>> &d2) const;
};

#include "Constraints.impl.hpp"

#endif /* CONSTRAINTS_HPP_ */
