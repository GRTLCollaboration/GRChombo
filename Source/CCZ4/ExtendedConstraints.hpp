/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This compute class calculates Hamiltonian and Momentum constraints

#ifndef EXT_CONSTRAINTS_HPP_
#define EXT_CONSTRAINTS_HPP_

#include "BSSNVars.hpp"
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
    /// CCZ4 variables
    template <class data_t> using Vars = BSSNVars::VarsNoGauge<data_t>;

    /// CCZ4 variables
    template <class data_t>
    using Diff2Vars = BSSNVars::Diff2VarsNoGauge<data_t>;

    template <class data_t> struct constraints_t
    {
        data_t Ham;
        // data_t Ham_K, Ham_rho, Ham_trA2, Ham_ricci;
        data_t rho, trA2, ricci_scalar, S, HamRel;
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

#include "ExtendedConstraints.impl.hpp"

#endif /* CONSTRAINTS_HPP_ */
