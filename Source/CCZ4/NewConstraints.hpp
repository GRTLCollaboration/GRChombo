/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This compute class calculates Hamiltonian and Momentum constraints

#ifndef NEWCONSTRAINTS_HPP_
#define NEWCONSTRAINTS_HPP_

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
    template <class data_t> using MetricVars = BSSNVars::VarsNoGauge<data_t>;

    /// CCZ4 variables
    template <class data_t>
    using Diff2Vars = BSSNVars::Diff2VarsNoGauge<data_t>;

    /// Vars object for Constraints
    template <class data_t> struct Vars
    {
        data_t Ham;
        data_t Ham_abs_terms;
        Tensor<1, data_t> Mom;
        Tensor<1, data_t> Mom_abs_terms;
    };

    // Constructor which allows specifying Ham and Mom vars
    // if the interval of a_c_Moms is of size 1, then
    // sqrt(Mom1^2 + Mom2^2 + Mom3^2) is stored in that variable
    // ...abs_terms stores the absolute value of the individual terms in the
    // conformally decomposed expressions which can be used in to normalize
    // the constraint violations
    // Any zero-length Interval or negative var is not calculated
    Constraints(double dx, int a_c_Ham, const Interval &a_c_Moms,
                int a_c_Ham_abs_terms = -1,
                const Interval &a_c_Moms_abs_terms = Interval(),
                double cosmological_constant = 0.0);

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const FourthOrderDerivatives m_deriv;
    const int m_c_Ham;
    const Interval m_c_Moms;
    const int m_c_Ham_abs_terms = -1;
    const Interval m_c_Moms_abs_terms;
    double m_cosmological_constant;

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    Vars<data_t> constraint_equations(const vars_t<data_t> &vars,
                                      const vars_t<Tensor<1, data_t>> &d1,
                                      const diff2_vars_t<Tensor<2, data_t>> &d2,
                                      const Tensor<2, data_t> &h_UU,
                                      const chris_t<data_t> &chris) const;

    template <class data_t>
    void store_vars(Vars<data_t> &out, Cell<data_t> &current_cell) const;
};

#include "NewConstraints.impl.hpp"

#endif /* NEWCONSTRAINTS_HPP_ */
