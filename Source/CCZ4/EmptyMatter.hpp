/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EMPTYMATTER_HPP_
#define EMPTYMATTER_HPP_

#include "CCZ4Geometry.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS, total num of components
#include "VarsTools.hpp"

//!  Calculates the matter type specific elements such as the EMTensor and
//   matter evolution
/*!
     This class is used when there is no matter but a matter template is needed
     e.g. in relaxation or in the Sommerfeld BCs.
*/
class EmptyMatter
{
  public:
    //!  Constructor of class EmptyMatter
    EmptyMatter() {}

    template <class data_t> struct Vars
    {
        // no vars
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
        }
    };

    template <class data_t> struct Diff2Vars
    {
        // no vars
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
        }
    };

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, including the potential
    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars,          //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const Tensor<2, data_t> &h_UU, //!< the inverse metric (raised indices)
        const Tensor<3, data_t>
            &chris_ULL) //!< the conformal christoffel symbol
        const
    {
        // all components zero
        emtensor_t<data_t> out;
        out.rho = 0.0;
        out.S = 0.0;
        FOR1(i)
        {
            out.Si[i] = 0.0;
            FOR1(j) { out.Sij[i][j] = 0.0; }
        }
        return out;
    }

    //! The function which adds in the RHS for the matter field vars,
    //! including the potential
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void add_matter_rhs(
        rhs_vars_t<data_t> &total_rhs,       //!< value of the RHS for all vars
        const vars_t<data_t> &vars,          //!< value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< value of the 2nd derivs
        const vars_t<data_t> &advec) //!< the value of the advection terms
        const
    {
    }

    //! The function which adds in the Sommerfeld RHS at the boundary for
    //! the matter field vars \sa SommerfeldBoundary() note this function
    //! will be called by SommerfeldBoundary which does not support simd
    //!  Currently the asymptotic vals are zero
    template <class vars_t>
    static void add_asymptotic_vals(vars_t &boundary_rhs, const double radius,
                                    const double time)
    {
    }
};

#endif /* EMPTYMATTER_HPP_ */
