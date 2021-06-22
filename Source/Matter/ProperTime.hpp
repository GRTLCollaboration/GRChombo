/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PROPERTIME_HPP_
#define PROPERTIME_HPP_

#include "CCZ4.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the rhs for the ProperTime evolution based on the lapse
class ProperTime
{
  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        data_t lapse;
        data_t tau;
        Tensor<1, data_t> shift;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_lapse, lapse);
            VarsTools::define_enum_mapping(mapping_function, c_tau, tau);
            VarsTools::define_enum_mapping(
                mapping_function, GRInterval<c_shift1, c_shift3>(), shift);
        }
    };

  public:
    ProperTime(double a_dx) : m_deriv(a_dx) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables
        const auto vars = current_cell.template load_vars<Vars>();
        const auto advec =
            m_deriv.template advection<Vars>(current_cell, vars.shift);
        data_t tau = vars.lapse + advec.tau;

        // assign values of density in output box
        current_cell.store_vars(tau, c_tau);
    }
};

#endif /* PROPERTIME_HPP_ */
