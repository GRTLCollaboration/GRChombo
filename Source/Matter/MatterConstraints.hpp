/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MATTERCONSTRAINTS_HPP_
#define MATTERCONSTRAINTS_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Constraints.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "simd.hpp"
#include <array>

//!  Calculates the Hamiltonian and Momentum constraints with matter fields
/*!
     The class calculates the Hamiltonian and Momentum constraints at each point
   in a box. It inherits from the Constraints class which calculates the
   constraints without the matter terms. It adds in the matter terms for a given
   matter class matter_t, which must provide it with the Energy Momentum Tensor.
   For an example of a matter_t class see ScalarField. \sa Constraints(),
   ScalarField()
*/
template <class matter_t>
class [[deprecated("Use new MatterConstraints class in "
                   "NewMatterConstraints.hpp")]] MatterConstraints
    : public Constraints
{
  public:
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    // Inherit the variable definitions from CCZ4 + matter_t
    template <class data_t>
    struct BSSNMatterVars : public Constraints::MetricVars<data_t>,
                            public MatterVars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            Constraints::MetricVars<data_t>::enum_mapping(mapping_function);
            MatterVars<data_t>::enum_mapping(mapping_function);
        }
    };

    //!  Constructor of class MatterConstraints
    /*!
         Takes in the grid spacing, and matter object plus
         optionally the value of Newton's constant, which is set to one by
       default.
    */
    MatterConstraints(const matter_t a_matter, double dx,
                      double G_Newton = 1.0);

    //! The compute member which calculates the constraints at each point in the
    //! box
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    matter_t my_matter; //!< The matter object, e.g. a scalar field
    double m_G_Newton;  //!< Newton's constant, set to one by default.
};

#include "MatterConstraints.impl.hpp"

#endif /* MATTERCONSTRAINTS_HPP_ */
