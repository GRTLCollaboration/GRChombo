/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MATTERCCZ4_HPP_
#define MATTERCCZ4_HPP_

#include "CCZ4.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//!  Calculates RHS using CCZ4 including matter terms, and matter variable
//!  evolution
/*!
     The class calculates the RHS evolution for all the variables. It inherits
   from the CCZ4 class, which it uses to do the non matter evolution of
   variables. It then adds in the additional matter terms to the CCZ4 evolution
   (those including the stress energy tensor), and calculates the evolution of
   the matter variables. It does not assume a specific form of matter but is
   templated over a matter class matter_t. Please see the class ScalarField as
   an example of a matter_t. \sa CCZ4(), ScalarField()
*/

template <class matter_t> class MatterCCZ4 : public CCZ4
{
  public:
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    template <class data_t>
    using MatterDiff2Vars = typename matter_t::template Diff2Vars<data_t>;

    // Inherit the variable definitions from CCZ4 + matter_t
    template <class data_t>
    struct Vars : public CCZ4::Vars<data_t>, public MatterVars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            CCZ4::Vars<data_t>::enum_mapping(mapping_function);
            MatterVars<data_t>::enum_mapping(mapping_function);
        }
    };

    template <class data_t>
    struct Diff2Vars : public CCZ4::Diff2Vars<data_t>,
                       public MatterDiff2Vars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            CCZ4::Diff2Vars<data_t>::enum_mapping(mapping_function);
            MatterDiff2Vars<data_t>::enum_mapping(mapping_function);
        }
    };

    //!  Constructor of class MatterCCZ4
    /*!
       Inputs are the grid spacing, plus the CCZ4 evolution parameters and a
       matter object. It also takes the dissipation parameter sigma, and allows
       the formulation to be toggled between CCZ4 and BSSN. The default is CCZ4.
       It allows the user to set the value of Newton's constant, which is set to
       one by default.
    */
    MatterCCZ4(matter_t a_matter, params_t params, double dx, double sigma,
               int formulation = CCZ4::USE_CCZ4, double G_Newton = 1.0);

    //!  The compute member which calculates the RHS at each point in the box
    //!  \sa matter_rhs_equation()
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    //! The function which adds in the EM Tensor terms to the CCZ4 rhs \sa
    //! compute()
    template <class data_t>
    void add_emtensor_rhs(
        Vars<data_t>
            &matter_rhs, //!< the RHS data for each variable at that point.
        const Vars<data_t> &vars, //!< the value of the variables at the point.
        const Vars<Tensor<1, data_t>>
            &d1 //!< the value of the first derivatives of the variables.
    ) const;

    // Class members
    matter_t my_matter;      //!< The matter object, e.g. a scalar field.
    const double m_G_Newton; //!< Newton's constant, set to one by default.
};

#include "MatterCCZ4.impl.hpp"

#endif /* MATTERCCZ4_HPP_ */
