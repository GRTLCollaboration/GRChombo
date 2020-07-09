/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MATTERONLY_HPP_
#define MATTERONLY_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//!  Calculates RHS of matter variables only, gravity vars assumed static
/*!
     The class calculates the RHS evolution for the matter variables.
     It does not assume a specific form of matter but is templated over a matter
   class matter_t. Please see the class SFMatter as an example of a matter_t.
     \sa CCZ4(), ScalarField(), VectorField()
*/

template <class matter_t> class MatterOnly
{
  public:
    //! Inherit the variable definitions from the Matter vars
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    template <class data_t> struct Vars : public MatterVars<data_t>
    {
        // BSSN vars needed in matter only rhs (ok for Proca and SF)
        data_t chi;
        Tensor<2, data_t> h;
        data_t K;
        data_t lapse;
        Tensor<1, data_t> shift;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            MatterVars<data_t>::enum_mapping(mapping_function);
            using namespace VarsTools;
            define_enum_mapping(mapping_function, c_K, K);
            define_enum_mapping(mapping_function, c_lapse, lapse);
            define_enum_mapping(mapping_function, c_chi, chi);
            define_enum_mapping(mapping_function,
                                GRInterval<c_shift1, c_shift3>(), shift);
            define_symmetric_enum_mapping(mapping_function,
                                          GRInterval<c_h11, c_h33>(), h);
        }
    };

    //! Inherit the variable definitions from the Matter vars
    //  Should not need any d2 for the metric vars (Proca and SF)
    template <class data_t>
    using MatterDiff2Vars = typename matter_t::template Diff2Vars<data_t>;

    template <class data_t> struct Diff2Vars : public MatterDiff2Vars<data_t>
    {
    };

    //!  Constructor of class MatterOnly
    MatterOnly(matter_t a_matter, double sigma, double dx)
        : my_matter(a_matter), m_sigma(sigma), m_deriv(dx), m_dx(dx)
    {
    }

    //!  The compute member which calculates the RHS at each point in the box
    //!  \sa matter_rhs_equation()
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variable and calculate
        // derivs
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
        const auto advec =
            m_deriv.template advection<Vars>(current_cell, vars.shift);

        // the RHS
        Vars<data_t> matter_rhs;
        VarsTools::assign(matter_rhs, 0.); // All components that are not
                                           // explicitly set in rhs_equation are
                                           // 0

        // add evolution of matter fields and dissipation
        my_matter.add_matter_rhs(matter_rhs, vars, d1, d2, advec);
        m_deriv.add_dissipation(matter_rhs, current_cell, m_sigma);

        // Write the rhs into the output FArrayBox
        current_cell.store_vars(matter_rhs);
    }

  protected:
    const matter_t my_matter;             //!< The matter object
    const FourthOrderDerivatives m_deriv; //!< An object for calculating
                                          //!< derivatives of the vars at the
                                          //!< point
    const double m_sigma;                 //!< Sigma for dissipation
    const double m_dx;                    //!< Grid spacing
};

#endif /* MATTERONLY_HPP_ */
