/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGPROCAFIELDTEST_HPP_
#define FIXEDBGPROCAFIELDTEST_HPP_

#include "ADMFixedBGVars.hpp"
#include "CCZ4Geometry.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Potential.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS
#include "VarsTools.hpp"
#include "simd.hpp"

//!  Calculates the matter type specific elements such as the EMTensor and
//!  matter evolution
/*!
     This class is an example of a matter_t object which calculates the matter
   type specific elements for the RHS update and the evaluation of the
   constraints. This includes the Energy Momentum Tensor, and the matter
   evolution terms. In this case, a vector field, the matter elements are the
   vector field A_\mu, which is decomposed into phi and Avec, and (minus) the
   conjugate momentum of Avec, Evec. It assumes minimal coupling of the field to
   gravity. \sa CCZ4Matter(), ConstraintsMatter()
*/
template <class potential_t> class FixedBGProcaFieldTest
{
  protected:
    double m_vector_mass;    //!< The local copy of the matter param - the mass
    double m_vector_damping; //!< The local copy of the matter param - the
                             //!< damping param
    const potential_t m_potential;

  public:
    //!  Constructor of class FixedBGProcaFieldTest, inputs are the matter
    //!  parameters.
    FixedBGProcaFieldTest(double a_vector_mass, double a_vector_damping,
                          const potential_t potential)
        : m_vector_mass(a_vector_mass), m_vector_damping(a_vector_damping),
          m_potential(potential)
    {
    }

    //! Structure containing all the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        // Vector fields
        data_t phi;
        data_t Z; // Auxilliary variable
        Tensor<1, data_t> Avec;
        Tensor<1, data_t> Evec;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_phi, phi);
            define_enum_mapping(mapping_function, c_Z, Z);
            define_enum_mapping(mapping_function,
                                GRInterval<c_Avec1, c_Avec3>(), Avec);
            define_enum_mapping(mapping_function,
                                GRInterval<c_Evec1, c_Evec3>(), Evec);
        }
    };

    //! Structure containing the rhs variables for the matter fields requiring
    //! 2nd derivs
    template <class data_t> struct Diff2Vars
    {
        Tensor<1, data_t> Avec;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(
                mapping_function, GRInterval<c_Avec1, c_Avec3>(), Avec);
        }
    };

    // Struct for the non grid ADM vars
    template <class data_t>
    using MetricVars = typename ADMFixedBGVars::template Vars<data_t>;

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives
    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars, //!< the value of the variables
        const MetricVars<data_t>
            &metric_vars, //!< the value of the metric variables
        const vars_t<Tensor<1, data_t>> &d1,    //!< the 1st derivs
        const Tensor<2, data_t> &gamma_UU,      //!< the inverse metric (raised)
        const Tensor<3, data_t> &chris_phys_ULL //!< the conformal christoffel
        ) const;

    //! The function which adds in the matter field RHS, given the vars and
    //! derivatives
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void matter_rhs(
        rhs_vars_t<data_t> &total_rhs, //!< RHS terms for all vars.
        const vars_t<data_t> &vars,    //!< the value of the variables
        const MetricVars<data_t>
            &metric_vars, //!< the value of the metric variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< the 2nd derivs
        const vars_t<data_t> &advec) const; //!< the value of beta^i d_i(var).
};

#include "FixedBGProcaFieldTest.impl.hpp"

#endif /* FIXEDBGPROCAFIELDTEST_HPP_ */
