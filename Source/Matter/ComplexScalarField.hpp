/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXSCALARFIELD_HPP_
#define COMPLEXSCALARFIELD_HPP_

#include "CCZ4Geometry.hpp"
#include "DefaultComplexPotential.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS, total num of components
#include "VarsTools.hpp"

//!  Calculates the matter type specific elements such as
//!  the EMTensor and matter evolution for a complex scalar field
/*!
     This class is an example of a matter_t object which calculates the matter
     type specific elements for the RHS update and the evaluation of the
     constraints. This includes the Energy Momentum Tensor, and the matter
     evolution terms. In this case, a complex scalar field, the matter elements
     are phi and (minus) its conjugate momentum, Pi (re and im).
     It is templated over a potential function potential_t which the user must
     specify in a class, although a default is provided which sets dVdphi and
     V_of_phi to zero. It assumes minimal coupling of the field to gravity.
     \sa CCZ4Matter(), ScalarField(), ConstraintsMatter()
*/
template <class potential_t = DefaultComplexPotential> class ComplexScalarField
{
  protected:
    //! The local copy of the potential
    potential_t my_potential;

  public:
    //!  Constructor of class ComplexScalarField
    ComplexScalarField(const potential_t a_potential)
        : my_potential(a_potential)
    {
    }

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        data_t phi_Re;
        data_t Pi_Re;
        data_t phi_Im;
        data_t Pi_Im;

        /// Defines the mapping between members of Vars and Chombo
        /// grid variables (enum in UserVariables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi_Re, phi_Re);
            VarsTools::define_enum_mapping(mapping_function, c_Pi_Re, Pi_Re);
            VarsTools::define_enum_mapping(mapping_function, c_phi_Im, phi_Im);
            VarsTools::define_enum_mapping(mapping_function, c_Pi_Im, Pi_Im);
        }
    };

    //! Structure containing the rhs variables for fields requiring 2nd derivs
    template <class data_t> struct Diff2Vars
    {
        data_t phi_Re;
        data_t phi_Im;

        /// Defines the mapping between members of Vars and
        /// Chombo grid variables (enum in UserVariables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi_Re, phi_Re);
            VarsTools::define_enum_mapping(mapping_function, c_phi_Im, phi_Im);
        }
    };

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives
    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t>
            &vars, //!< the value of the variables at the point.
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const Tensor<2, data_t> &h_UU, //!< the inverse metric (raised indices)
        const Tensor<3, data_t> &chris_ULL //!< the conformal christoffel symbol
        ) const;

    //! The function which adds in the RHS for the matter field vars
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    void add_matter_rhs(
        vars_t<data_t> &total_rhs,  //!< contains the value of the RHS vars.
        const vars_t<data_t> &vars, //!< the value of the vars at the point.
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs.
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< the 2nd derivs of vars
        const vars_t<data_t> &advec                //!< the advection terms
        ) const;
};

#include "ComplexScalarField.impl.hpp"

#endif /* COMPLEXSCALARFIELD_HPP_ */
