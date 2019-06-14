/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PERFECTFLUID_HPP_
#define PERFECTFLUID_HPP_

#include "CCZ4Geometry.hpp"
#include "DefaultEquationOfState.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS, total num of components
#include "VarsTools.hpp"

//!  Calculates the matter type specific elements such as the EMTensor and
//   matter evolution
/*!
     This class is an example of a matter_t object which calculates the
     matter type specific elements for the RHS update and the evaluation
     of the constraints. This includes the Energy Momentum Tensor, and
     the matter evolution terms. In this case, a Perfect Fluid,
     the primary matter elements are: rho (rest mass density), epsilon
     (internal energy) and u-ith (the 4-velocity of the fluid), both pressure and
     enthalpy depend on rho and epsilon. For convinence, other derived
     intermidiete derived variable W, D, Vi are defined for the time evolution.
     This class is templated over a EOS function EOS_t which the
     user must specify in a class, although a default is provided which
     sets trivial pressure and enthalpy as for pressurless matter (dust).
     It assumes minimal coupling of the field to gravity.
     \sa MatterCCZ4(), ConstraintsMatter()
*/


/* TODO: FIXME:  Here I annotate replacements
scalar_field --> perfect_field
ScalarField --> PerfectFluid
potential_t --> eos_t
DefaultPotential --> DefaultEquationOfState
my_potential --> my_eos
SFObject --> FluidObject

*/



template <class eos_t = DefaultEquationOfState> class PerfectFluid
{
  protected:
    //! The local copy of the potential
    eos_t my_eos;

  public:
    //!  Constructor of class PerfectFluid, inputs are the matter parameters.
    PerfectFluid(const eos_t a_eos) : my_eos(a_eos) {}

    //! Structure containing the variables for the matter fields
    template <class data_t> struct FluidObject                                    //FIXME: needed?
    {
        data_t density;
        data_t energy;
        data_t pressure;
        data_t enthalpy;
        data_t u0;
        Tensor<1, data_t> u;

        data_t W;
        data_t Z0;
        Tensor<1, data_t> V;

        data_t D;
        data_t E;
        Tensor<1, data_t> Z;
    };

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        data_t density;
        data_t energy;
        data_t pressure;
        data_t enthalpy;
        data_t u0;
        Tensor<1, data_t> u;

        data_t W;
        data_t Z0;
        Tensor<1, data_t> V;

        data_t D;
        data_t E;
        Tensor<1, data_t> Z;                                                              //FIXME: Stopped  coding here!

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi, phi);
            VarsTools::define_enum_mapping(mapping_function, c_Pi, Pi);
        }
    };

    //! Structure containing the rhs variables for the matter fields requiring
    //!  2nd derivs
    template <class data_t> struct Diff2Vars
    {
        data_t phi;

        /// Defines the mapping between members of Vars and Chombo grid
        ///  variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi, phi);
        }
    };

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, including the potential
    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars,          //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const Tensor<2, data_t> &h_UU, //!< the inverse metric (raised indices)
        const Tensor<3, data_t> &chris_ULL) //!< the conformal chris. symbol
        const;

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, excluding the potential
    template <class data_t, template <typename> class vars_t>
    static void emtensor_excl_potential(
        emtensor_t<data_t> &out,         //!< the em tensor output
        const vars_t<data_t> &vars,      //!< the value of the variables
        const FluidObject<data_t> &vars_sf, //!< the value of the sf variables
        const Tensor<1, data_t>
            &d1_phi,                   //!< the value of the first deriv of phi
        const Tensor<2, data_t> &h_UU, //!< the inverse metric (raised indices).
        const Tensor<3, data_t>
            &chris_ULL); //!< the conformal christoffel symbol

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
        const vars_t<data_t> &advec)
        const; //!< the value of the advection terms

    //! The function which calculates the RHS for the matter field vars
    //! excluding the potential
    template <class data_t, template <typename> class vars_t>
    static void matter_rhs_excl_potential(
        FluidObject<data_t>
            &rhs_sf, //!< the value of the RHS terms for the sf vars
        const vars_t<data_t> &vars,      //!< the values of all the variables
        const FluidObject<data_t> &vars_sf, //!< the value of the sf variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const Tensor<1, data_t> &d1_phi, //!< the value of the 1st derivs of phi
        const Tensor<2, data_t> &d2_phi, //!< the value of the 2nd derivs of phi
        const FluidObject<data_t> &advec_sf); //!< advection terms for the sf vars
};

#include "PerfectFluid.impl.hpp"

#endif /* PERFECTFLUID_HPP_ */
