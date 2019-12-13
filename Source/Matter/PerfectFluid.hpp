/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PERFECTFLUID_HPP_
#define PERFECTFLUID_HPP_

#include "CCZ4Geometry.hpp"
#include "DefaultEOS.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"

#include "Cell.hpp" // added for update_fluid_vars()

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



template <class eos_t = DefaultEOS> class PerfectFluid
{
  protected:
    //! The local copy of the potential
    eos_t my_eos;

  public:
    //!  Constructor of class PerfectFluid, inputs are the matter parameters.
    PerfectFluid(const eos_t a_eos) : my_eos(a_eos) {}

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        // fluid variables
        data_t density;
        data_t energy;
        data_t pressure;
        data_t enthalpy;
        data_t u0;
        Tensor<1, data_t> u;

        // defined vars for evol. eq.
        data_t W;
        data_t Z0;
        Tensor<1, data_t> V;

        // evolving vars
        data_t D;
        data_t E;
        Tensor<1, data_t> Z;

        data_t newlapse;                                                             //FIXME:  if one uses lapse it doesn't work, I donno why

       /// Defines the mapping between members of Vars and Chombo grid
       /// variables (enum in User_Variables)
       template <typename mapping_function_t>
       void enum_mapping(mapping_function_t mapping_function)
       {

           VarsTools::define_enum_mapping(mapping_function, c_density,
              density);
           VarsTools::define_enum_mapping(mapping_function, c_energy,
              energy);
           VarsTools::define_enum_mapping(mapping_function, c_pressure,
              pressure);
           VarsTools::define_enum_mapping(mapping_function, c_enthalpy,
              enthalpy);
           VarsTools::define_enum_mapping(mapping_function, c_u0, u0);
           VarsTools::define_enum_mapping(mapping_function,
             GRInterval<c_u1, c_u3>(), u);

           VarsTools::define_enum_mapping(mapping_function, c_W, W);
           VarsTools::define_enum_mapping(mapping_function, c_Z0, Z0);
           VarsTools::define_enum_mapping(mapping_function,
             GRInterval<c_V1, c_V3>(), V);

           VarsTools::define_enum_mapping(mapping_function, c_D, D);
           VarsTools::define_enum_mapping(mapping_function, c_E, E);
           VarsTools::define_enum_mapping(mapping_function,
             GRInterval<c_Z1, c_Z3>(), Z);


           VarsTools::define_enum_mapping(mapping_function, c_lapse,
                newlapse);                                                          //FIXME:  if one uses lapse it doesn't work, I donno why

       }
   };

   //! Structure containing the rhs variables for the matter fields
   template <class data_t> struct GeoVars
   {
       // geometric variables
       data_t lapse;
       data_t chi;
       Tensor<2, data_t> h;
       Tensor<1, data_t> shift;


      /// Defines the mapping between members of Vars and Chombo grid
      /// variables (enum in User_Variables)
      template <typename mapping_function_t>
      void enum_mapping(mapping_function_t mapping_function)
      {
          VarsTools::define_enum_mapping(mapping_function, c_lapse,
            lapse);

            VarsTools::define_enum_mapping(mapping_function, c_chi,
              chi);

          VarsTools::define_enum_mapping(mapping_function,
            GRInterval<c_shift1, c_shift3>(), shift);

          VarsTools::define_symmetric_enum_mapping(mapping_function,
            GRInterval<c_h11, c_h33>(), h);

      }
  };



   //! Structure containing the rhs variables for the matter fields requiring
   // //!  2nd derivs
   template <class data_t> struct Diff2Vars
   {
       // /*  Commented out as no variables needed so far

       data_t D;

       /// Defines the mapping between members of Vars and Chombo grid
       ///  variables (enum in User_Variables)
       template <typename mapping_function_t>
       void enum_mapping(mapping_function_t mapping_function)
       {
           VarsTools::define_enum_mapping(mapping_function, c_D, D);
       }

       // */
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

   //! The compute member which update the non-evolving fluid vars at each point
   //! in the box
   template <class data_t>
   void compute(Cell<data_t> current_cell) const;

};


#include "PerfectFluid.impl.hpp"

#endif /* PERFECTFLUID_HPP_ */
