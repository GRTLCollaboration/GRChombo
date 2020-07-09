/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WARPFIELD_HPP_
#define WARPFIELD_HPP_

#include "CCZ4Geometry.hpp"
#include "DefaultPotential.hpp"
#include "DimensionDefinitions.hpp"
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
     the matter evolution terms. In this case, a warp drive.
     \sa MatterCCZ4(), ConstraintsMatter()
*/
class WarpField
{
  public:
    struct params_t
    {
        double a1;
        double a2;
        double a3;
    };

  protected:
    params_t m_params;
    double m_time;
    double m_xs;
    double m_vs;

  public:
    //!  Constructor of class WarpField, input is time so vs can update
    WarpField(params_t a_params, double a_time, double a_xs, double a_vs)
        : m_time(a_time), m_xs(a_xs), m_vs(a_vs), m_params(a_params)
    {
    }

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        Tensor<2, data_t> S_tensor;
        Tensor<1, data_t> S_vector;
        data_t rho;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_rho, rho);
            VarsTools::define_enum_mapping(mapping_function,
                                           GRInterval<c_S1, c_S3>(), S_vector);
            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_S11, c_S33>(), S_tensor);
        }
    };

    //! Structure containing the rhs variables for the matter fields requiring
    //!  2nd derivs
    template <class data_t> struct Diff2Vars
    {
        data_t rho;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_rho, rho);
        }
    };

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives
    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars,          //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const Tensor<2, data_t> &h_UU, //!< the inverse metric (raised indices)
        const Tensor<3, data_t> &chris_ULL)
        const; //!< the conformal christoffel symbol

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
};

#include "WarpField.impl.hpp"

#endif /* WARPFIELD_HPP_ */
