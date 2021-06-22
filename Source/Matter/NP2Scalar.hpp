/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef NP2SCALAR_HPP_
#define NP2SCALAR_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "ProcaField.hpp"
#include "SphericalHarmonics.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS
#include "VarsTools.hpp"
#include "simd.hpp"

//! Struct for the null tetrad
template <class data_t> struct Tetrad_t
{
    Tensor<1, data_t> u; //!< the vector u^i
    Tensor<1, data_t> v; //!< the vector v^i
    Tensor<1, data_t> w; //!< the vector w^i
};

//! Struct for the Newman Penrose scalar
template <class data_t> struct NPScalar_t
{
    data_t Real;      // Real component
    data_t Im;        // Imaginary component
    data_t magnitude; // Abs value, ie sqrt(Re*Re + Im*Im)
};

//!  Calculates the Newman Penrose scalars for a Proca field.
/*!
     The class calculates the Newman Penrose Scalars for a Proca field. At the
   moment it only does Phi_2 but we can add to this in future.
*/
class NP2Scalar
{
  public:
    //! Structure containing the necessary variables for the fields
    template <class data_t> struct Vars
    {
        data_t chi;
        Tensor<1, data_t> Evec;
        Tensor<2, data_t> h;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools;
            define_enum_mapping(mapping_function, c_chi, chi);
            define_enum_mapping(mapping_function,
                                GRInterval<c_Evec1, c_Evec3>(), Evec);
            define_symmetric_enum_mapping(mapping_function,
                                          GRInterval<c_h11, c_h33>(), h);
        }
    };

    //! Structure containing the rhs variables for the matter fields requiring
    //! derivs
    template <class data_t> struct Diff1Vars
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

    //! Constructor of class NP2Scalar
    NP2Scalar(double dx, double L) : m_dx(dx), m_L(L), m_deriv(dx) {}

    //! The compute member which calculates the wave quantities at each point on
    //! the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const double m_dx; //!< the grid spacing
    const double m_L;  //!< The length of the grid, used for midpoint

    //! Calculation of the Newman Penrose scalar Psi_2 for a vector field Avec
    /*!
        Phi_2 = F_{\mu \nu} * l^\mu \bar(m)^\nu, see e.g. arxiv 1212.0551 eqn 32
        NB this assumes the field is purely real
    */
    template <class data_t>
    NPScalar_t<data_t>
    compute_NP_scalar_2(const Vars<data_t> &vars,
                        const Diff1Vars<Tensor<1, data_t>> &d1,
                        const Coordinates<data_t> &coords) const;

    //! Calculation of the null tetrad
    template <class data_t>
    Tetrad_t<data_t>
    compute_null_tetrad(const Vars<data_t> &vars,
                        const Coordinates<data_t> &coords) const;
};

#include "NP2Scalar.impl.hpp"

#endif /* NP2SCALAR_HPP_ */
