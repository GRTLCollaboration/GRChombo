/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGPROCACONSTRAINT_HPP_
#define FIXEDBGPROCACONSTRAINT_HPP_

#include "ADMFixedBGVars.hpp"
#include "CCZ4Geometry.hpp"
#include "FixedBGProcaField.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Potential.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS
#include "VarsTools.hpp"
#include "simd.hpp"

//!  Calculates the constraint violation of the proca field
/*!
    This class takes as input the current state of the hypersurface and
    calculates the constraint violation from it.

*/
template <class potential_t, class background_t> class FixedBGProcaConstraint
{
  public:
    // Use the variable definition in the proca matter class
    template <class data_t>
    using MatterVars =
        typename FixedBGProcaField<potential_t>::template Vars<data_t>;

    // Now the non grid ADM vars
    template <class data_t>
    using MetricVars = typename ADMFixedBGVars::template Vars<data_t>;

    //!  Constructor of class FixedBGProcaConstraint, inputs are the matter
    //!  parameters.
    FixedBGProcaConstraint(background_t a_background, double dx,
                           double a_vector_mass, double a_vector_damping,
                           const potential_t potential);

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    double m_vector_mass;    //!< The local copy of the matter param - the mass
    double m_vector_damping; //!< The local copy of the matter param - the
                             //!< damping param
    const FourthOrderDerivatives m_deriv; //!< The derivates of the fields
    const background_t m_background;      //!< The metric background

    const potential_t m_potential; // !< The potential of the Proca Field

    //!  Function that calculates the constraint equations using
    //!  the proca field variables, 1st derivatives and the values of the
    //!  matter background
    template <class data_t, template <typename> class vars_t>
    data_t constraint_equations(
        const vars_t<data_t> &vars, //!< the value of the variables
        const MetricVars<data_t>
            &metric_vars, //!< the value of the metric variables
        const vars_t<Tensor<1, data_t>> &d1 //!< the value of 1st derivatives
        ) const;
};

#include "FixedBGProcaConstraint.impl.hpp"

#endif /* FIXEDBGPROCACONSTRAINT_HPP_ */
