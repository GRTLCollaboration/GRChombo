/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXSCALARBUBBLE_HPP_
#define COMPLEXSCALARBUBBLE_HPP_

#include "Cell.hpp"
#include "ComplexScalarField.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates a bubble of a ComplexScalar field - re and im out of
//! phase
class ComplexScalarBubble
{
  public:
    //! A structure for the input params for scalar field initial conditions
    struct params_t
    {
        double amplitudeSF; //!< Amplitude of bump in initial SF bubble
        std::array<double, CH_SPACEDIM> centerSF; //!< Centre of perturbation
        double widthSF; //!< Width of bump in initial SF bubble
        double r_zero;  //!< Position of bump relative to centre
    };

    //! The constructor for the ComplexScalarField class
    ComplexScalarBubble(params_t a_params, double a_dx);

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const double m_dx;
    const params_t m_params; //!< The matter initial condition params

    //! Function to compute the value of phi at each point
    template <class data_t>
    data_t compute_gaussian(Coordinates<data_t> coords) const;
};

#include "ComplexScalarBubble.impl.hpp"

#endif /* COMPLEXSCALARBUBBLE_HPP_ */
