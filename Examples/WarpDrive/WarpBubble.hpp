/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WARPBUBBLE_HPP_
#define WARPBUBBLE_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "WarpField.hpp"
#include "simd.hpp"

//! Class which creates a bubble of a Warp field given params for initial
//! matter config
class WarpBubble
{
  public:
    //! A structure for the input params for Warp field properties and initial
    //! conditions
    struct params_t
    {
        double warp_speed;         //!< Amplitude of bump in initial warp bubble
        double acceleration = 0.0; //!< Amplitude of bump in initial warp bubble
        std::array<double, CH_SPACEDIM>
            bubble_center;  //!< Centre of perturbation in initial warp bubble
        double bubble_size; //!< Width of initial warp bubble
        double sigma_wall;  //!< Width of wall in initial warp bubble
    };

    //! The constructor
    WarpBubble(params_t a_params, double a_dx) : m_dx(a_dx), m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const double m_dx;
    const params_t m_params; //!< The matter initial condition params
};

#include "WarpBubble.impl.hpp"

#endif /* WARPBUBBLE_HPP_ */
