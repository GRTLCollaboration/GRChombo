/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SCALARSIN3D_HPP_
#define SCALARSIN3D_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates a bubble of a scalar field given params for initial
//! matter config
class ScalarSin3d
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double amplitudeSF; //!< Amplitude of bump in initial SF bubble
        std::array<double, CH_SPACEDIM>
            centerSF;   //!< Centre of perturbation in initial SF bubble
        double widthSF; //!< Width of bump in initial SF bubble
        double r_zero;  //!< Position of bump relative to centre
  	    double bkgSF;  // CJ:  mean homogeneous SF (background)
        double kfactor = 0.0;  // CJ:  mean homogeneous SF (background)
    };

    //! The constructor
    ScalarSin3d(params_t a_params, double a_dx);

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params

    //! Function to compute the value of phi at each point
    template <class data_t>
    data_t compute_phi(Coordinates<data_t> coords) const;

    //! Function to compute the value of K at each point (constant for now)
    template <class data_t>
    data_t compute_K(Coordinates<data_t> coords) const;  // CJ New
};

#include "ScalarSin3d.impl.hpp"

#endif /* SCALARGSIN3D_HPP_ */
