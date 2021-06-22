/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"

// Problem specific includes:
#include "SimulationParametersBase.hpp"
#include "WarpBubble.hpp"
#include "WarpField.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // Initial and SF data
        pp.load("G_Newton", G_Newton, 1.0);
        pp.load("warp_speed", warp_params.warp_speed, 0.5);
        pp.load("warp_acceleration", warp_params.acceleration, 0.0);
        pp.load("bubble_center", warp_params.bubble_center, center);
        pp.load("bubble_size", warp_params.bubble_size, 10.0);
        pp.load("sigma_wall", warp_params.sigma_wall, 1.0);
        pp.load("a1", warpfield_params.a1, 0.0);
        pp.load("a2", warpfield_params.a2, 0.0);
        pp.load("a3", warpfield_params.a3, 0.0);
    }

    // Collection of parameters necessary for the Warp drive
    double G_Newton;
    WarpBubble::params_t warp_params;
    WarpField::params_t warpfield_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
