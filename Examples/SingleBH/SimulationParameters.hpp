/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "BoostedBH.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        readParams(pp);
    }

    /// Read parameters from the parameter file
    void readParams(GRParmParse &pp)
    {
        // Initial data
        pp.load("mass", bh_params.mass);
        pp.load("momentum", bh_params.momentum);
        pp.load("activate_extraction", activate_extraction, 0);

        // Get the centers of the BHs either explicitly or as
        // an offset (not both, or they will be offset from center
        // provided)
        std::array<double, CH_SPACEDIM> offset;
        pp.load("offset", offset, {0.0, 0.0, 0.0});
        FOR1(idir) { bh_params.center[idir] = center[idir] + offset[idir]; }

#ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess, 0.5 * bh_params.mass);
#endif
    }

    // Initial data
    int activate_extraction;

    // Collection of parameters necessary for initial conditions
    BoostedBH::params_t bh_params;

#ifdef USE_AHFINDER
    double AH_initial_guess;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
