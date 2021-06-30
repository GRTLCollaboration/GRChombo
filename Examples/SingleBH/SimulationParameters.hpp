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

// define to use IsotropicBoostedBH initial data instead of SingleBH ID
#define USE_ISOTROPIC_BOOSTED_BH

#ifdef USE_ISOTROPIC_BOOSTED_BH
#ifdef USE_AHFINDER
#include "AHInitialGuess.hpp"
#endif
#endif

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
        FOR(idir) { bh_params.center[idir] = center[idir] + offset[idir]; }

#ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess, 0.5 * bh_params.mass);
#ifdef USE_ISOTROPIC_BOOSTED_BH
        double r_x = AH_initial_guess, r_y = AH_initial_guess,
               r_z = AH_initial_guess;

        bool ah_use_ellipsoid;
        pp.load("AH_use_ellipsoid", ah_use_ellipsoid, true);
        if (ah_use_ellipsoid)
        {
            double vel = bh_params.momentum[0];
            double contraction = sqrt(1. - vel * vel);
            r_x *= contraction;
        }

        AH_initial_guess_ellipsoid.set_params(r_x, r_y, r_z);
#endif
#endif
    }

    // Initial data
    int activate_extraction;

    // Collection of parameters necessary for initial conditions
    BoostedBH::params_t bh_params;

#ifdef USE_AHFINDER
    double AH_initial_guess;
#ifdef USE_ISOTROPIC_BOOSTED_BH
    AHInitialGuessEllipsoid AH_initial_guess_ellipsoid;
#endif
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
