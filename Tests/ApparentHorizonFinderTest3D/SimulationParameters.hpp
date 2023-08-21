/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"

#include "KerrBH.hpp"

#ifdef USE_AHFINDER
#include "AHFinder.hpp"
#endif

class SimulationParameters : public ChomboParameters
{
  public:
    SimulationParameters(GRParmParse &pp) : ChomboParameters(pp)
    {
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // Initial Kerr data
        pp.load("kerr_mass", kerr_params.mass);
        pp.load("kerr_spin", kerr_params.spin);
        pp.load("kerr_center", kerr_params.center, center);
        pp.load("kerr_spin_direction", kerr_params.spin_direction,
                {0., 0., 1.});

#ifdef USE_AHFINDER
        pp.load("AH_initial_guess", initial_guess, kerr_params.mass * 0.5);
        AH_params.read_params(pp, *this);
#endif
    }
    KerrBH::params_t kerr_params;

#ifdef USE_AHFINDER
    double initial_guess;
    AHParams_t<AHFunction> AH_params;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
