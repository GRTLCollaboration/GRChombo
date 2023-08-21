/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"

#ifdef USE_AHFINDER
#include "AHFinder.hpp"
#endif

class SimulationParameters : public ChomboParameters
{
  public:
    SimulationParameters(GRParmParse &pp) : ChomboParameters(pp)
    {
#ifdef USE_AHFINDER
        AH_params.read_params(pp, *this);
#endif
    }

#ifdef USE_AHFINDER
    AHParams_t<AHFunction> AH_params;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
