/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "FixedBGSimulationParametersBase.hpp"
#include "GRParmParse.hpp"
// Problem specific includes:
#include "BoostedBHFixedBG.hpp"
#include "ComplexPotential.hpp"

class SimulationParameters : public FixedBGSimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : FixedBGSimulationParametersBase(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // for regridding
        pp.load("regrid_length", regrid_length, L);

        // BH data
        pp.load("bh_mass", bg_params.mass);
        pp.load("bh_velocity", bg_params.velocity);
        pp.load("bh_center", bg_params.center, center);

        // Initial SF
        pp.load("scalar_amplitude", scalar_amplitude);
        pp.load("scalar_mass", scalar_mass);

        // Volume extraction radii
        pp.load("inner_r", inner_r, 5.0);
        pp.load("outer_r", outer_r, 100.0 / scalar_mass);
    }

    // Problem specific parameters
    double scalar_amplitude, scalar_mass, regrid_length;
    double inner_r, outer_r;
    // Collection of parameters necessary for the sims
    BoostedBHFixedBG::params_t bg_params;
    SphericalExtraction::params_t extraction_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
