/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "FixedBGSimulationParametersBase.hpp"

// Problem specific includes:
#include "InitialScalarData.hpp"
#include "BoostedBHFixedBG.hpp"

class SimulationParameters : public FixedBGSimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : FixedBGSimulationParametersBase(pp)
    {
        // read the problem specific params
        readParams(pp);
        check_params();
    }

    void readParams(GRParmParse &pp)
    {
        // Initial SF
        pp.load("scalar_amplitude", initial_params.amplitude, 0.1);
        pp.load("scalar_mass", initial_params.mass, 0.1);

	// BH data
        pp.load("bh_mass", bg_params.mass, 1.0);
        pp.load("bh_velocity", bg_params.velocity);
        pp.load("bh_center", bg_params.center, center);

        // Volume extraction radii
        pp.load("inner_r", inner_r, 5.0);
        pp.load("outer_r", outer_r, 100.0 / initial_params.mass);
    }

    void check_params()
    {
        warn_parameter("scalar_mass", initial_params.mass,
                       initial_params.mass <
                           0.2 / coarsest_dx / dt_multiplier,
                       "oscillations of scalar field do not appear to be "
                       "resolved on coarsest level");
        warn_parameter("bh_mass", bg_params.mass, bg_params.mass >= 0.0,
                       "should be >= 0.0");
        FOR(idir)
        {
            std::string name = "bh_center[" + std::to_string(idir) + "]";
            warn_parameter(
                name, bg_params.center[idir],
                (bg_params.center[idir] >= 0) &&
                    (bg_params.center[idir] <= (ivN[idir] + 1) * coarsest_dx),
                "should be within the computational domain");
        }
    }
  
    // Problem specific parameters
    double inner_r, outer_r;
    InitialScalarData::params_t initial_params;
    // Collection of parameters necessary for the sims
    BoostedBHFixedBG::params_t bg_params;
    SphericalExtraction::params_t extraction_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
