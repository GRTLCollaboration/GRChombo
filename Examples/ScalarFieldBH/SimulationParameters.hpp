/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "InitialScalarData.hpp"
#include "Potential.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        // Initial scalar field data
        initial_params.center =
            center; // already read in SimulationParametersBase
        pp.load("G_Newton", G_Newton, 1.0);
        pp.load("scalar_amplitude", initial_params.amplitude, 0.1);
        pp.load("bh_mass", initial_params.bh_mass, 1.0);
        pp.load("scalar_mass", potential_params.scalar_mass, 0.1);

#ifdef USE_AHFINDER
        double AH_guess =
            8. * initial_params.amplitude * initial_params.amplitude;
        pp.load("AH_initial_guess", AH_initial_guess, AH_guess);
#endif
    }

    void check_params()
    {
        warn_parameter("scalar_mass", potential_params.scalar_mass,
                       potential_params.scalar_mass <
                           0.2 / coarsest_dx / dt_multiplier,
                       "oscillations of scalar field do not appear to be "
                       "resolved on coarsest level");
        warn_parameter("bh_mass", initial_params.bh_mass,
                       initial_params.bh_mass < 0.1 * L,
                       "gives a BH size greater than 0.1 times the domain L");
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    InitialScalarData::params_t initial_params;
    Potential::params_t potential_params;

#ifdef USE_AHFINDER
    double AH_initial_guess;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
