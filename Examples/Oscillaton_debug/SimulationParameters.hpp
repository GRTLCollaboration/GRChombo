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
#include "CCZ4.hpp"
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
        pp.load("G_Newton", G_Newton, 1.0);
        pp.load("scalar_mass", potential_params.scalar_mass, 1.0);
        pp.load("initial_data_prefix", initial_data_prefix);
    }

    void check_params()
    {
        warn_parameter("scalar_mass", potential_params.scalar_mass,
                       potential_params.scalar_mass < 0.2 / coarsest_dx /
                                                          dt_multiplier *
                                                          pow(2.0, max_level),
                       "oscillations of scalar field do not appear to be "
                       "resolved on finest level");
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    Potential::params_t potential_params;
    std::string initial_data_prefix;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
