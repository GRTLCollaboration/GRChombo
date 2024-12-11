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
#include "KerrBH.hpp"
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
        pp.load("scalar_mass", potential_params.scalar_mass, 0.1);

        // Lineout params
        pp.load("lineout_num_points", lineout_num_points, 10);

        // Tagging params
        pp.load("tagging_center", tagging_center, center);
        pp.load("tagging_radius", tagging_radius, L);

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
    }

    // Initial data for matter and potential and BH
    double G_Newton, tagging_radius;
    int lineout_num_points;
    std::array<double, CH_SPACEDIM> tagging_center;
    InitialScalarData::params_t initial_params;
    Potential::params_t potential_params;
    KerrBH::params_t kerr_params;

#ifdef USE_AHFINDER
    double AH_initial_guess;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
