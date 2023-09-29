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
#include "MeansVars.hpp"

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
        grid_params.center = center;
        pp.load("G_Newton", G_Newton,
                0.0); // for now the example neglects backreaction
        //pp.load("max_level", grid_params.max_level, 0);
        pp.load("scalar_amplitude", initial_params.amplitude, 10.0);
        pp.load("scalar_velocity", initial_params.velocity, -0.001);
        pp.load("scalar_width", initial_params.width, 1.0);
        pp.load("scalar_mass", potential_params.scalar_mass, 0.1);
        pp.load("scalar_mass", initial_params.mass, 1e-6);
        pp.load("N_full", initial_params.N_init, 0);

        // Initial Kerr data
        /*pp.load("kerr_mass", kerr_params.mass, 1.0);
        pp.load("kerr_spin", kerr_params.spin, 0.0);
        pp.load("kerr_center", kerr_params.center, center);*/
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
    double G_Newton;
    double m_pl;
    //CalcMeans::params_t grid_params;
    InitialScalarData::params_t initial_params;
    Potential::params_t potential_params;
    //KerrBH::params_t kerr_params;
    MeansVars::params_t grid_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
