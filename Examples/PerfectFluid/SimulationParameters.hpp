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
#include "EquationOfState.hpp"

// #include "ScalarGauss.hpp" FIXME: needed?

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // for regridding
        pp.load("regrid_threshold_K", regrid_threshold_K);
        // pp.load("regrid_threshold_phi", regrid_threshold_phi);  TODO

        // Initial and fluid data
        // TODO initial_params.centerSF = center; // read in SimulationParametersBase

        // Cosmological parameters
        pp.load("G_Newton", G_Newton, 1.0);

        // Predefine ideal Fluids as eos = omega
        pp.load("omega", eos_params.omega, 0.0);


        // FIXME: I think this is useless.  (commented out)
        // // Bool choise on type of parameters: dynamic in time (1) or static (0)
        // pp.load("dyn_rho", eos_params.evo_rho, 1);
        // pp.load("dyn_epsilon", eos_params.evo_epsilon, 0);
        // pp.load("dyn_pressure", eos_params.evo_pressure, 0);
        // pp.load("dyn_u0", eos_params.scalar_mass, 1);
        // pp.load("dyn_u1", eos_params.scalar_mass, 1);
        // pp.load("dyn_u2", eos_params.scalar_mass, 1);
        // pp.load("dyn_u3", eos_params.scalar_mass, 1);


        // Relaxation params  // TODO: implement
        pp.load("relaxtime", relaxtime);
        pp.load("relaxspeed", relaxspeed);
    }

    // Regrid parameters
    Real regrid_threshold_K;
    // Real regrid_threshold_phi;  TODO

    // Initial data for matter and potential
    double G_Newton;
    EquationOfState::params_t eos_params;

    // Relaxation params
    Real relaxtime, relaxspeed;
    
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
