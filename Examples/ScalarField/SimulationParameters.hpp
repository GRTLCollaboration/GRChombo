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
#include "ScalarBubble.hpp"

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
        // Reads the problem specific params
        auto_read_params(pp);

        // Fill in the Matter Parameters
        initial_params.amplitudeSF = amplitudeSF;
        initial_params.centerSF = centerSF;
        initial_params.widthSF = widthSF;
        initial_params.r_zero = r_zero;

        // Fill in the potential parameters
        potential_params.scalar_mass = scalar_mass;
    }

    void auto_read_params(GRParmParse &pp)
    {
        pp.load("regrid_threshold_K", regrid_threshold_K);
        pp.load("regrid_threshold_phi", regrid_threshold_phi);

        // Initial and SF data
        pp.load("G_Newton", G_Newton, 1.0);
        pp.load("amplitudeSF", amplitudeSF);
        pp.load("widthSF", widthSF);
        pp.load("r_zero", r_zero);
        pp.load("scalar_mass", scalar_mass);

        // Relaxation params
        pp.load("relaxtime", relaxtime);
        pp.load("relaxspeed", relaxspeed);
    }

    // Problem specific parameters
    Real regrid_threshold_K, regrid_threshold_phi;
    // Initial data for matter and potential
    double G_Newton;
    Real amplitudeSF, widthSF, r_zero, scalar_mass;
    std::array<double, CH_SPACEDIM> centerSF;
    // Relaxation params
    Real relaxtime, relaxspeed;

    // Collection of parameters necessary for the problem
    ScalarBubble::params_t initial_params;
    Potential::params_t potential_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
