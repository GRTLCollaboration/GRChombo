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
        // for regridding
        pp.load("regrid_threshold_chi", regrid_threshold_chi);
        pp.load("regrid_threshold_phi", regrid_threshold_phi);

        // Initial and SF data
        initial_params.centerSF =
            center; // already read in SimulationParametersBase
        pp.load("G_Newton", G_Newton, 1.0);
        pp.load("amplitudeSF", initial_params.amplitudeSF);
        pp.load("widthSF", initial_params.widthSF);
        pp.load("r_zero", initial_params.r_zero);
        pp.load("scalar_mass", potential_params.scalar_mass);

        // Relaxation params
        pp.load("relaxtime", relaxtime);
        pp.load("relaxspeed", relaxspeed);

        // whether to do extraction
        pp.load("activate_extraction", activate_extraction, false);
    }

    // Problem specific parameters
    Real regrid_threshold_chi, regrid_threshold_phi;

    // Initial data for matter and potential
    double G_Newton;
    ScalarBubble::params_t initial_params;
    Potential::params_t potential_params;

    // Relaxation params
    Real relaxtime, relaxspeed;

    // do extraction?
    bool activate_extraction;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
