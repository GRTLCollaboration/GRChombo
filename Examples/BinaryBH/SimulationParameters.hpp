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
#include "BoostedBH.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        readParams(pp);
    }

    /// Read parameters from the parameter file
    void readParams(GRParmParse &pp)
    {
        // The automatically generated read parameters code defined in
        // SimulationParameters.inc
        auto_read_params(pp);

        // Fill in BinaryBHParameters
        bh1_params.mass = massA;
        bh1_params.center = centerA;
        bh1_params.momentum = momentumA;
        bh2_params.mass = massB;
        bh2_params.center = centerB;
        bh2_params.momentum = momentumB;
    }

    void auto_read_params(GRParmParse &pp)
    {
        // Initial data
        pp.load("massA", massA);
        pp.load("centerA", centerA);
        pp.load("momentumA", momentumA);
        pp.load("massB", massB);
        pp.load("centerB", centerB);
        pp.load("momentumB", momentumB);
        pp.load("activate_extraction", activate_extraction, 0);

        // Fill in BinaryBHParameters
        bh1_params.mass = massA;
        bh1_params.center = centerA;
        bh1_params.momentum = momentumA;
        bh2_params.mass = massB;
        bh2_params.center = centerB;
        bh2_params.momentum = momentumB;
    }

    // Initial data
    int activate_extraction;
    Real massA, massB;
    std::array<double, CH_SPACEDIM> centerA, centerB;
    std::array<double, CH_SPACEDIM> momentumA, momentumB;

    // Collection of parameters necessary for initial conditions
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
