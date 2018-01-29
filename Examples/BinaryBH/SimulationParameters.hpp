/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"

// Problem specific includes:
#include "BoostedBH.hpp"
#include "CCZ4.hpp"

class SimulationParameters
{
  public:
    SimulationParameters(GRParmParse &pp) { readParams(pp); }

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

        // Fill in he ccz4Parameters
        ccz4_params.kappa1 = kappa1;
        ccz4_params.kappa2 = kappa2;
        ccz4_params.kappa3 = kappa3;
        ccz4_params.shift_Gamma_coeff = shift_Gamma_coeff;
        ccz4_params.shift_advec_coeff = shift_advec_coeff;
        ccz4_params.eta = eta;
        ccz4_params.lapse_advec_coeff = lapse_advec_coeff;
        ccz4_params.lapse_power = lapse_power;
        ccz4_params.lapse_coeff = lapse_coeff;
    }

// SimulationParameters.inc declares all variables and defines
// auto_read_params(GRParmParse& pp)
#include "SimulationParameters.inc"

    // Collection of parameters necessary for the CCZ4 RHS
    CCZ4::params_t ccz4_params;
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
