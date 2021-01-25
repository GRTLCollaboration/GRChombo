/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "FixedBGSimulationParametersBase.hpp"
#include "GRParmParse.hpp"
// Problem specific includes:
#include "ComplexScalarPotential.hpp"
#include "NewtonianBHFixedBG.hpp"

class SimulationParameters : public FixedBGSimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : FixedBGSimulationParametersBase(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // for regridding
        pp.load("regrid_length", regrid_length, L);
        pp.load("integral_filename", integral_filename);

        // Initial data
        pp.load("bh_mass_1", bg_params1.mass);
        pp.load("bh_center_1", bg_params1.center, center);
        pp.load("bh_mass_2", bg_params2.mass);
        pp.load("bh_center_2", bg_params2.center, center);
        pp.load("omega_binary", omega_binary);
        pp.load("separation", separation, 0.0);
        bg_params1.center[0] += separation;
        bg_params2.center[0] -= separation;
        pp.load("field_amplitude_re", field_amplitude_re);
        pp.load("field_amplitude_im", field_amplitude_im);
        pp.load("scalar_mass", potential_params.scalar_mass);
        pp.load("inner_r", inner_r, 1.0);
        pp.load("outer_r", outer_r, L / 2.0);
    }

    // Problem specific parameters
    double field_amplitude_re, field_amplitude_im, regrid_length;
    double proca_mass, proca_damping;
    double inner_r, outer_r, omega_binary, separation;
    std::string integral_filename;

    // Collection of parameters necessary for the sims
    NewtonianBHFixedBG::params_t bg_params1;
    NewtonianBHFixedBG::params_t bg_params2;
    ComplexScalarPotential::params_t potential_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
