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
    }

    void read_params(GRParmParse &pp)
    {
        // Initial scalar field data
        initial_params.center =
            center; // already read in SimulationParametersBase
        pp.load("G_Newton", G_Newton, 1.0);
        pp.load("scalar_amplitude", initial_params.amplitude);
        pp.load("scalar_width", initial_params.width);
        pp.load("scalar_mass", potential_params.scalar_mass);

        // Initial Kerr data
        pp.load("kerr_mass", kerr_params.mass, 1.0);
        pp.load("kerr_spin", kerr_params.spin, 0.0);
        pp.load("kerr_center", kerr_params.center, center);
    }

    void check_params()
    {
        warn_parameter("kerr_mass", kerr_params.mass, kerr_params.mass >= 0.0,
                       "should be >= 0.0");
        check_parameter("kerr_spin", kerr_params.spin,
                        std::abs(kerr_params.spin) <= kerr_params.mass,
                        "must satisfy |a| <= M = " +
                            std::to_string(kerr_params.mass));
        FOR1(idir)
        {
            std::string name = "kerr_center[" + std::to_string(idir) + "]";
            warn_parameter(
                name, kerr_params.center[idir],
                (kerr_params.center[idir] >= 0) &&
                    (kerr_params.center[idir] <= (ivN[idir] + 1) * coarsest_dx),
                "should be within the computational domain");
        }
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    InitialScalarData::params_t initial_params;
    Potential::params_t potential_params;
    KerrBH::params_t kerr_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
