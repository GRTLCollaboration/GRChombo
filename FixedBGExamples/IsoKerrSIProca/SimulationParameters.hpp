/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"
// Problem specific includes:
#include "IsotropicKerrFixedBG.hpp"
#include "Potential.hpp"
#include "SpheroidalExtraction.hpp"

class SimulationParameters : public ChomboParameters
{
  public:
    SimulationParameters(GRParmParse &pp) : ChomboParameters(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // for regridding
        pp.load("nan_check", nan_check, 1);
        pp.load("sigma", sigma, 0.1);
        pp.load("integral_filename", integral_filename);

        // Initial and Kerr data
        pp.load("bh_mass", bg_params.mass);
        pp.load("bh_spin", bg_params.spin);
        pp.load("bh_center", bg_params.center, center);
        pp.load("proca_mass", potential_params.mass);
        pp.load("proca_self_interaction", potential_params.self_interaction);
        pp.load("field_amplitude", field_amplitude);
        pp.load("proca_damping", proca_damping);
        pp.load("r_min", r_min);
        pp.load("r_max", r_max);

        // Extraction params
        extraction_params.num_extraction_radii = 2;
        extraction_params.extraction_radii = {r_min, r_max};
        pp.load("num_points_phi", extraction_params.num_points_phi, 2);
        pp.load("num_points_t", extraction_params.num_points_t, 5);
        if (extraction_params.num_points_t % 2 == 0)
        {
            extraction_params.num_points_t += 1;
            pout() << "Parameter: num_points_t incompatible with Simpson's "
                   << "rule so increased by 1.\n";
        }
        pp.load("extraction_center", extraction_params.center, center);
        pp.load("zaxis_over_xaxis", extraction_params.zaxis_over_xaxis, 1.0);
        pp.load("write_extraction", extraction_params.write_extraction, false);
    }

    // Problem specific parameters
    double field_amplitude;
    double sigma, proca_damping, r_min, r_max;
    int nan_check;
    std::string integral_filename;
    // Collection of parameters necessary for the sims
    IsotropicKerrFixedBG::params_t bg_params;
    Potential::params_t potential_params;
    SpheroidalExtraction::params_t extraction_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
