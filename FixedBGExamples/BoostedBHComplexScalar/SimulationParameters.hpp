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
#include "BoostedBHFixedBG.hpp"
#include "ComplexScalarPotential.hpp"
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
        pp.load("regrid_length", regrid_length, L);
        pp.load("integral_filename", integral_filename);

        // Initial and Kerr data
        pp.load("bh_mass", bg_params.mass);
        pp.load("bh_velocity", bg_params.velocity);
        pp.load("bh_center", bg_params.center, center);
        pp.load("field_amplitude_re", field_amplitude_re);
        pp.load("field_amplitude_im", field_amplitude_im);
        pp.load("scalar_mass", potential_params.scalar_mass);

        // extraction params
        //        dx.fill(coarsest_dx);
        //        origin.fill(coarsest_dx / 2.0);

        // Extraction params
        pp.load("num_extraction_radii", extraction_params.num_extraction_radii,
                1);
        // Check for multiple extraction radii, otherwise load single
        // radius/level (for backwards compatibility).
        if (pp.contains("extraction_levels"))
        {
            pp.load("extraction_levels", extraction_params.extraction_levels,
                    extraction_params.num_extraction_radii);
        }
        else
        {
            pp.load("extraction_level", extraction_params.extraction_levels, 1,
                    0);
        }
        if (pp.contains("extraction_radii"))
        {
            pp.load("extraction_radii", extraction_params.extraction_radii,
                    extraction_params.num_extraction_radii);
        }
        else
        {
            pp.load("extraction_radius", extraction_params.extraction_radii, 1,
                    0.1);
        }
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
    double field_amplitude_re, field_amplitude_im, regrid_length;
    double sigma, proca_mass, proca_damping, excision_width;
    int nan_check;
    //    std::array<double, CH_SPACEDIM> origin,
    //        dx; // location of coarsest origin and dx
    std::string integral_filename;
    // Collection of parameters necessary for the sims
    BoostedBHFixedBG::params_t bg_params;
    SpheroidalExtraction::params_t extraction_params;
    ComplexScalarPotential::params_t potential_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
