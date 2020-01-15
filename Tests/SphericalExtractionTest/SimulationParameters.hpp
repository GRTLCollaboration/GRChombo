/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"

// Problem specific includes
#include "SphericalExtraction.hpp"

class SimulationParameters : public ChomboParameters
{
  public:
    SimulationParameters(GRParmParse &pp) : ChomboParameters(pp)
    {
        read_params(pp);
    }

  private:
    void read_params(GRParmParse &pp)
    {
        // extraction params
        dx.fill(coarsest_dx);
        origin.fill(coarsest_dx / 2.0);

        // Extraction params
        pp.load("num_extraction_radii", extraction_params.num_extraction_radii,
                1);

        if (pp.contains("extraction_radii"))
        {
            pp.load("extraction_radii", extraction_params.extraction_radii,
                    extraction_params.num_extraction_radii);
        }
        else
        {
            pp.load("extraction_radius", extraction_params.extraction_radii, 1,
                    L / 4);
        }
        pp.load("num_points_phi", extraction_params.num_points_phi, 8);
        pp.load("num_points_theta", extraction_params.num_points_theta, 16);
        pp.load("extraction_center", extraction_params.center, center);
        pp.load("write_extraction", extraction_params.write_extraction, false);
    }

  public:
    std::array<double, CH_SPACEDIM> origin, dx;
    SphericalExtraction::params_t extraction_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
