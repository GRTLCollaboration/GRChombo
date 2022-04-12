/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGSIMULATIONPARAMETERSBASE_HPP_
#define FIXEDBGSIMULATIONPARAMETERSBASE_HPP_

// General includes
#include "BoundaryConditions.hpp"
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"
#include "SphericalExtraction.hpp"

// add this type alias here for backwards compatibility
using extraction_params_t = SphericalExtraction::params_t;

class FixedBGSimulationParametersBase : public ChomboParameters
{
  public:
    FixedBGSimulationParametersBase(GRParmParse &pp) : ChomboParameters(pp)
    {
        read_params(pp);
        check_params();
    }

  private:
    void read_params(GRParmParse &pp)
    {
        // Dissipation
        pp.load("sigma", sigma, 0.1);

        // Nan Check and min chi and lapse values
        pp.load("nan_check", nan_check, 1);

        // directory to store data
        pp.load("data_subpath", data_path, std::string(""));
        if (!data_path.empty() && data_path.back() != '/')
            data_path += "/";
        if (output_path != "./" && !output_path.empty())
            data_path = output_path + data_path;

        // Extraction params
        pp.load("activate_extraction", activate_extraction, false);

        if (activate_extraction)
        {
            pp.load("num_extraction_radii",
                    extraction_params.num_extraction_radii, 1);
            // Check for multiple extraction radii, otherwise load single
            // radius/level (for backwards compatibility).
            if (pp.contains("extraction_levels"))
            {
                pp.load("extraction_levels",
                        extraction_params.extraction_levels,
                        extraction_params.num_extraction_radii);
            }
            else
            {
                pp.load("extraction_level", extraction_params.extraction_levels,
                        1, 0);
            }
            if (pp.contains("extraction_radii"))
            {
                pp.load("extraction_radii", extraction_params.extraction_radii,
                        extraction_params.num_extraction_radii);
            }
            else
            {
                pp.load("extraction_radius", extraction_params.extraction_radii,
                        1, 0.1);
            }

            pp.load("num_points_phi", extraction_params.num_points_phi, 2);
            pp.load("num_points_theta", extraction_params.num_points_theta, 5);
            if (extraction_params.num_points_theta % 2 == 0)
            {
                extraction_params.num_points_theta += 1;
                pout() << "Parameter: num_points_theta incompatible with "
                          "Simpson's "
                       << "rule so increased by 1.\n";
            }
            pp.load("extraction_center", extraction_params.center, center);

            if (pp.contains("modes"))
            {
                pp.load("num_modes", extraction_params.num_modes);
                std::vector<int> extraction_modes_vect(
                    2 * extraction_params.num_modes);
                pp.load("modes", extraction_modes_vect,
                        2 * extraction_params.num_modes);
                extraction_params.modes.resize(extraction_params.num_modes);
                for (int i = 0; i < extraction_params.num_modes; ++i)
                {
                    extraction_params.modes[i].first =
                        extraction_modes_vect[2 * i];
                    extraction_params.modes[i].second =
                        extraction_modes_vect[2 * i + 1];
                }
            }
            else
            {
                // by default extraction (l,m) = (2,0), (2,1) and (2,2)
                extraction_params.num_modes = 3;
                extraction_params.modes.resize(3);
                for (int i = 0; i < 3; ++i)
                {
                    extraction_params.modes[i].first = 2;
                    extraction_params.modes[i].second = i;
                }
            }

            pp.load("write_extraction", extraction_params.write_extraction,
                    false);
        }
    }

    void check_params()
    {
        check_parameter("dt_multiplier", dt_multiplier, dt_multiplier < 1.0,
                        "must be < 1.0 for stability");
        warn_parameter("dt_multiplier", dt_multiplier, dt_multiplier <= 0.5,
                       "is unlikely to be stable for > 0.5");

        check_parameter("sigma", sigma,
                        (sigma >= 0.0) && (sigma <= 2.0 / dt_multiplier),
                        "must be >= 0.0 and <= 2 / dt_multiplier for stability "
                        "(see Alcubierre p344)");
        warn_parameter("nan_check", nan_check, nan_check,
                       "should not normally be disabled");

        if (activate_extraction)
        {
            check_parameter(
                "num_extraction_radii", extraction_params.num_extraction_radii,
                extraction_params.num_extraction_radii > 0,
                "must be bigger than 0 when activate_extraction = 1");

            FOR(idir)
            {
                std::string center_name =
                    "extraction_center[" + std::to_string(idir) + "]";
                double center_in_dir = extraction_params.center[idir];
                check_parameter(
                    center_name, center_in_dir,
                    (center_in_dir >= reflective_domain_lo[idir]) &&
                        (center_in_dir <= reflective_domain_hi[idir]),
                    "must be in the computational domain after "
                    "applying reflective symmetry");
                for (int iradius = 0;
                     iradius < extraction_params.num_extraction_radii;
                     ++iradius)
                {
                    std::string radius_name =
                        "extraction_radii[" + std::to_string(iradius) + "]";
                    double radius = extraction_params.extraction_radii[iradius];
                    if (idir == 0)
                        check_parameter(radius_name, radius, radius >= 0.0,
                                        "must be >= 0.0");
                    check_parameter(
                        radius_name, radius,
                        (center_in_dir - radius >=
                         reflective_domain_lo[idir]) &&
                            (center_in_dir + radius <=
                             reflective_domain_hi[idir]),
                        "extraction sphere must lie within the computational "
                        "domain after applying reflective symmetry");
                }
            }
            for (int imode = 0; imode < extraction_params.num_modes; ++imode)
            {
                auto &mode = extraction_params.modes[imode];
                int l = mode.first;
                int m = mode.second;
                std::string mode_name = "modes[" + std::to_string(imode) + "]";
                std::string value_str = "(" + std::to_string(mode.first) +
                                        ", " + std::to_string(mode.second) +
                                        ")";
                check_parameter(
                    mode_name, value_str, (l >= 2) && (abs(m) <= l),
                    "l must be >= 2 and m must satisfy -l <= m <= l");
            }
        }
    }

  public:
    double sigma; // Kreiss-Oliger dissipation parameter
    int nan_check;

    // Collection of parameters necessary for the extraction
    bool activate_extraction;
    SphericalExtraction::params_t extraction_params;

    std::string data_path;
};

#endif /* FIXEDBGSIMULATIONPARAMETERSBASE_HPP_ */
