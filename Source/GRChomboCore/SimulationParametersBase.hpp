/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERSBASE_HPP_
#define SIMULATIONPARAMETERSBASE_HPP_

// General includes
#include "BoundaryConditions.hpp"
#include "CCZ4.hpp"
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"
#include "SphericalExtraction.hpp"

// add this type alias here for backwards compatibility
using extraction_params_t = SphericalExtraction::params_t;

class SimulationParametersBase : public ChomboParameters
{
  public:
    SimulationParametersBase(GRParmParse &pp) : ChomboParameters(pp)
    {
        read_params(pp);
    }

  private:
    void read_params(GRParmParse &pp)
    {
        // Lapse evolution
        pp.load("lapse_advec_coeff", ccz4_params.lapse_advec_coeff, 1.0);
        pp.load("lapse_coeff", ccz4_params.lapse_coeff, 2.0);
        pp.load("lapse_power", ccz4_params.lapse_power, 1.0);

        // Shift Evolution
        pp.load("shift_advec_coeff", ccz4_params.shift_advec_coeff, 0.0);
        pp.load("shift_Gamma_coeff", ccz4_params.shift_Gamma_coeff, 0.75);
        pp.load("eta", ccz4_params.eta, 1.0);

        // CCZ4 parameters
        pp.load("formulation", formulation, 0);
        pp.load("kappa1", ccz4_params.kappa1, 0.1);
        pp.load("kappa2", ccz4_params.kappa2, 0.0);
        pp.load("kappa3", ccz4_params.kappa3, 1.0);

        // Dissipation
        pp.load("sigma", sigma, 0.1);

        // Nan Check and min chi and lapse values
        pp.load("nan_check", nan_check, 1);
        pp.load("min_chi", min_chi, 1e-4);
        pp.load("min_lapse", min_lapse, 1e-4);

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
        pp.load("num_points_theta", extraction_params.num_points_theta, 5);
        if (extraction_params.num_points_theta % 2 == 0)
        {
            extraction_params.num_points_theta += 1;
            pout() << "Parameter: num_points_theta incompatible with Simpson's "
                   << "rule so increased by 1.\n";
        }
        pp.load("extraction_center", extraction_params.center, center);

        check_radii();

        if (pp.contains("modes"))
        {
            pp.load("num_modes", extraction_params.num_modes);
            std::vector<int> extraction_modes_vect(2 *
                                                   extraction_params.num_modes);
            pp.load("modes", extraction_modes_vect,
                    2 * extraction_params.num_modes);
            extraction_params.modes.resize(extraction_params.num_modes);
            for (int i = 0; i < extraction_params.num_modes; ++i)
            {
                extraction_params.modes[i].first = extraction_modes_vect[2 * i];
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

        pp.load("write_extraction", extraction_params.write_extraction, false);
    }

    void check_radii()
    {
        // check center of extraction and extraction radii are compatible with
        // box size
        for (auto &radius : extraction_params.extraction_radii)
        {
            std::array<double, CH_SPACEDIM> axis_distance_to_boundary;

            // upper boundary
            FOR1(i)
            {
                axis_distance_to_boundary[i] =
                    (boundary_params.hi_boundary[i] ==
                             BoundaryConditions::REFLECTIVE_BC
                         ? 2.
                         : 1.) *
                        (ivN[i] + 1) * coarsest_dx -
                    extraction_params.center[i];
            }
            if (radius >= *std::min_element(axis_distance_to_boundary.begin(),
                                            axis_distance_to_boundary.end()))
                MayDay::Error(
                    "Extraction radii go being the box's upper boundary");

            // lower boundary
            FOR1(i)
            {
                axis_distance_to_boundary[i] =
                    extraction_params.center[i] +
                    (boundary_params.lo_boundary[i] ==
                             BoundaryConditions::REFLECTIVE_BC
                         ? 1.
                         : 0.) *
                        (ivN[i] + 1) * coarsest_dx;
            }
            if (radius >= *std::min_element(axis_distance_to_boundary.begin(),
                                            axis_distance_to_boundary.end()))
                MayDay::Error(
                    "Extraction radii go being the box's lower boundary");
        }
    }

  public:
    double sigma; // Kreiss-Oliger dissipation parameter

    int nan_check;

    double min_chi, min_lapse;

    int formulation; // Whether to use BSSN or CCZ4

    // Collection of parameters necessary for the CCZ4 RHS and extraction
    CCZ4::params_t ccz4_params;
    SphericalExtraction::params_t extraction_params;
};

#endif /* SIMULATIONPARAMETERSBASE_HPP_ */
