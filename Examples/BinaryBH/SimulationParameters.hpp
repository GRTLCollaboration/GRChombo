/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"
#include "BoostedBH.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        read_params(pp);
        check_params();
    }

    /// Read parameters from the parameter file
    void read_params(GRParmParse &pp)
    {
        // Initial data
        pp.load("massA", bh1_params.mass);
        pp.load("momentumA", bh1_params.momentum);
        pp.load("massB", bh2_params.mass);
        pp.load("momentumB", bh2_params.momentum);

        // Get the centers of the BHs either explicitly or as
        // an offset (not both, or they will be offset from center
        // provided)
        std::array<double, CH_SPACEDIM> centerA, centerB;
        std::array<double, CH_SPACEDIM> offsetA, offsetB;
        pp.load("centerA", centerA, center);
        pp.load("centerB", centerB, center);
        pp.load("offsetA", offsetA, {0.0, 0.0, 0.0});
        pp.load("offsetB", offsetB, {0.0, 0.0, 0.0});
        FOR1(idir)
        {
            bh1_params.center[idir] = centerA[idir] + offsetA[idir];
            bh2_params.center[idir] = centerB[idir] + offsetB[idir];
        }

        // Do we want Weyl extraction, puncture tracking and constraint norm
        // calculation?
        pp.load("activate_extraction", activate_extraction, false);
        pp.load("track_punctures", track_punctures, false);
        pp.load("puncture_tracking_level", puncture_tracking_level, max_level);
        pp.load("calculate_constraint_norms", calculate_constraint_norms,
                false);
    }

    void check_params()
    {
        warn_parameter("massA", bh1_params.mass, bh1_params.mass >= 0,
                       "should be >= 0");
        warn_parameter("massB", bh2_params.mass, bh2_params.mass >= 0,
                       "should be >= 0");
        warn_array_parameter(
            "momentumA", bh1_params.momentum,
            std::sqrt(ArrayTools::norm2(bh1_params.momentum)) <
                0.3 * bh1_params.mass,
            "approximation used for boosted BH only valid for small boosts");
        warn_array_parameter(
            "momentumB", bh2_params.momentum,
            std::sqrt(ArrayTools::norm2(bh2_params.momentum)) <
                0.3 * bh1_params.mass,
            "approximation used for boosted BH only valid for small boosts");
        FOR1(idir)
        {
            std::string nameA = "centerA[" + std::to_string(idir) + "]";
            std::string nameB = "centerB[" + std::to_string(idir) + "]";
            double center_A_dir = bh1_params.center[idir];
            double center_B_dir = bh2_params.center[idir];
            warn_parameter(nameA, center_A_dir,
                           (center_A_dir >= 0.0) &&
                               (center_A_dir <= (ivN[idir] + 1) * coarsest_dx),
                           "should be within the computational domain");
            warn_parameter(nameB, center_B_dir,
                           (center_B_dir >= 0.0) &&
                               (center_B_dir <= (ivN[idir] + 1) * coarsest_dx),
                           "should be within the computational domain");
        }
        check_parameter("puncture_tracking_level", puncture_tracking_level,
                        (puncture_tracking_level >= 0) &&
                            (puncture_tracking_level <= max_level),
                        "must be between 0 and max_level (inclusive)");
    }

    // Initial data
    bool activate_extraction, track_punctures, calculate_constraint_norms;
    int puncture_tracking_level;
    // Collection of parameters necessary for initial conditions
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;
};

#endif /* SIMULATIONPARAMETERS_HPP */
