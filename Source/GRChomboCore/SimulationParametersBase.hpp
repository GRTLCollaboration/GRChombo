/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERSBASE_HPP_
#define SIMULATIONPARAMETERSBASE_HPP_

// General includes
#include "CCZ4.hpp"
#include "GRParmParse.hpp"

struct extraction_params_t
{
    double extraction_radius;
    std::array<double, CH_SPACEDIM> extraction_center;
    int num_points_phi;
    int num_points_theta;
    int extraction_level;
};

class SimulationParametersBase
{
  public:
    SimulationParametersBase(GRParmParse &pp) { readBaseParams(pp); }

    void readBaseParams(GRParmParse &pp)
    {
        // The read parameters code defined below
        auto_read_base_params(pp);

        // Fill in the ccz4Parameters
        ccz4_params.kappa1 = kappa1;
        ccz4_params.kappa2 = kappa2;
        ccz4_params.kappa3 = kappa3;
        ccz4_params.shift_Gamma_coeff = shift_Gamma_coeff;
        ccz4_params.shift_advec_coeff = shift_advec_coeff;
        ccz4_params.eta = eta;
        ccz4_params.lapse_power = lapse_power;
        ccz4_params.lapse_coeff = lapse_coeff;
        ccz4_params.lapse_advec_coeff = lapse_advec_coeff;

        // Fill in the params
        extraction_params.extraction_radius = extraction_radius;
        extraction_params.extraction_center = extraction_center;
        extraction_params.num_points_phi = num_points_phi;
        extraction_params.num_points_theta = num_points_theta;
        extraction_params.extraction_level = extraction_level;
    }

    void auto_read_base_params(GRParmParse &pp)
    {
        pp.load("verbosity", verbosity, 0);
        // Grid setup
        pp.load("L", L, 1.0);
        pp.load("center", center,
                {0.5 * L, 0.5 * L, 0.5 * L}); // default to center
        pp.load("regrid_threshold", regrid_threshold, 0.5);
        pp.load("isPeriodic", isPeriodic, {true, true, true});
        pp.load("num_ghosts", num_ghosts, 3);
        pp.load("tag_buffer_size", tag_buffer_size, 3);
        pp.load("dt_multiplier", dt_multiplier, 0.25);
        pp.load("fill_ratio", fill_ratio, 0.7);

        // Lapse evolution
        pp.load("lapse_advec_coeff", lapse_advec_coeff, 1.0);
        pp.load("lapse_coeff", lapse_coeff, 2.0);
        pp.load("lapse_power", lapse_power, 1.0);

        // Shift Evolution
        pp.load("shift_advec_coeff", shift_advec_coeff, 0.0);
        pp.load("shift_Gamma_coeff", shift_Gamma_coeff, 0.75);
        pp.load("eta", eta, 1.0);

        // CCZ4 parameters
        pp.load("formulation", formulation, 0);
        pp.load("kappa1", kappa1, 0.1);
        pp.load("kappa2", kappa2, 0.0);
        pp.load("kappa3", kappa3, 1.0);
        pp.load("covariantZ4", covariantZ4, 1);

        // Dissipation
        pp.load("sigma", sigma, 0.1);

        // Misc
        pp.load("nan_check", nan_check, 1);
        pp.load("ignore_checkpoint_name_mismatch",
                ignore_checkpoint_name_mismatch, false);

        // Setup the grid size
        ivN = IntVect::Unit;
        int max_N = 0;
        for (int dir = 0; dir < CH_SPACEDIM; ++dir)
        {
            char dir_str[20];
            sprintf(dir_str, "N%d", dir + 1);
            int N;
            pp.load(dir_str, N);
            ivN[dir] = N - 1;
            max_N = max(ivN[dir], max_N);
        }

        // the reference ratio is hard coded to 2
        // in principle it can be ther values but this is
        // not recommended!
        pp.load("max_level", max_level, 0);
        ref_ratios.resize(max_level + 1);
        ref_ratios.assign(2);
        pp.getarr("regrid_interval", regrid_interval, 0, max_level + 1);

        // time stepping outputs and regrid data
        pp.load("checkpoint_interval", checkpoint_interval, 1);
        pp.load("chk_prefix", checkpoint_prefix);
        pp.load("plot_interval", plot_interval, 0);
        pp.load("plot_prefix", plot_prefix);
        pp.load("stop_time", stop_time, 1.0);
        pp.load("max_steps", max_steps, 1000000);

        // alias the weird chombo names to something more descriptive
        // for these box params, and default to some reasonable values
        if (pp.contains("max_grid_size"))
        {
            pp.load("max_grid_size", max_grid_size);
        }
        else
        {
            pp.load("max_box_size", max_grid_size, 64);
        }
        if (pp.contains("block_factor"))
        {
            pp.load("block_factor", block_factor);
        }
        else
        {
            pp.load("min_box_size", block_factor, 8);
        }

        // extraction params
        double dx_scalar = L / max_N;
        dx.fill(dx_scalar);
        origin.fill(dx_scalar / 2.0);

        // Extraction params
        pp.load("extraction_level", extraction_level, 3);
        pp.load("extraction_radius", extraction_radius, 120.0);
        pp.load("num_points_phi", num_points_phi, 8);
        pp.load("num_points_theta", num_points_theta, 16);
        pp.load("extraction_center", extraction_center,
                {0.5 * L, 0.5 * L, 0.5 * L});
    }

    // General parameters
    int verbosity;
    double L;                               // Physical sidelength of the grid
    std::array<double, CH_SPACEDIM> center; // grid center
    IntVect ivN;         // The number of grid cells in each dimension
    int max_level;       // the max number of regriddings to do
    int num_ghosts;      // must be at least 3 for KO dissipation
    int tag_buffer_size; // Amount the tagged region is grown by
    std::array<bool, CH_SPACEDIM> isPeriodic; // periodicity
    Vector<int> ref_ratios;                   // ref ratios between levels
    Vector<int> regrid_interval; // steps between regrid at each level
    int nan_check, max_steps;
    bool ignore_checkpoint_name_mismatch;   // ignore mismatch of variable names
                                            // between restart file and program
    double dt_multiplier, stop_time;        // The Courant factor and stop time
    int checkpoint_interval, plot_interval; // Steps between outputs
    int max_grid_size, block_factor;        // max and min box sizes
    double fill_ratio; // determines how fussy the regridding is about tags
    std::string checkpoint_prefix, plot_prefix; // naming of files

    // Extraction parameters
    std::array<double, CH_SPACEDIM> origin,
        dx; // location of coarsest origin and dx
    double extraction_radius;
    std::array<double, CH_SPACEDIM> extraction_center;
    int num_points_phi, num_points_theta, extraction_level;

    // For tagging
    double regrid_threshold;
    // Lapse evolution
    double lapse_power, lapse_coeff, lapse_advec_coeff;
    // Shift Evolution
    double shift_advec_coeff, shift_Gamma_coeff, eta;
    // CCZ4 parameters
    int formulation;
    double kappa1, kappa2, kappa3;
    int covariantZ4;
    // Kreiss Oliger dissipation
    double sigma;
    // Collection of parameters necessary for the CCZ4 RHS and extraction
    CCZ4::params_t ccz4_params;
    extraction_params_t extraction_params;
};

#endif /* SIMULATIONPARAMETERSBASE_HPP_ */
