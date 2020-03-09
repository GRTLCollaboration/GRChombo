/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHOMBOPARAMETERS_HPP_
#define CHOMBOPARAMETERS_HPP_

// General includes
#include "BoundaryConditions.hpp"
#include "GRParmParse.hpp"
#include "UserVariables.hpp"
#include <algorithm>

class ChomboParameters
{
  public:
    ChomboParameters(GRParmParse &pp) { read_params(pp); }

    void read_params(GRParmParse &pp)
    {
        pp.load("verbosity", verbosity, 0);
        // Grid setup
        pp.load("L", L, 1.0);
        pp.load("regrid_threshold", regrid_threshold, 0.5);
        pp.load("num_ghosts", num_ghosts, 3);
        pp.load("tag_buffer_size", tag_buffer_size, 3);
        pp.load("dt_multiplier", dt_multiplier, 0.25);
        pp.load("fill_ratio", fill_ratio, 0.7);

        // Setup the grid size
        ivN = IntVect::Unit;
        int max_N = 0;
        for (int dir = 0; dir < CH_SPACEDIM; ++dir)
        {
            char dir_str[20];
            sprintf(dir_str, "N%d", dir + 1);
            int N;
            pp.load(dir_str, N);
            N_vect.push_back(N);
            ivN[dir] = N - 1;
            max_N = max(N, max_N);
        }
        coarsest_dx = L / max_N;

        pp.load("center", center,
                {0.5 * N_vect[0] * coarsest_dx, 0.5 * N_vect[1] * coarsest_dx,
                 0.5 * N_vect[2] * coarsest_dx}); // default to center

        // Periodicity and boundaries
        pp.load("isPeriodic", isPeriodic, {true, true, true});
        int bc = BoundaryConditions::STATIC_BC;
        pp.load("hi_boundary", boundary_params.hi_boundary, {bc, bc, bc});
        pp.load("lo_boundary", boundary_params.lo_boundary, {bc, bc, bc});
        // set defaults, then override them where appropriate
        boundary_params.vars_parity.fill(BoundaryConditions::EVEN);
        boundary_params.vars_asymptotic_values.fill(0.0);
        boundary_params.is_periodic.fill(true);
        boundary_params.extrapolation_order = 0;
        nonperiodic_boundaries_exist = false;
        boundary_solution_enforced = false;
        FOR1(idir)
        {
            if (isPeriodic[idir] == false)
            {
                nonperiodic_boundaries_exist = true;
                boundary_params.is_periodic[idir] = false;

                // read in relevent params - note that no defaults are set so as
                // to force the user to specify them where the relevant BCs are
                // selected
                if ((boundary_params.hi_boundary[idir] ==
                     BoundaryConditions::REFLECTIVE_BC) ||
                    (boundary_params.lo_boundary[idir] ==
                     BoundaryConditions::REFLECTIVE_BC))
                {
                    boundary_solution_enforced = true;
                    pp.load("vars_parity", boundary_params.vars_parity);

                    if ((boundary_params.hi_boundary[idir] ==
                         BoundaryConditions::REFLECTIVE_BC) &&
                        (boundary_params.lo_boundary[idir] !=
                         BoundaryConditions::REFLECTIVE_BC))
                    {
                        center[idir] = N_vect[idir] * coarsest_dx;
                    }

                    if ((boundary_params.lo_boundary[idir] ==
                         BoundaryConditions::REFLECTIVE_BC) &&
                        (boundary_params.hi_boundary[idir] !=
                         BoundaryConditions::REFLECTIVE_BC))
                    {
                        center[idir] = 0;
                    }
                }
                if ((boundary_params.hi_boundary[idir] ==
                     BoundaryConditions::EXTRAPOLATING_BC) ||
                    (boundary_params.lo_boundary[idir] ==
                     BoundaryConditions::EXTRAPOLATING_BC))
                {
                    boundary_solution_enforced = true;
                    pp.load("extrapolation_order",
                            boundary_params.extrapolation_order, 1);
                }
                if ((boundary_params.hi_boundary[idir] ==
                     BoundaryConditions::SOMMERFELD_BC) ||
                    (boundary_params.lo_boundary[idir] ==
                     BoundaryConditions::SOMMERFELD_BC))
                {
                    pp.load("vars_asymptotic_values",
                            boundary_params.vars_asymptotic_values);
                }
            }
        }
        pout() << "Center has been set to: ";
        FOR1(idir) { pout() << center[idir] << " "; }
        pout() << endl;

        if (nonperiodic_boundaries_exist)
        {
            // write out boundary conditions where non periodic - useful for
            // debug
            BoundaryConditions::write_boundary_conditions(boundary_params);
        }

        // Misc
        pp.load("ignore_checkpoint_name_mismatch",
                ignore_checkpoint_name_mismatch, false);

        pp.load("max_level", max_level, 0);
        // the reference ratio is hard coded to 2 on all levels
        // in principle it can be set to other values, but this is
        // not recommended since we do not test GRChombo with other
        // refinement ratios - use other values at your own risk
        ref_ratios.resize(max_level + 1);
        ref_ratios.assign(2);
        // read in frequency of regrid on each levels, needs
        // max_level + 1 entries (although never regrids on max_level+1)
        pp.getarr("regrid_interval", regrid_interval, 0, max_level + 1);

        // time stepping outputs and regrid data
        pp.load("checkpoint_interval", checkpoint_interval, 1);
        pp.load("chk_prefix", checkpoint_prefix);
        pp.load("plot_interval", plot_interval, 0);
        pp.load("plot_prefix", plot_prefix);
        pp.load("stop_time", stop_time, 1.0);
        pp.load("max_steps", max_steps, 1000000);
        pp.load("write_plot_ghosts", write_plot_ghosts, false);

        // load vars to write to plot files
        pp.load("num_plot_vars", num_plot_vars, 0);
        std::vector<std::string> plot_var_names(num_plot_vars, "");
        pp.load("plot_vars", plot_var_names, num_plot_vars, plot_var_names);
        for (std::string var_name : plot_var_names)
        {
            int var = variable_name_to_enum(var_name);
            if (var >= 0 && var < NUM_VARS)
            {
                plot_vars.push_back(var);
            }
        }
        num_plot_vars = plot_vars.size();

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
    }

    /// Takes a string and returns the variable enum number if the string
    /// matches one of those in UserVariables::variable_names, or returns -1
    /// otherwise
    int variable_name_to_enum(const std::string &a_var_name)
    {
        using namespace UserVariables;

        // std::find did not work very well with the char const* array type of
        // UserVariables::variable_names so here convert to a
        // std::array of std::strings first. This is quite inefficient but this
        // function isn't used much so doesn't matter.
        std::array<std::string, NUM_VARS> variable_names_array;
        for (int ivar = 0; ivar < NUM_VARS; ++ivar)
        {
            variable_names_array[ivar] = variable_names[ivar];
        }

        const auto var_name_it =
            std::find(variable_names_array.begin(), variable_names_array.end(),
                      a_var_name);

        int var = std::distance(variable_names_array.begin(), var_name_it);
        if (var != NUM_VARS)
            return var;
        else
        {
            pout() << "Variable with name " << a_var_name << " not found."
                   << std::endl;
            return -1;
        }
    }

    // General parameters
    int verbosity;
    std::vector<int> N_vect;
    double L;                               // Physical sidelength of the grid
    std::array<double, CH_SPACEDIM> center; // grid center
    IntVect ivN;                 // The number of grid cells in each dimension
    double coarsest_dx;          // The coarsest resolution
    int max_level;               // the max number of regriddings to do
    int num_ghosts;              // must be at least 3 for KO dissipation
    int tag_buffer_size;         // Amount the tagged region is grown by
    Vector<int> ref_ratios;      // ref ratios between levels
    Vector<int> regrid_interval; // steps between regrid at each level
    int max_steps;
    bool ignore_checkpoint_name_mismatch;   // ignore mismatch of variable names
                                            // between restart file and program
    double dt_multiplier, stop_time;        // The Courant factor and stop time
    int checkpoint_interval, plot_interval; // Steps between outputs
    int max_grid_size, block_factor;        // max and min box sizes
    double fill_ratio; // determines how fussy the regridding is about tags
    std::string checkpoint_prefix, plot_prefix; // naming of files
    bool write_plot_ghosts;
    int num_plot_vars;
    std::vector<int> plot_vars; // vars to write to plot file

    // Boundary conditions
    std::array<bool, CH_SPACEDIM> isPeriodic;     // periodicity
    BoundaryConditions::params_t boundary_params; // set boundaries in each dir
    bool nonperiodic_boundaries_exist;
    bool boundary_solution_enforced;

    // For tagging
    double regrid_threshold;
};

#endif /* CHOMBOPARAMETERS_HPP_ */
