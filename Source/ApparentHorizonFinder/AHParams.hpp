/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHPARAMS_HPP_
#define _AHPARAMS_HPP_

#include "AHInitialGuess.hpp"
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"

#include "SPMD.H"

// Chombo namespace
#include "UsingNamespace.H"

// PETSc Solver related params
// not very relevant to change, but this allows us to in case we want
// possible the only relevant one is 'snes_maxit' ('AH_SNES_max_iterations')
// which for highly ellipsoidal AHs sometimes needs to be increased as
// convergence is rather slow
// (as of 06/2021, the code outputs this as a recommendation when that happens)
struct PETSc_params
{
    /*
    (defaults show according to a certain PETSc version,
    can't assure it is the same for all)
    atol    - absolute convergence tolerance (default 1e-50)
    rtol    - relative convergence tolerance (default 1e-08)
    stol    - convergence tolerance in terms of the norm of the change in the
            solution between steps (default 1e-08)
    maxit   - maximum number of iterations (default 50)
    maxf    - maximum number of function evaluations (default 10000)
    divtol  - divergence tolerance (default 10000)
    */
    PetscReal snes_atol = -1., snes_rtol = -1., snes_stol = -1.;
    PetscInt snes_maxit = -1, snes_maxf = -1;
    PetscReal snes_divtol = -1.;

    /*
    rtol    - the relative convergence tolerance (default 1e-05)
    abstol  - the absolute convergence tolerance (default 1e-50)
    dtol    - the divergence tolerance (default 10000)
    maxits  - maximum number of iterations (default 10000)
    */
    PetscReal ksp_rtol = -1., ksp_abstol = -1., ksp_dtol = -1.;
    PetscInt ksp_maxits = -1;

    void read_params(GRParmParse &pp)
    {
        // variables to prevent ParmParse to complain about types (only happens
        // in some clusters)
        int tmp_int;
        float tmp_float;
        if (pp.contains("AH_SNES_absolute_tolerance"))
        {
            pp.load("AH_SNES_absolute_tolerance", tmp_float);
            snes_atol = tmp_float;
        }
        if (pp.contains("AH_SNES_relative_tolerance"))
        {
            pp.load("AH_SNES_relative_tolerance", tmp_float);
            snes_rtol = tmp_float;
        }
        if (pp.contains("AH_SNES_step_change_tolerance"))
        {
            pp.load("AH_SNES_step_change_tolerance", tmp_float);
            snes_stol = tmp_float;
        }
        if (pp.contains("AH_SNES_max_iterations"))
        {
            pp.load("AH_SNES_max_iterations", tmp_int);
            snes_maxit = tmp_int;
        }
        if (pp.contains("AH_SNES_max_evaluations"))
        {
            pp.load("AH_SNES_max_evaluations", tmp_int);
            snes_maxf = tmp_int;
        }
        if (pp.contains("AH_SNES_divergence_tolerance"))
        {
            pp.load("AH_SNES_divergence_tolerance", tmp_float);
            snes_divtol = tmp_float;
        }

        if (pp.contains("AH_KSP_relative_tolerance"))
        {
            pp.load("AH_KSP_relative_tolerance", tmp_float);
            ksp_rtol = tmp_float;
        }
        if (pp.contains("AH_KSP_absolute_tolerance"))
        {
            pp.load("AH_KSP_absolute_tolerance", tmp_float);
            ksp_abstol = tmp_float;
        }
        if (pp.contains("AH_KSP_divergence_tolerance"))
        {
            pp.load("AH_KSP_divergence_tolerance", tmp_float);
            ksp_dtol = tmp_float;
        }
        if (pp.contains("AH_KSP_max_evaluations"))
        {
            pp.load("AH_KSP_max_evaluations", tmp_int);
            ksp_maxits = tmp_int;
        }
    }
};

// prepend with 'AH_' in params file
template <class AHFunction> struct AHParams_t
{
    int num_ranks; //!< number of ranks for PETSc sub-communicator (default
                   //!< 0, which is 'all')

    int num_points_u; //!< number of points for 2D coordinate grid
#if CH_SPACEDIM == 3
    int num_points_v; //!< number of points for 2D coordinate grid
#endif
    int solve_interval;  //!< same as checkpoint_interval, for
                         //!< ApparentHorizon::solve (default 1)
    int print_interval;  //!< same as solve_interval, but for prints (default
                         //!< 1)
    bool track_center;   //!< whether or not to update the center
                         //!< (set to false if you know it won't move)
                         //!< (default true)
    bool predict_origin; //!< whether or not to estimate where the next
                         //!< origin will be at (default = track_center)

    int level_to_run; // if negative, it will count backwards (e.g. -1 is
                      // 'last level') (default 0)

    double start_time;   //!< time after which AH can start (default 0.)
                         //!< Useful for ScalarField collapse
    double give_up_time; //!< stop if at this time nothing was found
                         //!< (<0 to never, which is default)
                         //!< Useful for ScalarField collapse

    bool allow_re_attempt;    //!< re-attempt with initial guess if
                              //!< previous convergence failed (default false)
    int max_fails_after_lost; //!< number of time steps to try again after
                              //!< (-1 to never) the AH was lost
                              //!< (default is 0)

    int verbose; //!< print information during execution (default is 1)

    bool print_geometry_data; //!< print metric and extrinsic
                              //!< curvature of the AHs (default false)

    bool re_solve_at_restart; //!< whether or not to re-run even if AH
                              //!< already exists (useful in order to be
                              //!< able to provide an initial guess and
                              //!< re-run the AH on that time step)
                              //!< (default false)

    bool stop_if_max_fails; //! breaks the run if AH doesn't converge
                            //! 'max_fails_after_lost' times or if
                            //! 'give_up_time' is reached without
                            //! convergence (default is 'false')

    std::map<std::string, std::tuple<int, VariableType, int>>
        extra_vars;     //! extra vars to write to coordinates file (<enum,
                        //! evolution or diagnostic, int for local|d1|d2>)
    int num_extra_vars; // total number of extra vars (!=extra_vars.size()
                        // as derivative count for multiple vars)

    int extra_contain_diagnostic; // not a parameter (set internally);
                                  // counts how many

    std::string stats_path = "",
                stats_prefix = "stats_AH"; //!< name for stats file with
                                           //!< area, spin and AH origin/center
    std::string coords_path = "",
                coords_prefix = "coords_AH"; //!< name for coords file with AH
                                             //!< coordinates at each time step

    //! mergers will be searched when distance between 'parent' BHs is
    //! distance < merger_search_factor * 4. * (AH_initial_guess_1 +
    //! AH_initial_guess_2) should be roughly '2M=2(m1+m2)' for initial
    //! guess at m/2 (set to non-positive to 'always search')
    double merger_search_factor; // see note above (default is 1)
    //! initial guess for merger is 'merger_pre_factor * 4. *
    //! (AH_initial_guess_1 + AH_initial_guess_2)'
    //! set to something bigger to avoid finding the inner AH
    double merger_pre_factor; // see note above (default to 1.)

    typename AHFunction::params func_params;

    PETSc_params petsc_params;

    enum verbose_level
    {
        NONE,
        MIN,  // minimal
        SOME, // a bit technical
        MORE  // debug
    };

    void read_params(GRParmParse &pp, const ChomboParameters &a_p);
};

template <class AHFunction>
void AHParams_t<AHFunction>::read_params(GRParmParse &pp,
                                         const ChomboParameters &a_p)
{
    pp.load("AH_num_ranks", num_ranks, 0); // 0 means "all"

    pp.load("AH_num_points_u", num_points_u);
#if CH_SPACEDIM == 3
    pp.load("AH_num_points_v", num_points_v);
#endif
    int actual_num_ranks = (num_ranks != 0) ? num_ranks : numProc();

    // Let's enforce the below as Katy had problems otherwise
    ChomboParameters::check_parameter("AH_num_ranks", num_ranks,
                                      num_ranks <= numProc(),
                                      "must be <= total number of ranks");

    ChomboParameters::check_parameter("AH_num_points_u", num_points_u,
                                      num_points_u > 0, "must be >0");
#if CH_SPACEDIM == 3
    ChomboParameters::check_parameter("AH_num_points_v", num_points_v,
                                      num_points_v > 0, "must be >0");
    // if box division ends up with side length less than 3, PETSc will
    // complain - estimate side length using sqrt(num_ranks)
    ChomboParameters::check_parameter(
        "AH_num_points_u", num_points_u,
        num_points_u / sqrt(actual_num_ranks) >= 3,
        "(num_points_u/sqrt(num_ranks)) must be >3"); // make sure for size 'u'
    ChomboParameters::check_parameter(
        "AH_num_points_v", num_points_v,
        num_points_v / sqrt(actual_num_ranks) >= 3,
        "(num_points_v/sqrt(num_ranks)) must be >3"); // make sure for size 'v'

    // Not sure if PETSc requires the below but let's just enforce it just
    // in case
    ChomboParameters::check_parameter(
        "AH_num_ranks", num_ranks,
        (num_points_u * num_points_v) % num_ranks == 0,
        "must divide AH_num_points_u * AH_num_points_v");
#elif CH_SPACEDIM == 2
    ChomboParameters::check_parameter(
        "AH_num_points_u", num_points_u, num_points_u / actual_num_ranks >= 3,
        "(num_points_u/size) must be >0"); // make sure for size 'u'
#endif

    pp.load("AH_solve_interval", solve_interval, 1);
    pp.load("AH_print_interval", print_interval, 1);
    // sanity checks
    ChomboParameters::check_parameter("AH_solve_interval", solve_interval,
                                      solve_interval > 0, "must be >0");
    ChomboParameters::check_parameter("AH_print_interval", print_interval,
                                      print_interval > 0, "must be >0");
    pp.load("AH_track_center", track_center, true);
    pp.load("AH_predict_origin", predict_origin, track_center);
    // can't predict if center is not being tracked
    if (predict_origin)
    {
        ChomboParameters::check_parameter(
            "AH_track_center", track_center, track_center == 1,
            "can't predict if center is not being tracked");
    }
    pp.load("AH_level_to_run", level_to_run, 0);
    ChomboParameters::check_parameter("AH_level_to_run", level_to_run,
                                      level_to_run <= a_p.max_level,
                                      "must be <= max_level");
    ChomboParameters::check_parameter("AH_level_to_run", level_to_run,
                                      level_to_run >= 0, "must be >= 0");
    ChomboParameters::check_parameter("AH_level_to_run", level_to_run,
                                      level_to_run > -(a_p.max_level + 1),
                                      "must be > -(max_level+1)");
    if (level_to_run < 0) // if negative, count backwards
        level_to_run += a_p.max_level + 1;

    pp.load("AH_start_time", start_time, 0.0);
    pp.load("AH_give_up_time", give_up_time, -1.0);

    pp.load("AH_allow_re_attempt", allow_re_attempt, false);
    pp.load("AH_max_fails_after_lost", max_fails_after_lost, 0);
    pp.load("AH_stop_if_max_fails", stop_if_max_fails, false);

    pp.load("AH_verbose", verbose, (int)MIN);
    pp.load("AH_print_geometry_data", print_geometry_data, false);
    pp.load("AH_re_solve_at_restart", re_solve_at_restart, false);
    // load vars to write to coord files
    num_extra_vars = 0;
    extra_contain_diagnostic = 0;

    int AH_num_extra_vars;
    pp.load("AH_num_extra_vars", AH_num_extra_vars, 0);
    if (AH_num_extra_vars > 0)
    {
        std::vector<std::string> AH_extra_var_names(AH_num_extra_vars, "");
        pp.load("AH_extra_vars", AH_extra_var_names, AH_num_extra_vars);
        for (const std::string &full_name : AH_extra_var_names)
        {
            std::string var_name = full_name;

            // variable names might start with "d1_" or "d2_" to indicate
            // the user wants derivatives
            int der_type = 0;
            std::string der = var_name.substr(0, 3);
            if (der == "d1_")
            {
                der_type = 1;
                var_name = var_name.substr(3);
            }
            else if (der == "d2_")
            {
                der_type = 2;
                var_name = var_name.substr(3);
            }

            // first assume extra_var is a normal evolution var
            int var = UserVariables::variable_name_to_enum(var_name);
            VariableType var_type = VariableType::evolution;
            if (var < 0)
            {
                // if not an evolution var check if it's a diagnostic var
                var = DiagnosticVariables::variable_name_to_enum(var_name);
                if (var < 0)
                {
                    // it's neither :(
                    pout() << "Variable with name " << var_name
                           << " not found.\n";
                }
                else
                {
                    var_type = VariableType::diagnostic;
                    ++extra_contain_diagnostic;
                }
            }
            if (var >= 0)
            {
                extra_vars[full_name] =
                    std::tuple<int, VariableType, int>(var, var_type, der_type);
                if (der_type == 0)
                    num_extra_vars += 1;
                else if (der_type == 1)
                    num_extra_vars += CH_SPACEDIM;
                else // if (der_type == 2)
                    num_extra_vars += CH_SPACEDIM * (CH_SPACEDIM + 1) / 2;
            }
        }
    }

    stats_path = a_p.data_path;

    if (pp.contains("AH_coords_subpath"))
    {
        pp.load("AH_coords_subpath", coords_path);
        if (!coords_path.empty() && coords_path.back() != '/')
            coords_path += "/";
        if (a_p.output_path != "./" && !a_p.output_path.empty())
            coords_path = a_p.output_path + coords_path;
    }
    else
        coords_path = stats_path;

    pp.load("AH_stats_prefix", stats_prefix, std::string("stats_AH"));
    pp.load("AH_coords_prefix", coords_prefix, std::string("coords_AH"));

    pp.load("AH_merger_search_factor", merger_search_factor, 1.);
    pp.load("AH_merger_pre_factor", merger_pre_factor, 1.);

    func_params.read_params(pp);
    petsc_params.read_params(pp);
}

#endif /* _AHPARAMS_HPP_ */
