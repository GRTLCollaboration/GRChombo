/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef USE_AHFINDER

#include "AHFinder.hpp"

#include "AHInterpolation.hpp"
#include "ApparentHorizon.hpp"

AHFinder::~AHFinder()
{
    // destroy horizon pointers and finalize PETSc
    for (auto &ah : m_apparent_horizons)
        delete ah;
    PETSc_finalize();
}

int AHFinder::add_ah(const AHSurfaceGeometry &a_coord_system,
                     double a_initial_guess, const AHFinder::params &a_params,
                     const std::string &a_stats, const std::string &a_coords,
                     bool solve_first_step)
{
    return add_ah(a_coord_system, a_initial_guess, a_params,
                  AHFunction::params(), a_stats, a_coords, solve_first_step);
}

int AHFinder::add_ah(const AHSurfaceGeometry &a_coord_system,
                     double a_initial_guess, const AHFinder::params &a_params,
                     const typename AHFunction::params &a_func_params,
                     const std::string &a_stats, const std::string &a_coords,
                     bool solve_first_step)
{
    if (!AHFinder::m_initialized)
        AHFinder::PETSc_initialize(a_params.num_ranks);

    // determine how many AH there are already
    const int num_ah = m_apparent_horizons.size();

    AHInterpolation<AHSurfaceGeometry, AHFunction> interp(a_coord_system,
                                                          m_interpolator);

    m_apparent_horizons.push_back(
        new ApparentHorizon<AHSurfaceGeometry, AHFunction>(
            interp, a_initial_guess, a_params, a_func_params,
            a_stats + std::to_string(num_ah + 1),
            a_coords + std::to_string(num_ah + 1) + "_", solve_first_step));
    m_merger_pairs.push_back({-1, -1});

    return num_ah;
}

int AHFinder::add_ah_merger(int ah1, int ah2, const params &a_params)
{
    CH_assert(a_params.track_center); // center tracking must be on for mergers

    int num_ah = m_apparent_horizons.size();
    CH_assert(ah1 >= 0 && ah1 < num_ah);
    CH_assert(ah2 >= 0 && ah2 < num_ah);
    CH_assert(ah1 != ah2);

    double initial_guess_merger;
    std::array<double, CH_SPACEDIM> origin_merger;
    bool do_solve = solve_merger(ah1, ah2, initial_guess_merger, origin_merger);

    // assume same coordiante system between merging horizons
    AHSurfaceGeometry coord_system =
        m_apparent_horizons[ah1]->get_ah_interp().get_coord_system();
    coord_system.set_origin(origin_merger);

    std::string stats = m_apparent_horizons[ah1]->m_stats;
    std::string coords = m_apparent_horizons[ah1]->m_coords;
    // remove number of AH1 from file names
    stats = stats.substr(0, stats.size() - std::to_string(ah1).size());
    coords = coords.substr(0, coords.size() - std::to_string(ah1).size() - 1);

    auto &function_to_optimize_params = m_apparent_horizons[ah1]->m_func_params;

    int num = add_ah(coord_system, initial_guess_merger, a_params,
                     function_to_optimize_params, stats, coords, do_solve);
    m_merger_pairs[num] = {ah1, ah2};

    return num;
}

void AHFinder::solve(double a_dt, double a_time, double a_restart_time)
{
    CH_TIME("AHFinder::solve");

    enum
    {
        NOT_SOLVED,
        TOO_FAR,
        SKIPPED,
        SOLVED
    };

    // in case this was called at specificPostTimeStep at t=0 due to
    // MultiLevelTask (solve already happened when AH was created)
    if (a_time == 0.)
        return;

    // solve if it has never converged before OR if has already converged in the
    // past this means it will stop solving as soon as it stops converging but
    // keeps trying (with same initial guess) if never found one
    // (e.g. in the BinaryBH case, AH1 and AH2 stop solving once they
    // disappear, leaving the merged AH alone from then on)
    // (e.g. in the ScalarField case, only after a few time steps does an AH
    // form, and from then on it keeps solving forever)

    const int num_ah = m_apparent_horizons.size();

    // skip if AH 0 is not supposed to be solved
    // assumes all AHs have same 'solve_interval'
    if (num_ah == 0 || !m_apparent_horizons[0]->do_solve(a_dt, a_time))
        return;

    // the mess below is done so that AH are searched for iteratively
    // taking into account merger priorities
    // (a merger will only be solved if it's 'parents' have already been solved)
    int ahs_solved = 0;
    // sets whether or not AH[i] has already been solved
    // 0 - not solved; 1 - too far away; 2 - skipped; 3 - solved
    int *ah_solved =
        new int[num_ah]; // with '-g' compilation flag was complaining that int
                         // ah_solved[num_ah] was invalid
    for (int i = 0; i < num_ah; ++i)
        ah_solved[i] = NOT_SOLVED;

    int max_attempts = 100, attempt = 0;
    while (ahs_solved < num_ah)
    {
        if (max_attempts < attempt++)
        {
            pout() << "Reached maximum number of attempts in AH solving loop."
                   << std::endl;
            break;
        }

        for (int i = 0; i < num_ah; ++i)
        {
            if (ah_solved[i] != NOT_SOLVED)
                continue; // continue if already solved

            // if is a merger and has not yet converged
            if (m_merger_pairs[i].first >= 0 &&
                !m_apparent_horizons[i]->has_been_found())
            {
                auto pair = m_merger_pairs[i];

                if (!((ah_solved[pair.first] == SKIPPED ||
                       ah_solved[pair.first] == SOLVED) &&
                      (ah_solved[pair.second] == SKIPPED ||
                       ah_solved[pair.second] == SOLVED)))
                    continue; // continue if one of the 'parents' was not
                              // solved/skipped yet

                if (!(m_apparent_horizons[pair.first]->has_been_found() &&
                      m_apparent_horizons[pair.second]->has_been_found()))
                {
                    ++ahs_solved;
                    ah_solved[i] = SKIPPED;
                    continue; // skip if one of the parents was never found yet
                }

                if (ah_solved[pair.first] == TOO_FAR ||
                    ah_solved[pair.second] == TOO_FAR)
                { // if one of the 'parents' is a merger too far away to be
                  // solved, skip children too
                    if (m_apparent_horizons[i]->m_params.verbose > NONE)
                    {
                        pout() << "Skipping BH #" << i << std::endl;
                    }
                    ++ahs_solved;
                    ah_solved[i] = TOO_FAR;
                }

                double initial_guess_merger;
                std::array<double, CH_SPACEDIM> center_merger;
                bool do_solve =
                    solve_merger(pair.first, pair.second, initial_guess_merger,
                                 center_merger);

                // if distance is too large, ignore this one and move on
                if (!do_solve)
                {
                    ++ahs_solved;
                    ah_solved[i] = TOO_FAR;
                    continue;
                }

                // if it already converged, don't do anything;
                // no need to change the initial guess, as it is the same since
                // the beginning;
                // if one was skipped don't update the center (because it will
                // not be central anymore)
                if (!m_apparent_horizons[i]->get_converged() &&
                    ah_solved[pair.first] == SOLVED &&
                    ah_solved[pair.second] == SOLVED)
                {
                    m_apparent_horizons[i]->set_origin(center_merger);
                    // m_apparent_horizons[i]->set_initial_guess(
                    // initial_guess_merger); // should still be the same
                }
            }

            // solve if it converged last time or if is has never been found in
            // the past
            if (m_apparent_horizons[i]->good_to_go(a_dt, a_time))
            {
                if (m_apparent_horizons[i]->m_params.verbose > NONE)
                    pout() << "Solving AH #" << i << std::endl;

                m_apparent_horizons[i]->solve(a_dt, a_time, a_restart_time);

                ++ahs_solved;
                ah_solved[i] = SOLVED;
            }
            else
            {
                ++ahs_solved;
                ah_solved[i] = SKIPPED;
            }
        }
    }

    delete ah_solved;
}

bool AHFinder::need_diagnostics(double a_dt, double a_time) const
{
    bool out = false;
    for (auto &ah : m_apparent_horizons)
        out |=
            ah->m_params.extra_contain_diagnostic && ah->do_print(a_dt, a_time);
    return out;
}

bool AHFinder::solve_merger(int ah1, int ah2, double &initial_guess_merger,
                            std::array<double, CH_SPACEDIM> &center_merger)
{
    // SKIP if 'parents' not yet close enough
    auto AH1 = m_apparent_horizons[ah1];
    auto AH2 = m_apparent_horizons[ah2];
    auto center1 = AH1->get_center();
    auto center2 = AH2->get_center();

    double initial_guess_sum =
        (AH1->get_initial_guess() + AH2->get_initial_guess());

    // some tests revealed it was about 1.2 * initial_guess_sum, but 1.5 to
    // ensure we don't catch the inner AH
    double merger_pre_factor = std::max(AH1->m_params.merger_pre_factor,
                                        AH2->m_params.merger_pre_factor);
    initial_guess_merger = merger_pre_factor * initial_guess_sum;

    // update center of merged, otherwise it does it by
    // itself in solve
    FOR1(a)
    {
        if (!AH1->get_ah_interp().get_interpolator()->get_boundary_reflective(
                Side::Lo, a))
            center_merger[a] = (center1[a] + center2[a]) / 2.;
        else
            center_merger[a] = 0.;
    }

    // check if distance is too big
    double distance =
        sqrt((center1[0] - center2[0]) * (center1[0] - center2[0]) +
             (center1[1] - center2[1]) * (center1[1] - center2[1])
#if CH_SPACEDIM == 3
             + (center1[2] - center2[2]) * (center1[2] - center2[2])
#endif
        );

    double merger_search_factor = std::max(AH1->m_params.merger_search_factor,
                                           AH2->m_params.merger_search_factor);

    double min_distance = merger_search_factor * 4.0 * initial_guess_sum;

    bool do_solve = false;

    if (AH1->has_been_found() && AH2->has_been_found())
    {
        do_solve = merger_search_factor <= 0. || distance <= min_distance;

        if (do_solve)
            // radius of guess is bigger than AH distance
            CH_assert(min_distance < initial_guess_merger);

        if (AH1->m_params.verbose > NONE)
        {
            pout() << "BHs #" << ah1 << " and #" << ah2
                   << " at distance = " << distance;
            // if distance is too large, ignore this one and move on
            if (!do_solve)
                pout() << " > minimum distance = " << min_distance
                       << ". Skipping solve for merger...";

            pout() << std::endl;
        }
    }

    return do_solve;
}

void AHFinder::params::read_params(GRParmParse &pp, const ChomboParameters &a_p)
{
    pp.load("AH_num_ranks", num_ranks, 0); // 0 means "all"
    pp.load("AH_num_points_u", num_points_u);
#if CH_SPACEDIM == 3
    pp.load("AH_num_points_v", num_points_v);
#endif
    pp.load("AH_solve_interval", solve_interval, 1);
    pp.load("AH_print_interval", print_interval, 1);
    pp.load("AH_track_center", track_center, true);
    pp.load("AH_predict_origin", predict_origin, track_center);
    // can't predict if center is not being tracked
    CH_assert(!(predict_origin && !track_center));

    pp.load("AH_level_to_run", level_to_run, 0);
    CH_assert(level_to_run <= a_p.max_level &&
              level_to_run > -(a_p.max_level + 1));
    if (level_to_run < 0) // if negative, count backwards
        level_to_run += a_p.max_level + 1;

    pp.load("AH_start_time", start_time, 0.0);
    pp.load("AH_give_up_time", give_up_time, -1.0);

    pp.load("AH_merger_search_factor", merger_search_factor, 1.);
    pp.load("AH_merger_pre_factor", merger_pre_factor, 1.);

    pp.load("AH_allow_re_attempt", allow_re_attempt, false);
    pp.load("AH_max_fails_after_lost", max_fails_after_lost, 0);
    pp.load("AH_stop_if_max_fails", stop_if_max_fails, false);

    pp.load("AH_verbose", verbose, (int)AHFinder::MIN);
    pp.load("AH_print_geometry_data", print_geometry_data, false);
    pp.load("AH_re_solve_at_restart", re_solve_at_restart, false);

    // sanity checks
    CH_assert(solve_interval > 0 && print_interval > 0);
    CH_assert(level_to_run >= 0 && level_to_run <= a_p.max_level);

    // if box division ends up with size less than 3, PETSc will
    // complain (this only gives an estimate of the box side)
    int size = 1;
#if CH_MPI
    MPI_Comm_size(Chombo_MPI::comm, &size);
#endif
    size = std::min(num_ranks, size);
#if CH_SPACEDIM == 3
    CH_assert(num_points_u > 0 && num_points_v > 0);
    CH_assert(num_points_u / sqrt(size) >= 3); // make sure for size 'u'
    CH_assert(num_points_v / sqrt(size) >= 3); // make sure for size 'v'
#elif CH_SPACEDIM == 2
    CH_assert(num_points_u > 0);
    CH_assert(num_points_u / size >= 3); // make sure for size 'u'
#endif

    // load vars to write to coord files
    num_extra_vars = 0;
    extra_contain_diagnostic = 0;

    int AH_num_write_vars;
    pp.load("AH_num_write_vars", AH_num_write_vars, 0);
    if (AH_num_write_vars > 0)
    {
        std::vector<std::string> AH_write_var_names(AH_num_write_vars, "");
        pp.load("AH_write_vars", AH_write_var_names, AH_num_write_vars);
        for (const std::string &full_name : AH_write_var_names)
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

            // first assume write_var is a normal evolution var
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
}

/////////////////////////////////////////////////////////
// PETSc control methods
/////////////////////////////////////////////////////////

void AHFinder::set_num_ranks(int a_num_ranks)
{
    if (m_initialized)
        return;

#ifdef CH_MPI
    if (m_mpi_group != MPI_GROUP_NULL)
        MPI_Group_free(&m_mpi_group);

    if (m_mpi_comm != MPI_COMM_NULL)
        MPI_Comm_free(&m_mpi_comm);

    if (a_num_ranks > 0) // otherwise use the whole Chombo communicator
    {
        // don't make more ranks than there are processes
        int size;
        MPI_Comm_size(Chombo_MPI::comm, &size);
        a_num_ranks = std::min(a_num_ranks, size);

        int rank;
        MPI_Comm_rank(Chombo_MPI::comm, &rank);
        if (rank == 0)
            std::cout << "Using PETSc with " << a_num_ranks << " ranks"
                      << std::endl;

        MPI_Group MPI_GROUP_WORLD;
        MPI_Comm_group(Chombo_MPI::comm, &MPI_GROUP_WORLD);

        // (TF): for things like 'SmallDataIO', rank 0 of Chombo_MPI::comm
        // is the one that prints (not rank 0 from PETSC_COMM_WORLD) so pay
        // attention before taking rank 0 of Chombo out of PETSc. As of
        // 26/05/2019, I tested with 'int range[3] = { 1, a_num_ranks , 1 }'
        // and all the AH stuff was compatible with doing so
        // Not true anymore (19/08/2020) due to how
        // ApparentHorizon::write_coords_file uses rank 0 of PETSc to transfer
        // data and to write, but SmallDataIO uses rank 0 of Chombo to write.
        int range[3] = {0, a_num_ranks - 1, 1};
        MPI_Group_range_incl(MPI_GROUP_WORLD, 1,
                             reinterpret_cast<int(*)[3]>(&range), &m_mpi_group);
        MPI_Group_free(&MPI_GROUP_WORLD);

        MPI_Comm_create(Chombo_MPI::comm, m_mpi_group, &m_mpi_comm);
    }
#endif
}

//! initialize PETSc and its MPI sub-communicator
PetscErrorCode AHFinder::PETSc_initialize(int a_num_ranks)
{
    set_num_ranks(a_num_ranks);

    PetscErrorCode err = 0;
    if (!m_initialized)
    {
#ifdef CH_MPI
        if (m_mpi_comm != MPI_COMM_NULL)
            PETSC_COMM_WORLD = m_mpi_comm;
        else // use Chombo communicator if no 'set_num_ranks' was called (or
             // if called with 'a_num_ranks'<=0)
            PETSC_COMM_WORLD = Chombo_MPI::comm;
#endif

        if (AHFinder::is_rank_active())
            err = PetscInitializeNoArguments();

        if (!err)
            m_initialized = true;
    }
    return err;
}

//! finalize PETSc
PetscErrorCode AHFinder::PETSc_finalize()
{
    PetscErrorCode err = 0;
    if (m_initialized)
    {
        if (AHFinder::is_rank_active())
            err = PetscFinalize();

        if (!err)
            m_initialized = false;
    }
    return err;
}

//! true if part of PETSc MPI sub-communicator
bool AHFinder::is_rank_active()
{
#ifdef CH_MPI
    if (m_mpi_group == MPI_GROUP_NULL)
        return true;
    else
    {
        int rank;
        MPI_Group_rank(m_mpi_group, &rank);
        return rank != MPI_UNDEFINED;
    }
#else
    return true;
#endif
}

/////////////////////////////////////////////////////////
// initialize some "static" variables of AHFinder class
/////////////////////////////////////////////////////////

bool AHFinder::m_initialized = false;

#ifdef CH_MPI
MPI_Group AHFinder::m_mpi_group = MPI_GROUP_NULL;
MPI_Comm AHFinder::m_mpi_comm = MPI_COMM_NULL;
#endif

#endif
