/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHFINDER_HPP_
#error "This file should only be included through AHFinder.hpp"
#endif

#ifndef _AHFINDER_IMPL_HPP_
#define _AHFINDER_IMPL_HPP_

#include "FilesystemTools.hpp"

template <class SurfaceGeometry, class AHFunction>
AHFinder<SurfaceGeometry, AHFunction>::~AHFinder()
{
    // destroy horizon pointers and finalize PETSc
    for (auto &ah : m_apparent_horizons)
        delete ah;
    PETScCommunicator::finalize();
}

template <class SurfaceGeometry, class AHFunction>
int AHFinder<SurfaceGeometry, AHFunction>::add_ah(
    const SurfaceGeometry &a_coord_system, double a_initial_guess,
    const AHParams &a_params, bool solve_first_step)
{
    return add_ah(a_coord_system, a_initial_guess, a_params,
                  typename AHFunction::params(), solve_first_step);
}

template <class SurfaceGeometry, class AHFunction>
int AHFinder<SurfaceGeometry, AHFunction>::add_ah(
    const SurfaceGeometry &a_coord_system, double a_initial_guess,
    const AHParams &a_params, const typename AHFunction::params &a_func_params,
    bool solve_first_step)
{
    PETScCommunicator::initialize(a_params.num_ranks);

    if (!FilesystemTools::directory_exists(a_params.stats_path))
        FilesystemTools::mkdir_recursive(a_params.stats_path);

    if (!FilesystemTools::directory_exists(a_params.coords_path))
        FilesystemTools::mkdir_recursive(a_params.coords_path);

    // determine how many AH there are already
    const int num_ah = m_apparent_horizons.size();

    AHInterpolation interp(a_coord_system, m_interpolator);

    m_apparent_horizons.push_back(
        new ApparentHorizon<SurfaceGeometry, AHFunction>(
            interp, a_initial_guess, a_params, a_func_params,
            a_params.stats_prefix + std::to_string(num_ah + 1),
            a_params.coords_prefix + std::to_string(num_ah + 1) + "_",
            solve_first_step));
    m_merger_pairs.push_back({-1, -1});

    return num_ah;
}

template <class SurfaceGeometry, class AHFunction>
int AHFinder<SurfaceGeometry, AHFunction>::add_ah_merger(
    int ah1, int ah2, const AHParams &a_params)
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
    SurfaceGeometry coord_system =
        m_apparent_horizons[ah1]->get_ah_interp().get_coord_system();
    coord_system.set_origin(origin_merger);

    auto &function_to_optimize_params =
        m_apparent_horizons[ah1]->get_petsc_solver().m_func_params;

    int num = add_ah(coord_system, initial_guess_merger, a_params,
                     function_to_optimize_params, do_solve);
    m_merger_pairs[num] = {ah1, ah2};

    return num;
}

template <class SurfaceGeometry, class AHFunction>
void AHFinder<SurfaceGeometry, AHFunction>::solve(double a_dt, double a_time,
                                                  double a_restart_time)
{
    CH_TIME("AHFinder<SurfaceGeometry, AHFunction>::solve");

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

                if (ah_solved[pair.first] == NOT_SOLVED ||
                    ah_solved[pair.second] == NOT_SOLVED)
                    continue; // continue if one of the 'parents' was not
                              // dealt with yet

                if (!(m_apparent_horizons[pair.first]->has_been_found() &&
                      m_apparent_horizons[pair.second]->has_been_found()))
                {
                    if (m_apparent_horizons[i]->m_params.verbose >
                        AHParams::NONE)
                    {
                        pout() << "Skipping BH #" << i << std::endl;
                    }
                    ++ahs_solved;
                    ah_solved[i] = SKIPPED;
                    continue; // skip if one of the parents was never found yet
                }

                if (ah_solved[pair.first] == TOO_FAR ||
                    ah_solved[pair.second] == TOO_FAR)
                { // if one of the 'parents' is a merger too far away to be
                  // solved, skip children too
                    if (m_apparent_horizons[i]->m_params.verbose >
                        AHParams::NONE)
                    {
                        pout() << "Skipping BH #" << i << std::endl;
                    }
                    ++ahs_solved;
                    ah_solved[i] = TOO_FAR;
                    continue;
                }

                double initial_guess_merger;
                std::array<double, CH_SPACEDIM> center_merger;
                bool do_solve =
                    solve_merger(pair.first, pair.second, initial_guess_merger,
                                 center_merger);

                // if distance is too large, ignore this one and move on
                if (!do_solve)
                {
                    // already printed 'Skipping' when calling 'solve_merger'
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
                if (m_apparent_horizons[i]->m_params.verbose > AHParams::NONE)
                    pout() << "Solving AH #" << i << std::endl;

                m_apparent_horizons[i]->solve(a_dt, a_time, a_restart_time);

                ++ahs_solved;
                ah_solved[i] = SOLVED;
            }
            else
            {
                if (m_apparent_horizons[i]->m_params.verbose > AHParams::NONE)
                {
                    pout() << "Skipping BH #" << i << ". Not good to go."
                           << std::endl;
                }
                ++ahs_solved;
                ah_solved[i] = SKIPPED;
            }
        }
    }

    delete ah_solved;
}

template <class SurfaceGeometry, class AHFunction>
bool AHFinder<SurfaceGeometry, AHFunction>::need_diagnostics(
    double a_dt, double a_time) const
{
    bool out = false;
    for (auto &ah : m_apparent_horizons)
        out |=
            ah->m_params.extra_contain_diagnostic && ah->do_print(a_dt, a_time);
    return out;
}

template <class SurfaceGeometry, class AHFunction>
bool AHFinder<SurfaceGeometry, AHFunction>::solve_merger(
    int ah1, int ah2, double &initial_guess_merger,
    std::array<double, CH_SPACEDIM> &center_merger)
{
    // SKIP if 'parents' not yet close enough
    auto AH1 = m_apparent_horizons[ah1];
    auto AH2 = m_apparent_horizons[ah2];
    auto center1 = AH1->get_center();
    auto center2 = AH2->get_center();

    auto initial_guess1 = AH1->get_petsc_solver().get_initial_guess();
    auto initial_guess2 = AH2->get_petsc_solver().get_initial_guess();

    if (m_merger_pairs[ah1].first >= 0)
    {
        double merger_pre_factor =
            std::max(m_apparent_horizons[m_merger_pairs[ah1].first]
                         ->m_params.merger_pre_factor,
                     m_apparent_horizons[m_merger_pairs[ah1].second]
                         ->m_params.merger_pre_factor);
        if (merger_pre_factor >
            0.) // undo the factor is this a merger from a merger
            initial_guess1 /= (merger_pre_factor * 4.);
    }
    if (m_merger_pairs[ah2].first >= 0)
    {
        double merger_pre_factor =
            std::max(m_apparent_horizons[m_merger_pairs[ah2].first]
                         ->m_params.merger_pre_factor,
                     m_apparent_horizons[m_merger_pairs[ah2].second]
                         ->m_params.merger_pre_factor);
        if (merger_pre_factor >
            0.) // undo the factor is this a merger from a merger
            initial_guess2 /= (merger_pre_factor * 4.);
    }

    double initial_guess_sum = (initial_guess1 + initial_guess2);

    // some tests revealed it was about 1.2 * initial_guess_sum, but 1.5 to
    // ensure we don't catch the inner AH
    double merger_pre_factor = std::max(AH1->m_params.merger_pre_factor,
                                        AH2->m_params.merger_pre_factor);
    initial_guess_merger = merger_pre_factor * 4. * initial_guess_sum;

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

    double min_distance = merger_search_factor * 4. * initial_guess_sum;

    bool do_solve = false;

    if (AH1->has_been_found() && AH2->has_been_found())
    {
        if (AH1->get_converged() && AH2->get_converged())
        {
            do_solve = merger_search_factor <= 0. || distance <= min_distance;

            if (do_solve)
                // make sure twice the radius of the guess is bigger than AH
                // distance
                CH_assert(min_distance <= 2. * initial_guess_merger);

            if (AH1->m_params.verbose > AHParams::NONE)
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
        else
        {
            // if AH1 or AH2 was lost, probably merger is already around
            do_solve = true;
        }
    }
    else
    {
        if (AH1->m_params.verbose > AHParams::NONE)
        {
            pout() << "Original BH not yet found. Skipping solve for merger."
                   << std::endl;
        }
    }

    return do_solve;
}

template <class SurfaceGeometry, class AHFunction>
void AHFinder<SurfaceGeometry, AHFunction>::set_origins(
    const std::vector<std::array<double, CH_SPACEDIM>> &origins,
    bool includes_mergers)
{
    const int num_ah = m_apparent_horizons.size();

    int non_merger_counter = 0; // for 'includes_mergers == false'
    for (int i = 0; i < num_ah; ++i)
    {
        if (includes_mergers)
        {
            m_apparent_horizons[i]->set_origin(origins[i]);
        }
        else
        {
            // if it is a merger that has not yet been found
            if (m_merger_pairs[i].first >= 0)
            {
                if (!m_apparent_horizons[i]->has_been_found())
                {
                    auto origin1 = m_apparent_horizons[m_merger_pairs[i].first]
                                       ->get_origin();
                    auto origin2 = m_apparent_horizons[m_merger_pairs[i].second]
                                       ->get_origin();
                    std::array<double, CH_SPACEDIM> origin;
                    FOR1(a) { origin[a] = (origin1[a] + origin2[a]) / 2.; }
                    m_apparent_horizons[i]->set_origin(origin);
                }
            }
            else
            {
                m_apparent_horizons[i]->set_origin(origins[non_merger_counter]);
                ++non_merger_counter;
            }
        }
    }
}

#endif // _AHFINDER_IMPL_HPP_
