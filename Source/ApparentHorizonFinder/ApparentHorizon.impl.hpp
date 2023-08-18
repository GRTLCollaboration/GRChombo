/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _APPARENTHORIZON_HPP_
#error "This file should only be included through ApparentHorizon.hpp"
#endif

#ifndef _APPARENTHORIZON_IMPL_HPP_
#define _APPARENTHORIZON_IMPL_HPP_

#include "FilesystemTools.hpp"
#include "GRAMR.hpp"
#include "PETScCommunicator.hpp"
#include "SmallDataIO.hpp"
#include "SmallDataIOReader.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"

// for 1D interpolation of restart centers
#include "Lagrange.hpp"
#include "SimpleArrayBox.hpp"
#include "SimpleInterpSource.hpp"

// Chombo MPI functions
#include "SPMD.H"

template <class SurfaceGeometry, class AHFunction>
ApparentHorizon<SurfaceGeometry, AHFunction>::ApparentHorizon(
    const AHInterpolation &a_interp, const AHInitialGuessPtr &a_initial_guess,
    const AHParams &a_params, const std::string &a_stats,
    const std::string &a_coords, bool solve_first_step)
    : m_params(a_params),

      m_stats(a_stats), m_coords(a_coords),

      m_printed_once(false),

      m_printed_after_restart(false),

      m_converged(0),

      m_has_been_found(false),

      m_num_failed_convergences(0),

      m_old_centers(3, a_interp.get_coord_system().get_origin()),

      m_max_F(0.), m_min_F(0.), m_ave_F(0.), m_std_F(0.),

      m_integration_methods({
          a_interp.get_coord_system().get_recommended_integration_method_u(
              a_params.num_points_u)
#if CH_SPACEDIM == 3
              ,
              a_interp.get_coord_system().get_recommended_integration_method_v(
                  a_params.num_points_v)
#endif
      }),

      m_area(NAN),
#if GR_SPACEDIM == 3
      m_mass(NAN), m_linear_momentum_P_norm(NAN),
#if CH_SPACEDIM == 3
      m_irreducible_mass(NAN), m_spin(NAN), m_spin_z_alt(NAN),
      m_dimensionless_spin_vector({NAN}), m_linear_momentum_P({NAN}),
#endif
#endif

      origin_already_updated(false),

      solver(a_interp, a_initial_guess, a_params)
{
    set_origin(a_interp.get_coord_system().get_origin());
    check_integration_methods();
    restart(solve_first_step);
}

template <class SurfaceGeometry, class AHFunction>
bool ApparentHorizon<SurfaceGeometry, AHFunction>::good_to_go(
    double a_dt, double a_time) const
{
    if (a_time < m_params.start_time - 1.e-7 /*just to be safe*/)
        return false;
    if (m_params.give_up_time >= 0. && a_time >= m_params.give_up_time &&
        !get_converged())
        return false;

    bool is_lost =
        (m_params.max_fails_after_lost >= 0
             ? m_num_failed_convergences > m_params.max_fails_after_lost
             : false);

    // stop if it has been found but was lost
    return do_solve(a_dt, a_time) && !(has_been_found() && is_lost);
}

template <class SurfaceGeometry, class AHFunction>
bool ApparentHorizon<SurfaceGeometry, AHFunction>::do_solve(double a_dt,
                                                            double a_time) const
{
    CH_assert(a_dt != 0); // Check if time was set!
    return !(((int)(std::round(a_time / a_dt))) % m_params.solve_interval);
}
template <class SurfaceGeometry, class AHFunction>
bool ApparentHorizon<SurfaceGeometry, AHFunction>::do_print(double a_dt,
                                                            double a_time) const
{
    CH_assert(a_dt != 0); // Check if time was set!
    return (get_converged() || !has_been_found()) &&
           !(((int)(std::round(a_time / a_dt))) %
             (m_params.solve_interval * m_params.print_interval));
}

template <class SurfaceGeometry, class AHFunction>
const std::array<double, CH_SPACEDIM> &
ApparentHorizon<SurfaceGeometry, AHFunction>::get_origin() const
{
    return solver.get_origin();
}

template <class SurfaceGeometry, class AHFunction>
const AHInterpolation_t<SurfaceGeometry, AHFunction> &
ApparentHorizon<SurfaceGeometry, AHFunction>::get_ah_interp() const
{
    return solver.m_interp;
}

template <class SurfaceGeometry, class AHFunction>
PETScAHSolver<SurfaceGeometry, AHFunction> &
ApparentHorizon<SurfaceGeometry, AHFunction>::get_petsc_solver()
{
    return solver;
}

template <class SurfaceGeometry, class AHFunction>
void ApparentHorizon<SurfaceGeometry, AHFunction>::set_origin(
    const std::array<double, CH_SPACEDIM> &a_origin)
{
    solver.set_origin(a_origin);

    origin_already_updated = true;
}

template <class SurfaceGeometry, class AHFunction>
const std::array<double, CH_SPACEDIM> &
ApparentHorizon<SurfaceGeometry, AHFunction>::get_center() const
{
    return m_old_centers[0];
}
template <class SurfaceGeometry, class AHFunction>
bool ApparentHorizon<SurfaceGeometry, AHFunction>::get_converged() const
{
    return m_converged;
}

template <class SurfaceGeometry, class AHFunction>
int ApparentHorizon<SurfaceGeometry, AHFunction>::get_failed_convergences()
    const
{
    return m_num_failed_convergences;
}
template <class SurfaceGeometry, class AHFunction>
bool ApparentHorizon<SurfaceGeometry, AHFunction>::has_been_found() const
{
    return m_has_been_found;
}
template <class SurfaceGeometry, class AHFunction>
double ApparentHorizon<SurfaceGeometry, AHFunction>::get_max_F() const
{
    if (m_max_F <= 0.)
        calculate_minmax_F();
    return m_max_F;
}
template <class SurfaceGeometry, class AHFunction>
double ApparentHorizon<SurfaceGeometry, AHFunction>::get_min_F() const
{
    if (m_min_F <= 0.)
        calculate_minmax_F();
    return m_min_F;
}
template <class SurfaceGeometry, class AHFunction>
double ApparentHorizon<SurfaceGeometry, AHFunction>::get_ave_F() const
{
    if (m_ave_F <= 0.)
        calculate_average_F();
    return m_ave_F;
}
template <class SurfaceGeometry, class AHFunction>
double ApparentHorizon<SurfaceGeometry, AHFunction>::get_std_F() const
{
    if (m_std_F <= 0.)
        calculate_average_F();
    return m_std_F;
}

template <class SurfaceGeometry, class AHFunction>
void ApparentHorizon<SurfaceGeometry, AHFunction>::predict_next_origin()
{
    std::array<double, CH_SPACEDIM> new_center = m_old_centers[0];
    if (m_converged >= 3) // add 2nd derivative
    {
        FOR(a)
        {
            new_center[a] += (m_old_centers[0][a] + m_old_centers[2][a] -
                              2. * m_old_centers[1][a]);
        }
        if (m_params.verbose > AHParams::SOME)
        {
            pout() << "OLD[-2]: (" << m_old_centers[2][0] << ","
                   << m_old_centers[2][1]
#if CH_SPACEDIM == 3
                   << "," << m_old_centers[2][2]
#endif
                   << ")" << std::endl;
        }
    }
    if (m_converged >= 2) // add 1st derivative
    {
        FOR(a) { new_center[a] += (m_old_centers[0][a] - m_old_centers[1][a]); }

        if (m_params.verbose > AHParams::SOME)
        {
            pout() << "OLD[-1]: (" << m_old_centers[1][0] << ","
                   << m_old_centers[1][1]
#if CH_SPACEDIM == 3
                   << "," << m_old_centers[1][2]
#endif
                   << ")" << std::endl;
            pout() << "OLD[0]: (" << m_old_centers[0][0] << ","
                   << m_old_centers[0][1]
#if CH_SPACEDIM == 3
                   << "," << m_old_centers[0][2]
#endif
                   << ")" << std::endl;
        }
        if (m_params.verbose > AHParams::MIN)
        {
            pout() << "Estimated center: (" << new_center[0] << ","
                   << new_center[1]
#if CH_SPACEDIM == 3
                   << "," << new_center[2]
#endif
                   << ")" << std::endl;
        }
    }

    set_origin(new_center);
}

template <class SurfaceGeometry, class AHFunction>
void ApparentHorizon<SurfaceGeometry, AHFunction>::update_old_centers(
    std::array<double, CH_SPACEDIM> new_center)
{
    FOR(a)
    {
        m_old_centers[2][a] = m_old_centers[1][a];
        m_old_centers[1][a] = m_old_centers[0][a];
        m_old_centers[0][a] = new_center[a];
    }
}

template <class SurfaceGeometry, class AHFunction>
void ApparentHorizon<SurfaceGeometry, AHFunction>::solve(double a_dt,
                                                         double a_time,
                                                         double a_restart_time)
{
    CH_TIME("ApparentHorizon::solve");
    if (!good_to_go(a_dt, a_time))
        return;

    solver.m_interp.refresh_interpolator(
        do_print(a_dt, a_time),
        m_params.extra_vars); // (ALL CHOMBO ranks do it!!)

    // estimate the next position where origin will be
    if (!origin_already_updated && get_converged())
    {
        if (m_params.predict_origin)
            predict_next_origin();
        else if (m_params.track_center)
            set_origin(get_center());
    }

    // PETSc processes go inside 'if', others "wait" until 'if' gets to
    // 'solver.m_interp.break_interpolation_loop()'
    if (solver.m_interp.keep_interpolating_if_inactive())
    {
        CH_TIME("ApparentHorizon::solve::solving");

        if (!get_converged())
            solver.reset_initial_guess(); // reset initial guess if diverged (or
                                          // in first attempt)

        solver.solve();
        solver.m_interp.break_interpolation_loop();
    }

    {
        if (m_params.verbose > AHParams::SOME)
        {
            pout() << "In [ApparentHorizon::solve::post-solving]" << std::endl;
        }

        CH_TIME("ApparentHorizon::solve::post-solving");

        bool save_converged = get_converged();

        // ask PETSc if it converged
        check_convergence();

        // retry with initial guess once if failed (and everything was ok
        // before)
        // on the 2nd iteration, 'save_converged' will be false and it will skip
        // the 'if'
        if (save_converged && !get_converged() && m_params.allow_re_attempt)
        {
            --m_num_failed_convergences; // reduce failed and try again
            if (m_params.verbose > AHParams::NONE)
            {
                pout() << "Re-attempting to solve using initial guess."
                       << std::endl;
            }
            solve(a_dt, a_time, a_restart_time);
            return; // do nothing else
        }

        // update center
        // note that after this, the point that corresponds to the (u,v,r)
        // coordinates is 'get_origin()', not 'get_center()'
        if (m_params.track_center)
            calculate_center();

        origin_already_updated = false;

#if GR_SPACEDIM == 3 // GR_SPACEDIM, not CH_SPACEDIM !!!
#if CH_SPACEDIM == 3
        Tensor<1, double> J;
        calculate_ah_quantities(m_area, m_linear_momentum_P, J);
        m_linear_momentum_P_norm =
            sqrt(m_linear_momentum_P[0] * m_linear_momentum_P[0] +
                 m_linear_momentum_P[1] * m_linear_momentum_P[1] +
                 m_linear_momentum_P[2] * m_linear_momentum_P[2]);
        double J_norm = sqrt(J[0] * J[0] + J[1] * J[1] + J[2] * J[2]);
        m_mass = calculate_mass(m_area, J_norm);
        m_irreducible_mass = calculate_irreducible_mass(m_area);
#elif CH_SPACEDIM == 2
        Tensor<1, double> P;
        calculate_ah_quantities(m_area, P);
        // assume that P_y = 0 in Cartoon code (if we calculate it it
        // will not be 0, but that is because in the Cartoon code the
        // AH is only half an horizon, so the P_y gets wrongly
        // calculated)
        m_linear_momentum_P_norm = P[0];
        double J_norm = 0.;
        m_mass = calculate_mass(m_area, J_norm);
#endif

#if CH_SPACEDIM == 3
        m_spin = J_norm / m_mass;
        FOR(a) { m_dimensionless_spin_vector[a] = J[a] / (m_mass * m_mass); }

        m_spin_z_alt = calculate_spin_dimensionless(m_area);
#endif

#else // GR_SPACEDIM != 3
        calculate_ah_quantities(m_area);
#endif

        if (m_params.verbose > AHParams::NONE)
        {
#if GR_SPACEDIM == 3 // GR_SPACEDIM, not CH_SPACEDIM !!!
            pout() << "mass = " << m_mass << std::endl;
#if CH_SPACEDIM == 3
            pout() << "spin = " << m_spin << std::endl;
#endif
            if (m_params.verbose > AHParams::MIN)
            {
#if CH_SPACEDIM == 3
                pout() << "irreducible mass = " << m_irreducible_mass
                       << std::endl;
                pout() << "dimensionless spin vector = ("
                       << m_dimensionless_spin_vector[0] << ", "
                       << m_dimensionless_spin_vector[1] << ", "
                       << m_dimensionless_spin_vector[2] << ")" << std::endl;
                pout() << "dimensionless spin in z (from equator-length "
                          "integral) = "
                       << m_spin_z_alt << std::endl;
#endif
                pout() << "linear momentum norm |P| = "
                       << m_linear_momentum_P_norm << std::endl;
#if CH_SPACEDIM == 3
                pout() << "linear momentum P = (" << m_linear_momentum_P[0]
                       << ", " << m_linear_momentum_P[1] << ", "
                       << m_linear_momentum_P[2] << ")" << std::endl;
#endif
            }
#else // GR_SPACEDIM != 3
            pout() << "area = " << m_area << std::endl;
#endif
        }

        // reset min and max F, to force re-calculation
        m_max_F = 0.;
        m_min_F = 0.;
        m_ave_F = 0.;
        m_std_F = 0.;
    }

    write_outputs(a_dt, a_time, a_restart_time);

    // break if 'AH_stop_if_max_fails == true'
    // only break after re-attempting
    // (if max fails reached or give up time reached)
    if (m_params.stop_if_max_fails &&
        ((m_params.max_fails_after_lost >= 0 &&
          m_num_failed_convergences > m_params.max_fails_after_lost) ||
         (m_params.give_up_time >= 0. && a_time >= m_params.give_up_time &&
          !get_converged())))
        MayDay::Error("Reached max fails. Stopping. Parameter "
                      "'stop_if_max_fails' is set to true.");

    if (m_params.verbose > AHParams::SOME)
    {
        pout() << "ApparentHorizon::solve finished successfully!" << std::endl;
    }
}

template <class SurfaceGeometry, class AHFunction>
void ApparentHorizon<SurfaceGeometry, AHFunction>::write_outputs(
    double a_dt, double a_time, double a_restart_time)
{
    // print step. Printing the step allows the user to modify params as
    // 'solve_interval', 'print_interval' or even 'dt_multiplier' (that change
    // the frequency of writing outputs) and the class will still be able to
    // restart from the correct file  (this only applies for frequency increase,
    // for which the file number will suddenly jump, as for frequency decrease
    // the AH may override old coord files)

    // print stats (area, spin, origin, center) and coordinates
    // stop printing if it stopped converging
    // then in 'write_coords_file' nothing is printed if not converged
    if (do_print(a_dt, a_time))
    {
        CH_TIME("ApparentHorizon::solve::printing");
        if (m_params.verbose > AHParams::MIN)
        {
            pout() << "Printing statistics and coordinates." << std::endl;
        }

        if (!m_printed_after_restart && m_printed_once)
        {
            // this forces 'remove_duplicate_time_data' to work (necessary for
            // example with mergers, when the merger doesn't start right at
            // restart but we still need to do 'remove_duplicate_time_data')
            a_restart_time = a_time;
            m_printed_after_restart = true;
        }

        // write stats
        double fake_dt =
            a_dt * m_params.solve_interval * m_params.print_interval;
        SmallDataIO file(m_params.stats_path + m_stats, fake_dt, a_time,
                         a_restart_time, SmallDataIO::APPEND, !m_printed_once);

        file.remove_duplicate_time_data();

        // std::string coords_filename = file.get_new_file_number(fake_dt,
        // a_time);

        CH_assert(CH_SPACEDIM == 3 || CH_SPACEDIM == 2);
        // first '1' corresponds to 'step'
        // area+mass+irred.mass+spin+spin_vec+spin_alt+|P|+ P vec OR
        // area+mass+|P| in 2D Cartoon OR only area for other cases (e.g.
        // 4D -> 2D cartoon)
        int stats = (GR_SPACEDIM == 3 ? (CH_SPACEDIM == 3 ? 12 : 3) : 1);
        std::vector<double> values(1 + stats +
                                   CH_SPACEDIM * (1 + m_params.track_center));

        auto origin = get_origin();

        int step = std::round(a_time / fake_dt);

        int idx = 0;
        values[idx++] = step;
        values[idx++] = m_area;
#if GR_SPACEDIM == 3 // GR_SPACEDIM, not CH_SPACEDIM !!!
        values[idx++] = m_mass;
#if CH_SPACEDIM == 3
        values[idx++] = m_irreducible_mass;
        values[idx++] = m_spin;
        FOR(a) { values[idx++] = m_dimensionless_spin_vector[a]; }
        values[idx++] = m_spin_z_alt;
#endif
        values[idx++] = m_linear_momentum_P_norm;
#if CH_SPACEDIM == 3
        FOR(a) { values[idx++] = m_linear_momentum_P[a]; }
#endif
#endif

        for (int i = 0; i < CH_SPACEDIM; ++i)
            values[idx++] = origin[i];
        if (m_params.track_center)
            for (int i = 0; i < CH_SPACEDIM; ++i)
                values[idx++] = m_old_centers[0][i];

        // print headers to stats file in the beginning of evolution
        if (!m_printed_once)
        {
            std::vector<std::string> headers(
                1 + stats + CH_SPACEDIM * (1 + m_params.track_center));

            idx = 0;
            headers[idx++] = "file";
            headers[idx++] = "area";
#if GR_SPACEDIM == 3 // GR_SPACEDIM, not CH_SPACEDIM !!!
            headers[idx++] = "mass";
#if CH_SPACEDIM == 3
            headers[idx++] = "irreducible mass";
            headers[idx++] = "spin";
            headers[idx++] = "dimless spin-x";
            headers[idx++] = "dimless spin-y";
            headers[idx++] = "dimless spin-z";
            headers[idx++] = "dimless spin-z-alt";
#endif
            headers[idx++] = "linear mom. |P|";
#if CH_SPACEDIM == 3
            headers[idx++] = "linear mom. Px";
            headers[idx++] = "linear mom. Py";
            headers[idx++] = "linear mom. Pz";
#endif

#endif
            headers[idx++] = "origin_x";
            headers[idx++] = "origin_y";
#if CH_SPACEDIM == 3
            headers[idx++] = "origin_z";
#endif
            if (m_params.track_center)
            {
                headers[idx++] = "center_x";
                headers[idx++] = "center_y";
#if CH_SPACEDIM == 3
                headers[idx++] = "center_z";
#endif
            }

            file.write_header_line(headers);

            m_printed_once = true;
        }

        file.write_time_data_line(values);

        // write coordinates
        solver.m_interp.interpolate_extra_vars(m_params.extra_vars);
        write_coords_file(a_dt, a_time, a_restart_time,
                          m_params.coords_path + m_coords,
                          m_params.print_geometry_data);
    }
}
template <class SurfaceGeometry, class AHFunction>
void ApparentHorizon<SurfaceGeometry, AHFunction>::check_convergence()
{
    CH_TIME("ApparentHorizon::check_convergence");

    // Check SNES convergence reasons in:
    // https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESConvergedReason.html#SNESConvergedReason

    int result;
    if (PETScCommunicator::is_rank_active())
    {
        SNESConvergedReason reason = solver.getConvergedReason();

        // result will be 0 if any of the PETSc ranks says 'reason <=0' (PETSc
        // convergence error)
        int tmp = reason; // convert just to make sure it's 'int'
#ifdef CH_MPI
        MPI_Allreduce(&tmp, &result, 1, MPI_INT, MPI_MIN, Chombo_MPI::comm);
#else
        result = tmp;
#endif
    }
    else
    {
        int ALOT = 100; // smaller than all SNESConvergedReasons
#ifdef CH_MPI
        MPI_Allreduce(&ALOT, &result, 1, MPI_INT, MPI_MIN, Chombo_MPI::comm);
#else
        result = ALOT;
#endif
    }

    bool converged = (result > 0);

    if ((bool)converged)
    {
        ++m_converged;
        m_num_failed_convergences = 0;
    }
    else
    {
        m_converged = 0; // reset
        ++m_num_failed_convergences;
    }

    if (m_params.verbose > AHParams::NONE)
    {
        pout() << (m_converged ? "Solver converged. Horizon FOUND."
                               : "Solver diverged. Horizon NOT found.")
               << std::endl;

        if (m_params.verbose > AHParams::MIN)
            pout() << "SNESConvergedReason = " << result << std::endl;

        if (result == SNES_DIVERGED_MAX_IT)
            pout() << "(maximum iterations reached, try increasing "
                      "AH_SNES_max_iterations)"
                   << std::endl;
    }

    if (!m_has_been_found && m_converged)
        m_has_been_found = true; // finally found :D
}

template <class SurfaceGeometry, class AHFunction>
void ApparentHorizon<SurfaceGeometry, AHFunction>::restart(
    bool solve_first_step)
{
    CH_TIME("ApparentHorizon::restart");

    static double eps = 1.e-7;

    // GUIDING PRINCIPLE : restart should not depend on m_params, as these might
    // have changed! We should be able to figure everything out without them

    const GRAMR &gramr = (static_cast<const GRAMR &>(
        solver.m_interp.get_interpolator()->getAMR()));

    int current_step = AMR::s_step;
    int restart_step = gramr.get_restart_step();

    bool solve_time_0_now = (!m_params.re_solve_at_restart &&
                             current_step == 0 && restart_step < 0);

    // solve if at t=0 and not a restart
    if (solve_time_0_now)
    {
        // 'dt' doesn't matter for 1st time step, so set to 1
        if (solve_first_step)
            solve(1., 0., 0.);
        return;
    }

    double current_time = gramr.getCurrentTime();
    double coarse_dt = (current_step == 0 ? 0. : current_time / current_step);
    double level_dt = coarse_dt * pow(2., -m_params.level_to_run);
    double current_solve_dt = level_dt * m_params.solve_interval;

    // READ STATS

    // get centre from stats file
    std::string file = m_params.stats_path + m_stats + ".dat";
    std::vector<std::vector<double>> stats;
    int idx = 0;
    double old_print_dt = 0.;

    // first check the file exists
    bool stats_file_exists = FilesystemTools::file_exists(file);

    if (!stats_file_exists)
    { // case when the file doesn't exist
        if (m_params.verbose > AHParams::NONE && current_step != 0)
        {
            pout() << "Stats file '" << file
                   << "' does not exist. Assuming AHFinder was not run yet for "
                      "this AH."
                   << std::endl;
        }
        m_printed_once = false;
        m_has_been_found = false;
        m_converged = 0;
        m_num_failed_convergences = 0;
    }
    else
    {
        SmallDataIOReader stats_reader;
        stats_reader.open(file);

        if (!stats_reader.contains_data())
        { // case when it never ran the AH
            if (m_params.verbose > AHParams::NONE && current_step != 0)
            {
                pout() << "Empty stats file '" << file
                       << "'. Assuming AHFinder was not run yet for this AH."
                       << std::endl;
            }
            m_printed_once = false;
            m_has_been_found = false;
            m_converged = 0;
            m_num_failed_convergences = 0;
        }
        else
        {
            stats = stats_reader.get_all_columns();
            // look for stats line with time 'current_time', as there may be
            // further output messed up in the file after the last checkpoint
            // was made
            idx = stats[0].size();
            while ((--idx) >= 0 && stats[0][idx] - current_time > eps)
                ;

            // figure out what was the old 'dt' used
            // this may have changed if 'level_to_run' changed of if
            // 'dt_multiplier' changed, or 'solve_interval', or 'print_interval'
            if (idx > 0)
                old_print_dt = stats[0][idx] - stats[0][idx - 1];
            else if (idx == 0 && stats[0].size() > 1)
                old_print_dt =
                    stats[0][idx + 1] -
                    stats[0][idx]; // there's still time steps forward
            else // if we can't know, just default to the current one
                old_print_dt = current_solve_dt * m_params.print_interval;

            if (m_params.verbose > AHParams::SOME)
            {
                if (idx >= 0)
                    pout() << "stats[0][idx] = " << stats[0][idx] << std::endl;
                pout() << "stats[0].size() = " << stats[0].size() << std::endl;
                pout() << "current_time = " << current_time << std::endl;
                pout() << "old_print_dt = " << old_print_dt << std::endl;
            }

            m_printed_once = true; // don't delete old

            if (idx < 0)
            { // case when job was restarted at a point before the AH was first
              // ever found
                if (m_params.verbose > AHParams::NONE)
                {
                    pout() << "First time step is after restart in '" << file
                           << "'. Assuming AHFinder started after current "
                              "step, but "
                              "before "
                              "PETSc had never converged."
                           << std::endl;
                }
                m_has_been_found = false;
                m_converged = 0;
                m_num_failed_convergences = 0;
            }
            else if (stats[0][idx] - current_time <
                     -(old_print_dt == 0. ? eps : old_print_dt - eps))
            { // case when the PETSc stopped converging and stopped printing
              // so there is a mismatch with the last time in the file
                if (m_params.verbose > AHParams::NONE)
                {
                    pout() << "Last time step not found in '" << file
                           << "'. Assuming AH stopped converging." << std::endl;
                    pout() << "(number of failed convergences will be set to 1)"
                           << std::endl;
                }
                m_has_been_found = true;
                m_converged = 0;
                m_num_failed_convergences = 1;
            }
            else if (std::isnan(stats[2][idx]))
            { // case when the simulation was still going without ever
              // having converged
              // OR when job was restarted at a point before the AH was first
              // ever found
                if (m_params.verbose > AHParams::NONE)
                {
                    pout() << "Last time step is NAN in '" << file
                           << "'. Assuming AH wasn't found and PETSc didn't "
                              "converged."
                           << std::endl;
                    pout() << "(number of failed convergences will be set to 1)"
                           << std::endl;
                }
                m_has_been_found = false;
                m_converged = 0;
                m_num_failed_convergences = 1;
            }
            else
            { // case when the AH was found and there is a coordinate file
                if (m_params.verbose > AHParams::NONE)
                {
                    pout() << "Last time step found in '" << file
                           << "'. Reading coordinates file." << std::endl;
                }
                m_has_been_found = true;
                m_converged = 1;
                m_num_failed_convergences = 0;
            }
        }
    }
    // force restart if interpolation was needed (because #points changed)
    // or if time step found is not the current time
    // 'int' in order to use MPI_INT later
    int force_restart = false;

    // if not converged, leave
    // origin as the initial guessed origin
    // and 'solve' will reset coords to initial guess
    if (m_converged)
    {
        ////////////////////
        // READ ORIGIN
        ////////////////////

        int cols = stats.size();

        std::array<double, CH_SPACEDIM> origin;
#if GR_SPACEDIM == 3 // GR_SPACEDIM, not CH_SPACEDIM !!!
#if CH_SPACEDIM == 3
        bool was_center_tracked =
            (cols >
             14 +
                 CH_SPACEDIM); // 14 for time + file + area + mass + irreducible
                               // + spin + spin_vector + spin_z_alt + |P| + P
                               // vec
#else
        // 3D -> 2D Cartoon method
        bool was_center_tracked =
            (cols > 5 + CH_SPACEDIM); // 5 for time + file + area + mass + |P|
#endif
#else
        bool was_center_tracked =
            (cols > 3 + CH_SPACEDIM); // 3 for time + file + area
#endif

        int offset = CH_SPACEDIM + was_center_tracked * CH_SPACEDIM;
        FOR(a) { origin[a] = stats[cols - offset + a][idx]; }

        if (m_params.verbose > AHParams::NONE)
        {
            pout() << "Setting origin from stats file '"
                   << m_params.stats_path + m_stats << "' at (" << origin[0]
                   << "," << origin[1]
#if CH_SPACEDIM == 3
                   << "," << origin[2]
#endif
                   << ")" << std::endl;
        }

        set_origin(origin);

        ////////////////////
        // READ CENTERS
        ////////////////////

        // this is only relevant if track_center was on
        if (was_center_tracked)
        {
            // SET OLD CENTERS directly if old 'dt' matches current 'dt'
            if (std::abs(old_print_dt - current_solve_dt) < eps)
            {
                if (m_params.verbose > AHParams::NONE &&
                    m_params.predict_origin)
                {
                    pout() << "Reading old centers to predict next origin."
                           << std::endl;
                }

                // skip i=0, this one is directly in the file and is done after
                for (int i = m_old_centers.size() - 1; i > 0; --i)
                {
                    if (idx >= i && !std::isnan(stats[2][idx - 1]))
                    {
                        m_converged++;
                        update_old_centers({
#if CH_SPACEDIM == 3
                            stats[cols - 3][idx - i],
#endif
                                stats[cols - 2][idx - i],
                                stats[cols - 1][idx - i]
                        });
                    }
                }
            }
            else
            {
                // calculate last NAN to make sure we don't calculate points of
                // when AH didn't converge
                int last_nan_idx = stats[0].size();
                while ((--last_nan_idx) >= 0 &&
                       !std::isnan(stats[2][last_nan_idx]))
                    ;

                // interpolate old centers when dt has changed
                std::vector<double> old_centers_time_index;
                // skip i=0, this one is directly in the file and is done after
                for (int i = 1; i < m_old_centers.size(); ++i)
                {
                    double index = idx - (current_solve_dt * i) / old_print_dt;

                    if (index >= last_nan_idx)
                        old_centers_time_index.push_back(index);
                }

                if (m_params.verbose > AHParams::NONE &&
                    m_params.predict_origin)
                {
                    if (old_centers_time_index.size() == 0)
                        pout() << "AH time step changed and no old centers "
                                  "able to "
                                  "use for prediction."
                               << std::endl;
                    else
                        pout()
                            << "Old AH time step, " << old_print_dt
                            << ", is different than current one, "
                            << current_solve_dt
                            << " (this might happen if 'AH_print_interval != "
                               "1').\nRecovering "
                            << old_centers_time_index.size()
                            << " old centers from interpolation, for accurate "
                               "prediction of next origin."
                            << std::endl;
                }

                int rows = stats[0].size();
                SimpleInterpSource<1> source({rows}, {old_print_dt});
                SimpleArrayBox<1> box_x({rows}, stats[cols - CH_SPACEDIM]);
                SimpleArrayBox<1> box_y({rows}, stats[cols - CH_SPACEDIM + 1]);
#if CH_SPACEDIM == 3
                SimpleArrayBox<1> box_z({rows}, stats[cols - CH_SPACEDIM + 2]);
#endif

                const int Order = 2;
                Lagrange<Order, 1> interpolator(source);
                for (int i = old_centers_time_index.size() - 1; i >= 0; --i)
                {
                    std::array<double, CH_SPACEDIM> old_center;

                    interpolator.setup({0}, {old_centers_time_index[i]});
                    old_center[0] = interpolator.interpData(box_x);
                    old_center[1] = interpolator.interpData(box_y);
#if CH_SPACEDIM == 3
                    old_center[2] = interpolator.interpData(box_z);
#endif

                    m_converged++;
                    update_old_centers(old_center);

                    if (m_params.verbose > AHParams::SOME)
                    {
                        pout() << "OLD[-" << i + 1 << "] = (" << old_center[0]
                               << "," << old_center[1]
#if CH_SPACEDIM == 3
                               << "," << old_center[2]
#endif
                               << ")" << std::endl;
                    }
                }
            }
        }

        // save current center no matter what
        // place back OLD[0] - THE CENTER
        update_old_centers({
#if CH_SPACEDIM == 3
            stats[cols - 3][idx],
#endif
                stats[cols - 2][idx], stats[cols - 1][idx]
        });

        ////////////////////
        // READ COORDINATES
        ////////////////////

        // determine last step in which there was an AH output
        int coords_file_number = stats[1][idx];
        std::string coords_filename =
            SmallDataIO::get_new_filename(m_params.coords_path + m_coords,
                                          1. /*fake dt*/, coords_file_number);

        SmallDataIOReader coords_reader;
        coords_reader.open(coords_filename);
        auto coords = coords_reader.get_all_columns();

        if (m_params.verbose > AHParams::NONE)
        {
            pout() << "Setting Initial Guess to previous file '"
                   << coords_filename << "' found." << std::endl;
        }

        // now write local points based on file (or interpolated file)
        if (PETScCommunicator::is_rank_active())
        {
            // doesn't change anything if number of points remained the same
            force_restart |= solver.interpolate_ah(coords);
        }

#ifdef CH_MPI
        int tmp = force_restart;
        MPI_Allreduce(&tmp, &force_restart, 1, MPI_INT, MPI_LOR,
                      Chombo_MPI::comm);
#endif
    }

    // if t=0, solve only if it's a restart and if solve_first_step = true (it
    // may be false for example for a merger of a binary, for which we don't
    // want to solve!)
    bool solve_no_matter_what =
        (m_params.re_solve_at_restart &&
         (current_step != 0 ||
          (current_step == 0 && solve_first_step && restart_step >= 0)));

    if (force_restart || solve_no_matter_what)
        solve(current_time == 0. ? 1. : level_dt, current_time, current_time);
}

template <class SurfaceGeometry, class AHFunction>
void ApparentHorizon<SurfaceGeometry, AHFunction>::write_coords_file(
    double a_dt, double a_time, double a_restart_time,
    const std::string &filename, bool write_geometry_data) const
{
    CH_TIME("ApparentHorizon::write_coords_file");

    if (!get_converged())
    {
        // old way: writing an empty file
        // new way: just don't write any coords file if no AH was founds
        // coords_file.write_header_line({"Horizon NOT found."}, "");
        return;
    }

    if (m_params.verbose > AHParams::NONE && write_geometry_data)
        pout() << "Writing geometry data." << std::endl;

    CH_assert(a_dt != 0); // Check if time was set!!

    double fake_dt = a_dt * m_params.solve_interval * m_params.print_interval;
    SmallDataIO file(filename, fake_dt, a_time, a_restart_time,
                     SmallDataIO::NEW, true);

    unsigned num_components; // only related to 'write_geometry_data'
    if (write_geometry_data)
        num_components = AHFunction::num_write_vars();
    else
        num_components = 0;

    // Write headers
    std::vector<std::string> components(num_components +
                                        m_params.num_extra_vars);
    int el = 0;

    // extra vars headers
    for (auto &var : m_params.extra_vars)
    {
        int var_enum = std::get<0>(var.second);
        VariableType var_type = std::get<1>(var.second);
        int der_type = std::get<2>(var.second);

        std::string var_name;
        if (var_type == VariableType::evolution)
            var_name = UserVariables::variable_names[var_enum];
        else
            var_name = DiagnosticVariables::variable_names[var_enum];

        static std::array<std::string, CH_SPACEDIM> xyz{D_DECL("x", "y", "z")};

        if (der_type == 0)
            components[el++] = var_name;
        else if (der_type == 1)
            for (int i = 0; i < CH_SPACEDIM; ++i)
                components[el++] = "d" + xyz[i] + "_" + var_name;
        else // if (der_type == 2)
            for (int i = 0; i < CH_SPACEDIM; ++i)
                for (int j = i; j < CH_SPACEDIM; ++j)
                    components[el++] =
                        "d" + xyz[i] + "d" + xyz[j] + "_" + var_name;
    }

    // geometry headers
    if (write_geometry_data)
        AHFunction::write_headers(&components[el]);
    file.write_header_line(components, solver.m_interp.get_labels());

    //////////////////
    // method:
    // 1) all PESc ranks send their coordinates to rank 0
    // (non-PETSc processes will do nothing)
    // 2) rank 0 receives and stores all coordinates
    // 3) non-blocking sends finish the communication ('MPI_Wait')
    // 4) rank 0 writes all the information
    // Note: this does NOT assume that rank 0 is part of the PETSc communicator!
    // (even though for now the PETSc processes always include rank 0)
    //////////////////

#ifdef CH_MPI
    // not needed for non-PETSc ranks, but they will be waiting anyway
    MPI_Request *requests = nullptr; // needed for non-blocking communication
#endif
    double *output = nullptr; // needed for non-blocking communication
    unsigned num_components_total =
        num_components + m_params.num_extra_vars + CH_SPACEDIM;

    // step 1
    if (PETScCommunicator::is_rank_active())
    {

#ifdef CH_MPI
        int rank_petsc;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank_petsc);

        // make sure rank 0 of Chombo is rank 0 of PETSc
        // this must be the case just because it is rank 0 of PETSc
        // that will be writing in the file, and SmallDataIO assumes it's rank 0
        // of Chombo (!) the one writing
        // Alternative 1: if rank 0 of PETSc is different than rank 0 of Chombo,
        // transfer all the data to Chombo's rank 0 before writing
        // Alternative 2: find a way to tell SmallDataIO to use another rank for
        // writing
        int rank = 0;
        MPI_Comm_rank(Chombo_MPI::comm, &rank);
        if ((rank != 0 && rank_petsc == 0) || (rank == 0 && rank_petsc != 0))
            MayDay::Error("PETSc's rank 0 should be Chombo's rank 0");
#endif

        int local_total = (solver.m_umax - solver.m_umin);
#if CH_SPACEDIM == 3
        local_total *= (solver.m_vmax - solver.m_vmin);
#endif

#ifdef CH_MPI
        requests = new MPI_Request[local_total];
#endif
        output = new double[local_total * num_components_total];

        int idx = 0;
#if CH_SPACEDIM == 3
        for (int v = solver.m_vmin; v < solver.m_vmax; ++v)
#endif
        {
            for (int u = solver.m_umin; u < solver.m_umax; ++u)
            {
                output[idx * num_components_total] = solver.m_u[idx];
#if CH_SPACEDIM == 3
                output[idx * num_components_total + 1] = solver.m_v[idx];
#endif
                output[idx * num_components_total + SpaceDim - 1] =
                    solver.m_F[idx];

                auto extra = solver.m_interp.get_extra_data(idx);

                int el = 0;

                // extra vars
                for (auto &var : m_params.extra_vars)
                {
                    // int var_enum = std::get<0>(var.second);
                    // VariableType var_type = std::get<1>(var.second);
                    int der_type = std::get<2>(var.second);

                    if (der_type == 0)
                    {
                        output[idx * num_components_total + CH_SPACEDIM +
                               el++] = extra.vars.at(var.first);
                    }
                    else if (der_type == 1)
                    {
                        for (int i = 0; i < CH_SPACEDIM; ++i)
                            output[idx * num_components_total + CH_SPACEDIM +
                                   el++] = extra.d1.at(var.first)[i];
                    }
                    else // if (der_type == 2)
                    {
                        for (int i = 0; i < CH_SPACEDIM * (CH_SPACEDIM + 1) / 2;
                             ++i)
                            output[idx * num_components_total + CH_SPACEDIM +
                                   el++] = extra.d2.at(var.first)[i];
                    }
                }

                // geometry vars
                if (write_geometry_data)
                {
                    const auto data = solver.m_interp.get_data(idx);
                    const auto coords = solver.m_interp.get_coords(idx);
                    const auto coords_cart =
                        solver.m_interp.get_cartesian_coords(idx);
                    AHFunction func(data, coords, coords_cart);
                    func.write_vars(
                        &output[idx * num_components_total + CH_SPACEDIM + el]);
                }

#ifdef CH_MPI
                // all processes send their 'output' to rank 0, who receives
                // and writes everything (to simplify, rank 0 also sends it
                // to itself, otherwise we wouldn't know which tags rank 0
                // has) sends are tagged by global index, so that receives
                // can be unique and writes indexed correctly
#if CH_SPACEDIM == 3
                int idx_global = v * solver.m_num_global_u + u;
#elif CH_SPACEDIM == 2
                int idx_global = u;
#endif
                MPI_Issend(&output[idx * num_components_total],
                           num_components_total, MPI_DOUBLE, 0, idx_global,
                           PETSC_COMM_WORLD, &(requests[idx]));
#else
                // only rank 0 is active -> serial
                file.write_data_line(
                    std::vector<double>(output + idx * num_components_total +
                                            CH_SPACEDIM,
                                        output + idx * num_components_total +
                                            num_components_total),
                    std::vector<double>(output + idx * num_components_total,
                                        output + idx * num_components_total +
                                            CH_SPACEDIM));
#endif

                idx++;
            }
        }

        // step 2
#ifdef CH_MPI
#if CH_SPACEDIM == 3
        int total = solver.m_num_global_u * solver.m_num_global_v;
#elif CH_SPACEDIM == 2
        int total = solver.m_num_global_u;
#endif
        double temp[total * num_components_total];
        if (rank_petsc == 0)
        {
            for (unsigned i = 0; i < total; ++i)
            {
                MPI_Status status;
                MPI_Recv(&temp[i * num_components_total], num_components_total,
                         MPI_DOUBLE, MPI_ANY_SOURCE,
                         i /* use i as tag to make them ordered*/,
                         PETSC_COMM_WORLD, &status);
            }
        }

        // step 3
        idx = 0;
#if CH_SPACEDIM == 3
        for (int v = solver.m_vmin; v < solver.m_vmax; ++v)
#endif
        {
            for (int u = solver.m_umin; u < solver.m_umax; ++u)
            {
                MPI_Wait(&(requests[idx]), MPI_STATUS_IGNORE);
                idx++;
            }
        }

        delete[] requests; // free allocated memory
        delete[] output;

        // step 4
        if (rank_petsc == 0)
        {
            for (unsigned i = 0; i < total; ++i)
            {
                file.write_data_line(
                    std::vector<double>(
                        &temp[i * num_components_total] + CH_SPACEDIM,
                        &temp[i * num_components_total] + num_components_total),
                    std::vector<double>(&temp[i * num_components_total],
                                        &temp[i * num_components_total] +
                                            CH_SPACEDIM));
            }
        }
#endif
    }
}

template <class SurfaceGeometry, class AHFunction>
void ApparentHorizon<SurfaceGeometry, AHFunction>::check_integration_methods()
{
    // check if integration methods are valid given periodicity and number of
    // points
    bool valid_u = m_integration_methods[0].is_valid(
        m_params.num_points_u,
        solver.m_interp.get_coord_system().is_u_periodic());

    const IntegrationMethod &method_default_u =
        solver.m_interp.get_coord_system().get_recommended_integration_method_u(
            solver.m_num_global_u);
    if (!valid_u)
    {
        std::string warn =
            "ApparentHorizon: Simpson IntegrationMethod for u is "
            "not valid with this num_points_u.\n"
            "Reverting to trapezium rule.";
        MayDay::Warning(warn.c_str());
        pout() << warn << std ::endl;
        m_integration_methods[0] = method_default_u;
    }

#if CH_SPACEDIM == 3
    bool valid_v = m_integration_methods[1].is_valid(
        m_params.num_points_v,
        solver.m_interp.get_coord_system().is_v_periodic());
    const IntegrationMethod &method_default_v =
        solver.m_interp.get_coord_system().get_recommended_integration_method_v(
            solver.m_num_global_v);
#endif
}

#if CH_SPACEDIM == 3
// ONLY FOR 3D
// OLD method to calculate the spin using the 1D equator length
// Adopted area integral 'calculate_angular_momentum_J' that allows to get the
// direction of the spin as well
template <class SurfaceGeometry, class AHFunction>
double
ApparentHorizon<SurfaceGeometry, AHFunction>::calculate_spin_dimensionless(
    double a_area)
{
    CH_assert(CH_SPACEDIM == 3);
    CH_TIME("ApparentHorizon::calculate_spin_dimensionless");

    if (!get_converged())
        return NAN;

    double equator_length;
    double integral = 0.; // temporary, but defined outside to be in scope for
                          // non-PETSc processes

    // PETSc processes go inside 'if', others "wait" until 'if' gets to
    // 'solver.m_interp.break_interpolation_loop()'
    if (solver.m_interp.keep_interpolating_if_inactive())
    {
        int idx = 0;

        Vec localF;
        dmda_arr_t in;
        solver.get_dmda_arr_t(localF, in);

        // solver.m_F is already set from solve

        int u_equator = std::round(M_PI / 2.0 / solver.m_du);

        for (int v = solver.m_vmin; v < solver.m_vmax; ++v)
        {
            for (int u = solver.m_umin; u < solver.m_umax; ++u)
            {
                if (u_equator == u)
                {
                    AHDerivData deriv = solver.diff(in, u, v);
                    const auto geometry_data =
                        solver.m_interp.get_geometry_data(idx);
                    const auto data = solver.m_interp.get_data(idx);
                    const auto coords = solver.m_interp.get_coords(idx);
                    const auto coords_cart =
                        solver.m_interp.get_cartesian_coords(idx);
                    AHFunction func(data, coords, coords_cart);
                    auto &g = func.get_metric();

                    double dxdv[3];
                    dxdv[0] = geometry_data.dxdv[0] +
                              deriv.dvF * geometry_data.dxdf[0];
                    dxdv[1] = geometry_data.dxdv[1] +
                              deriv.dvF * geometry_data.dxdf[1];
                    dxdv[2] = geometry_data.dxdv[2] +
                              deriv.dvF * geometry_data.dxdf[2];

                    double norm2 = 0.;
                    FOR(i, j) norm2 += g[i][j] * dxdv[i] * dxdv[j];
                    CH_assert(norm2 >= 0.);

                    double weight = m_integration_methods[1].weight(
                        v, m_params.num_points_v,
                        solver.m_interp.get_coord_system().is_v_periodic());

                    integral += sqrt(norm2) * weight * solver.m_dv;
                }

                idx++;
            }
        }

        solver.restore_dmda_arr_t(localF, in);
        solver.m_interp.break_interpolation_loop();
    }

    // reduction across all Chombo processes (note that 'integral' is 0 for
    // non-PETSc processes) because SmallDataIO uses rank 0 to write this
    // ensures rank 0 will have 'integral' (even though for now the PETSc
    // processes always include rank 0)
#ifdef CH_MPI
    MPI_Allreduce(&integral, &equator_length, 1, MPI_DOUBLE, MPI_SUM,
                  Chombo_MPI::comm);
#else // serial
    equator_length = integral;
#endif

    // mass could be estimated based on "equator_length / (4. * M_PI)"
    // but calculation using area as well is probably more precise

    double factor =
        ((2. * M_PI * a_area / (equator_length * equator_length)) - 1.);
    double spin_dimensionless =
        (factor > 1. ? 0
                     : sqrt(1. - factor * factor)); // factor>1 means numerical
                                                    // error with spin as 0

    return spin_dimensionless;
}
#endif

template <class SurfaceGeometry, class AHFunction>
void ApparentHorizon<SurfaceGeometry, AHFunction>::calculate_ah_quantities(
    double &area
#if GR_SPACEDIM == 3
    ,
    Tensor<1, double> &P
#if CH_SPACEDIM == 3
    ,
    Tensor<1, double> &J
#endif
#endif
)
{
    CH_TIME("ApparentHorizon::calculate_ah_quantities");

    if (!get_converged())
    {
        area = NAN;
#if GR_SPACEDIM == 3
        P = {NAN};
#if CH_SPACEDIM == 3
        J = {NAN};
#endif
#endif
        return;
    }

    // temporary, but defined outside to be
    // in scope for non-PETSc processes
    double integral_area = 0.;
#if GR_SPACEDIM == 3
    Tensor<1, double> integrals_P = {0.};
#if CH_SPACEDIM == 3
    Tensor<1, double> integrals_J = {0.};
#endif
#endif

    // PETSc processes go inside 'if', others "wait" until 'if' gets to
    // 'solver.m_interp.break_interpolation_loop()'
    if (solver.m_interp.keep_interpolating_if_inactive())
    {
        int idx = 0;

        Vec localF;
        dmda_arr_t in;
        solver.get_dmda_arr_t(localF, in);

        // solver.m_F is already set from solve

        const std::array<double, CH_SPACEDIM> &center = get_center();

#if CH_SPACEDIM == 3
        for (int v = solver.m_vmin; v < solver.m_vmax; ++v)
#endif
        {

            double inner_integral_area = 0.;
#if GR_SPACEDIM == 3
            Tensor<1, double> inner_integral_P = {0.};
#if CH_SPACEDIM == 3
            Tensor<1, double> inner_integral_J = {0.};
#endif
#endif

            for (int u = solver.m_umin; u < solver.m_umax; ++u)
            {
                AHDerivData deriv = solver.diff(in, u
#if CH_SPACEDIM == 3
                                                ,
                                                v
#endif
                );
                const auto geometric_data =
                    solver.m_interp.get_geometry_data(idx);
                const auto data = solver.m_interp.get_data(idx);
                const auto coords = solver.m_interp.get_coords(idx);
                const auto coords_cart =
                    solver.m_interp.get_cartesian_coords(idx);
                AHFunction func(data, coords, coords_cart);
                Tensor<2, double> g = func.get_metric();
#if GR_SPACEDIM == 3
                Tensor<2, double> K = func.get_extrinsic_curvature();
                Tensor<1, double> s_L =
                    func.get_level_function_derivative(geometric_data, deriv);
                Tensor<1, double> S_U = func.get_spatial_normal_U(s_L);

                Tensor<1, double> coords_cart_centered;
                FOR(i) { coords_cart_centered[i] = coords_cart[i] - center[i]; }

                // Linear Momentum P
                Tensor<1, double> p_integrand = {0.};
                // double directions[3][3] = {
                //     {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
                // FOR(a, b, c)
                // {
                //     p_integrand[c] += directions[c][a] * S_U[b] * K[a][b];
                // }
                // simplify the commented above to:
                FOR(a, b) { p_integrand[a] += S_U[b] * K[a][b]; }

#if CH_SPACEDIM == 3
                // Angular Momentum J
                Tensor<1, double> spin_integrand = {0.};
                double directions[3][3] = {
                    {0., -coords_cart_centered[2], coords_cart_centered[1]},
                    {coords_cart_centered[2], 0., -coords_cart_centered[0]},
                    {-coords_cart_centered[1], coords_cart_centered[0], 0.}};
                FOR(a, b, c)
                {
                    // as in arXiv:gr-qc/0206008 eq. 25
                    // but computing with 'killing vector' in directions 'x,y,z'
                    // ('x,y,z' killing vector are directions[c][a] =
                    // epsilon[c][j][a] * x[j], just like for ADM Momentum)
                    spin_integrand[c] += directions[c][a] * S_U[b] * K[a][b];
                }
#endif
#endif
                // now calculate sqrt(-g) area element

                // Calculate Jacobian matrix for transformation from Cartesian
                // to (f,u,v) coords
                Tensor<2, double> Jac;
                FOR(k)
                {
                    Jac[0][k] = geometric_data.dxdf[k];
                    Jac[1][k] = geometric_data.dxdu[k];
#if CH_SPACEDIM == 3
                    Jac[2][k] = geometric_data.dxdv[k];
#endif
                }

                // Now do the coordinate transformation
                Tensor<2, double> g_spherical = {0.};
                FOR(i, j, k, l)
                {
                    g_spherical[i][j] += Jac[i][k] * Jac[j][l] * g[k][l];
                }

                // Construct the 2-metric on the horizon in (u,v) coords
                // i.e. substitute df = (df/du)du + (df/dv)dv
                // into the spherical metric
                Tensor<2, double, CH_SPACEDIM - 1> g_horizon = {0.};
                g_horizon[0][0] = g_spherical[1][1] +
                                  g_spherical[0][0] * deriv.duF * deriv.duF +
                                  2.0 * g_spherical[0][1] * deriv.duF;
#if CH_SPACEDIM == 3
                g_horizon[1][1] = g_spherical[2][2] +
                                  g_spherical[0][0] * deriv.dvF * deriv.dvF +
                                  2.0 * g_spherical[0][2] * deriv.dvF;
                g_horizon[0][1] = g_horizon[1][0] =
                    g_spherical[1][2] +
                    g_spherical[0][0] * deriv.duF * deriv.dvF +
                    g_spherical[0][1] * deriv.dvF +
                    g_spherical[0][2] * deriv.duF;
#endif

                double det = TensorAlgebra::compute_determinant(g_horizon);

                double weight = m_integration_methods[0].weight(
                    u, m_params.num_points_u,
                    solver.m_interp.get_coord_system().is_u_periodic());

                double element_area = sqrt(det) * weight * solver.m_du;

// assume a (GR_SPACEDIM - CH_SPACEDIM)-sphere leftover
#if GR_SPACEDIM != CH_SPACEDIM
                double n_sphere = (GR_SPACEDIM - CH_SPACEDIM);
                double element_hd = pow(sqrt(func.get_metric_hd()) *
                                            coords_cart[CH_SPACEDIM - 1],
                                        n_sphere);
                element_area *= element_hd;
#endif

                inner_integral_area += element_area;

#if GR_SPACEDIM == 3
                // Linear Momentum P
                FOR(a)
                {
                    double element_P = p_integrand[a] / (8. * M_PI) *
                                       sqrt(det) * weight * solver.m_du;

// assume a (GR_SPACEDIM - CH_SPACEDIM)-sphere leftover
#if GR_SPACEDIM != CH_SPACEDIM
                    element_P *= element_hd;
#endif

                    // hack for the poles (theta=0,\pi)
                    // where nans appear in 'p_integrand', but
                    // 'det' should be 0
                    if (std::isnan(element_P))
                        element_P = 0.;

                    inner_integral_P[a] += element_P;
                }

#if CH_SPACEDIM == 3
                // Angular Momentum J
                FOR(a)
                {
                    double element_J = spin_integrand[a] / (8. * M_PI) *
                                       sqrt(det) * weight * solver.m_du;

                    // hack for the poles (theta=0,\pi)
                    // where nans appear in 'spin_integrand', but
                    // 'det' should be 0
                    if (std::isnan(element_J))
                        element_J = 0.;

                    inner_integral_J[a] += element_J;
                }
#endif
#endif

                idx++;
            }
#if CH_SPACEDIM == 3
            double weight = m_integration_methods[1].weight(
                v, m_params.num_points_v,
                solver.m_interp.get_coord_system().is_v_periodic());
            double weight_dv = weight * solver.m_dv;
            integral_area += weight_dv * inner_integral_area;
            FOR(a) { integrals_P[a] += weight_dv * inner_integral_P[a]; }
            FOR(a) { integrals_J[a] += weight_dv * inner_integral_J[a]; }
#elif CH_SPACEDIM == 2
            integral_area += inner_integral_area;
#if GR_SPACEDIM == 3
            FOR(a) { integrals_P[a] += inner_integral_P[a]; }
#endif
#endif
        }

// do this separate so that the gamma is only calculated once
#if GR_SPACEDIM != CH_SPACEDIM
        double n_sphere = (GR_SPACEDIM - CH_SPACEDIM);
        // this is 2pi for n=1, 4pi for n=2, 2pi^2 for n=3, ...
        double n_sphere_coeff = 2. * std::pow(M_PI, (n_sphere + 1.) / 2.) /
                                std::tgamma((n_sphere + 1.) / 2.);
        integral_area *= n_sphere_coeff;
#if GR_SPACEDIM == 3
        FOR(a) { integrals_P[a] *= n_sphere_coeff; }
#endif
#endif

        solver.restore_dmda_arr_t(localF, in);
        solver.m_interp.break_interpolation_loop();
    }

    // reduction across all Chombo processes (note that 'integral' is 0 for
    // non-PETSc processes) because SmallDataIO uses rank 0 to write this
    // ensures rank 0 will have 'integral' (even though for now the PETSc
    // processes always include rank 0)
#ifdef CH_MPI
    MPI_Allreduce(&integral_area, &area, 1, MPI_DOUBLE, MPI_SUM,
                  Chombo_MPI::comm);
#if GR_SPACEDIM == 3
    MPI_Allreduce(&integrals_P, &P, GR_SPACEDIM, MPI_DOUBLE, MPI_SUM,
                  Chombo_MPI::comm);
#if CH_SPACEDIM == 3
    MPI_Allreduce(&integrals_J, &J, GR_SPACEDIM, MPI_DOUBLE, MPI_SUM,
                  Chombo_MPI::comm);
#endif
#endif
#else // serial
    area = integral_area;
#if GR_SPACEDIM == 3
    FOR(a) { P[a] = integrals_P[a]; }
#if CH_SPACEDIM == 3
    FOR(a) { J[a] = integrals_J[a]; }
#endif
#endif
#endif
}

template <class SurfaceGeometry, class AHFunction>
std::array<double, CH_SPACEDIM>
ApparentHorizon<SurfaceGeometry, AHFunction>::calculate_center()
{
    CH_TIME("ApparentHorizon::calculate_center");

    if (!get_converged())
    {
#if CH_SPACEDIM == 3
        return {NAN, NAN, NAN};
#elif CH_SPACEDIM == 2
        return {NAN, NAN};
#endif
    }

    // [OLD] Method:
    // Calculate centroid of all the points, by summing the coordinates {x,y,z}
    // of all the points, which are spread across processors, and then reducing
    // them all with MPI. Finally, divide by total number of points summed, to
    // get the centroid
    // Above is not good. Replace by calculating the maximum and minimum point
    // of the surface and estimating the center by their average

    // std::array<double, CH_SPACEDIM> temp = {0., 0., 0.}; // old method
    std::array<double, CH_SPACEDIM> max_temp = get_origin();
    std::array<double, CH_SPACEDIM> min_temp = max_temp;

    int idx = 0;
    if (PETScCommunicator::is_rank_active()) // in principle this wouldn't be
                                             // needed, as non-PETSc processes
                                             // wouldn't even enter the loop
                                             // anyways
    {
#if CH_SPACEDIM == 3
        for (int v = solver.m_vmin; v < solver.m_vmax; ++v)
#endif
        {
            for (int u = solver.m_umin; u < solver.m_umax; ++u)
            {
                for (unsigned i = 0; i < CH_SPACEDIM; ++i)
                {
                    double coord_i =
                        solver.m_interp.get_coord_system().get_grid_coord(
                            i, D_DECL(solver.m_F[idx], solver.m_u[idx],
                                      solver.m_v[idx]));

                    // temp[i] += point[i]; // old method
                    if (coord_i > max_temp[i])
                        max_temp[i] = coord_i;
                    if (coord_i < min_temp[i])
                        min_temp[i] = coord_i;
                }

                idx++;
            }
        }
    }

    std::array<double, CH_SPACEDIM> max_point;
    std::array<double, CH_SPACEDIM> min_point;
    std::array<double, CH_SPACEDIM> center;
    // int idx_sum;

#ifdef CH_MPI
    MPI_Allreduce(&max_temp, &max_point, CH_SPACEDIM, MPI_DOUBLE, MPI_MAX,
                  Chombo_MPI::comm);
    MPI_Allreduce(&min_temp, &min_point, CH_SPACEDIM, MPI_DOUBLE, MPI_MIN,
                  Chombo_MPI::comm);
    // MPI_Allreduce(&center, &center, CH_SPACEDIM, MPI_DOUBLE, MPI_SUM,
    // Chombo_MPI::comm); // old method
    // MPI_Allreduce(&idx, &idx_sum, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
#else // serial
    max_point = max_temp;
    min_point = min_temp;
    // idx_sum = idx;
#endif

    // CH_assert(idx_sum != 0); // old method

    for (unsigned i = 0; i < CH_SPACEDIM; ++i)
    {
        // center[i] /= idx_sum; // old method
        center[i] = (max_point[i] + min_point[i]) / 2.;
    }

    if (m_params.verbose > AHParams::SOME)
    {
        pout() << "max_point: (" << max_point[0] << ", " << max_point[1]
#if CH_SPACEDIM == 3
               << "," << max_point[2]
#endif
               << ")" << std::endl;
        pout() << "min_point: (" << min_point[0] << ", " << min_point[1]
#if CH_SPACEDIM == 3
               << ", " << min_point[2]
#endif
               << ")" << std::endl;
    }

    update_old_centers(center); // set new

    if (m_params.verbose > AHParams::MIN)
    {
        pout() << "center: (" << center[0] << ", " << center[1]
#if CH_SPACEDIM == 3
               << ", " << center[2]
#endif
               << ")" << std::endl;
    }

    return center;
}

template <class SurfaceGeometry, class AHFunction>
void ApparentHorizon<SurfaceGeometry, AHFunction>::calculate_minmax_F() const
{
    if (!get_converged())
        return;

    double local_max = std::numeric_limits<int>::min();
    double local_min = std::numeric_limits<int>::max();
    if (PETScCommunicator::is_rank_active())
    {
        auto local_minmax =
            (std::minmax_element(solver.m_F.begin(), solver.m_F.end()));
        local_min = *(local_minmax.first);
        local_max = *(local_minmax.second);
    }
    double global_max = local_max;
    double global_min = local_min;

#ifdef CH_MPI
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN,
                  Chombo_MPI::comm);
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX,
                  Chombo_MPI::comm);
#endif

    m_max_F = global_max;
    m_min_F = global_min;
}

template <class SurfaceGeometry, class AHFunction>
void ApparentHorizon<SurfaceGeometry, AHFunction>::calculate_average_F() const
{
    if (!get_converged())
        return;

    int local_size = 0;
    double local_sum = 0.;
    double local_sum_sq = 0.;
    if (PETScCommunicator::is_rank_active())
    {
        std::pair<double, double> sums(0., 0.);
        auto lambda_sum = [](std::pair<double, double> sums, double r)
        {
            sums.first += r;      // radius
            sums.second += r * r; // radius^2
            return std::move(sums);
        };

        auto local_sums = std::accumulate(solver.m_F.begin(), solver.m_F.end(),
                                          sums, lambda_sum);

        local_size = solver.m_F.size();
        local_sum = local_sums.first;
        local_sum_sq = local_sums.second;
    }
    int global_size = local_size;
    double global_sum = local_sum;
    double global_sum_sq = local_sum_sq;

#ifdef CH_MPI
    MPI_Allreduce(&local_size, &global_size, 1, MPI_INT, MPI_SUM,
                  Chombo_MPI::comm);
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM,
                  Chombo_MPI::comm);
    MPI_Allreduce(&local_sum_sq, &global_sum_sq, 1, MPI_DOUBLE, MPI_SUM,
                  Chombo_MPI::comm);
#endif

    m_ave_F = global_sum / global_size;
    m_std_F = sqrt(global_sum_sq / global_size -
                   m_ave_F * m_ave_F); // sqrt( variance )
}

#endif // _APPARENTHORIZON_IMPL_HPP_