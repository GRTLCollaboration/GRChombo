/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "BinaryBHLevel.hpp"
#include "AMRReductions.hpp"
#include "BinaryBH.hpp"
#include "BoxLoops.hpp"
#include "CCZ4.hpp"
#include "ChiExtractionTaggingCriterion.hpp"
#include "ChiPunctureExtractionTaggingCriterion.hpp"
#include "ComputePack.hpp"
#include "Constraints.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "PunctureTracker.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"
#include "TraceARemoval.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

// Things to do during the advance step after RK4 steps
void BinaryBHLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck("NaNCheck in specific Advance: "), m_state_new,
                       m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// This initial data uses an approximation for the metric which
// is valid for small boosts
void BinaryBHLevel::initialData()
{
    CH_TIME("BinaryBHLevel::initialData");
    if (m_verbosity)
        pout() << "BinaryBHLevel::initialData " << m_level << endl;

    // Set up the compute class for the BinaryBH initial data
    BinaryBH binary(m_p.bh1_params, m_p.bh2_params, m_dx);

    // setup initial puncture coords for tracking
    // do puncture tracking, just set them once, so on level 0
    if (m_p.track_punctures == 1 && m_level == 0)
    {
        const double coarsest_dt = m_p.coarsest_dx * m_p.dt_multiplier;
        PunctureTracker my_punctures(m_time, m_restart_time, coarsest_dt,
                                     m_p.checkpoint_prefix);
        my_punctures.set_initial_punctures(m_bh_amr,
                                           m_p.initial_puncture_coords);
    }

    // First set everything to zero (to avoid undefinded values in constraints)
    // then calculate initial data
    BoxLoops::loop(make_compute_pack(SetValue(0.), binary), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS);
}

// Things to do after a restart
void BinaryBHLevel::postRestart()
{
    // do puncture tracking, just set them once, so on the top level
    if (m_p.track_punctures == 1 && m_level == m_p.max_level)
    {
        // need to set a temporary interpolator for finding the shift
        // as the happens in setupAMRObject() not amr.run()
        AMRInterpolator<Lagrange<4>> interpolator(m_bh_amr, m_p.origin, m_p.dx,
                                                  m_p.verbosity);
        m_bh_amr.set_interpolator(&interpolator);
        PunctureTracker my_punctures(m_time, m_restart_time, m_dt,
                                     m_p.checkpoint_prefix);
        my_punctures.restart_punctures(m_bh_amr, m_p.initial_puncture_coords);
    }
}

// Calculate RHS during RK4 substeps
void BinaryBHLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                    const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side and set constraints to zero to avoid
    // undefined values
    BoxLoops::loop(CCZ4(m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation),
                   a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

// enforce trace removal during RK4 substeps
void BinaryBHLevel::specificUpdateODE(GRLevelData &a_soln,
                                      const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

// specify the cells to tag
void BinaryBHLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                            const FArrayBox &current_state)
{
    if (m_p.track_punctures == true)
    {
        const vector<double> puncture_masses = {m_p.bh1_params.mass,
                                                m_p.bh2_params.mass};
        std::vector<std::array<double, CH_SPACEDIM>> puncture_coords =
            m_bh_amr.get_puncture_coords();
        BoxLoops::loop(ChiPunctureExtractionTaggingCriterion(
                           m_dx, m_level, m_p.max_level, m_p.extraction_params,
                           puncture_coords, m_p.activate_extraction,
                           m_p.track_punctures, puncture_masses),
                       current_state, tagging_criterion);
    }
    else
    {
        BoxLoops::loop(ChiExtractionTaggingCriterion(
                           m_dx, m_level, m_p.max_level, m_p.extraction_params,
                           m_p.activate_extraction),
                       current_state, tagging_criterion);
    }
}

void BinaryBHLevel::specificPostTimeStep()
{
    CH_TIME("BinaryBHLevel::specificPostTimeStep");

    bool first_step = (m_time == m_dt);

    if (m_p.activate_extraction == 1)
    {
        // Populate the Weyl Scalar values on the grid
        fillAllGhosts();
        BoxLoops::loop(Weyl4(m_p.extraction_params.center, m_dx), m_state_new,
                       m_state_diagnostics, EXCLUDE_GHOST_CELLS);

        // Do the extraction on the min extraction level
        if (m_level == m_p.extraction_params.min_extraction_level())
        {
            CH_TIME("WeylExtraction");
            // Now refresh the interpolator and do the interpolation
            m_gr_amr.m_interpolator->refresh();
            WeylExtraction my_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            my_extraction.execute_query(m_gr_amr.m_interpolator);
        }
    }

    if (m_p.calculate_constraint_norms)
    {
        fillAllGhosts();
        BoxLoops::loop(Constraints(m_dx), m_state_new, m_state_diagnostics,
                       EXCLUDE_GHOST_CELLS);
        if (m_level == 0)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            double L2_Ham = amr_reductions.norm(c_Ham);
            double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3));
            SmallDataIO constraints_file("constraint_norms", m_dt, m_time,
                                         m_restart_time, SmallDataIO::APPEND,
                                         first_step);
            constraints_file.remove_duplicate_time_data();
            if (first_step)
            {
                constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
            }
            constraints_file.write_time_data_line({L2_Ham, L2_Mom});
        }
    }

    // do puncture tracking on requested level
    if (m_p.track_punctures == 1 && m_level == m_p.puncture_tracking_level)
    {
        CH_TIME("PunctureTracking");
        // only do the write out for every coarsest level timestep
        bool write_punctures = false;
        const double coarsest_dt = m_p.coarsest_dx * m_p.dt_multiplier;
        const double remainder = fmod(m_time, coarsest_dt);
        PunctureTracker my_punctures(m_time, m_restart_time, m_dt,
                                     m_p.checkpoint_prefix);
        if (min(abs(remainder), abs(remainder - coarsest_dt)) < 1.0e-8)
        {
            write_punctures = true;
        }
        my_punctures.execute_tracking(m_bh_amr, write_punctures);
    }
}

// Things to do before a plot level - need to calculate the Weyl scalars
void BinaryBHLevel::prePlotLevel()
{
    fillAllGhosts();
    if (m_p.activate_extraction == 1)
    {
        BoxLoops::loop(
            make_compute_pack(Weyl4(m_p.extraction_params.center, m_dx),
                              Constraints(m_dx)),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    }
}
