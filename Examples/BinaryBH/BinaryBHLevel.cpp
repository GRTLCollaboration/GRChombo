/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "BinaryBHLevel.hpp"
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

    // First set everything to zero (to avoid undefinded values in constraints)
    // then calculate initial data
    BoxLoops::loop(make_compute_pack(SetValue(0.), binary), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS);
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
    if (m_p.track_punctures)
    {
        const vector<double> puncture_masses = {m_p.bh1_params.mass,
                                                m_p.bh2_params.mass};
        auto puncture_coords = m_bh_amr.puncture_tracker.get_puncture_coords();
        BoxLoops::loop(ChiPunctureExtractionTaggingCriterion(
                           m_dx, m_level, m_p.max_level, m_p.extraction_params,
                           puncture_coords, m_p.activate_extraction,
                           m_p.track_punctures, puncture_masses),
                       current_state, tagging_criterion);
    }
    else
    {
        BoxLoops::loop(ChiExtractionTaggingCriterion(m_dx, m_level,
                                                     m_p.extraction_params,
                                                     m_p.activate_extraction),
                       current_state, tagging_criterion);
    }
}

void BinaryBHLevel::specificPostTimeStep()
{
    CH_TIME("BinaryBHLevel::specificPostTimeStep");
    if (m_p.activate_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            BoxLoops::loop(Weyl4(m_p.extraction_params.center, m_dx),
                           m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);

            // Do the extraction on the min extraction level
            if (m_level == min_level)
            {
                CH_TIME("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                m_gr_amr.m_interpolator->refresh();
                WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                             m_time, m_restart_time);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }

    // do puncture tracking on requested level
    if (m_p.track_punctures && m_level == m_p.puncture_tracking_level)
    {
        CH_TIME("PunctureTracking");
        // only do the write out for every coarsest level timestep
        int coarsest_level = 0;
        bool write_punctures = at_level_timestep_multiple(coarsest_level);
        m_bh_amr.puncture_tracker.execute_tracking(m_time, m_restart_time, m_dt,
                                                   write_punctures);
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
