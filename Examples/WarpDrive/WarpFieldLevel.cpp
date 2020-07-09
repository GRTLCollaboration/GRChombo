/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "WarpFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4.hpp"

// For constraints calculation
#include "MatterConstraints.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "SetValue.hpp"
#include "WarpBubble.hpp"
#include "WarpField.hpp"
#include "WarpMatter.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

// Things to do at each advance step, after the RK4 is calculated
void WarpFieldLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    /*    fillAllGhosts();
        // Populate the constraints with the sources set to zero
        BoxLoops::loop(SetValue(0.0, Interval(c_rho, c_S3)),
                       m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
        WarpField warp_field(m_time, m_xs, m_vs);
        BoxLoops::loop(MatterConstraints<WarpField>(
                           warp_field, m_dx, m_p.G_Newton),
                       m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
        // run part to fix sources
        BoxLoops::loop(WarpBubble2(m_p.warp_params, m_dx),
                       m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
        fillAllGhosts();
    */
    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void WarpFieldLevel::initialData()
{
    CH_TIME("WarpFieldLevel::initialData");
    if (m_verbosity)
        pout() << "WarpFieldLevel::initialData " << m_level << endl;

    // set the velocity and position
    m_xs = m_p.warp_params.bubble_center[0];
    m_vs = m_p.warp_params.warp_speed;

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then  initial conditions for Warp field - here a
    // bubble
    BoxLoops::loop(
        make_compute_pack(SetValue(0.0), WarpBubble(m_p.warp_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // run part to fix sources, NB ghosts will be filled in pre plot level
    BoxLoops::loop(WarpMatter(m_p.warp_params, m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);

    // The constraints initially
    WarpField warp_field(m_p.warpfield_params, m_time, m_xs, m_vs);
    BoxLoops::loop(MatterConstraints<WarpField>(warp_field, m_dx, m_p.G_Newton),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

// Things to do before outputting a plot file
void WarpFieldLevel::prePlotLevel()
{
    fillAllGhosts();
    WarpField warp_field(m_p.warpfield_params, m_time, m_xs, m_vs);
    BoxLoops::loop(MatterConstraints<WarpField>(warp_field, m_dx, m_p.G_Newton),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void WarpFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                     const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = WarpField
    // We don't want undefined values floating around in the constraints so
    // zero these
    WarpField warp_field(m_p.warpfield_params, m_time, m_xs, m_vs);
    MatterCCZ4<WarpField> my_ccz4_matter(warp_field, m_p.ccz4_params, m_dx,
                                         m_p.sigma, m_p.formulation,
                                         m_p.G_Newton);
    BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

// Things to do at ODE update, after soln + rhs
void WarpFieldLevel::specificUpdateODE(GRLevelData &a_soln,
                                       const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void WarpFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                             const FArrayBox &current_state)
{
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.L, m_p.center),
                   current_state, tagging_criterion, disable_simd());
}

void WarpFieldLevel::specificPostTimeStep()
{
    // Populate the Weyl Scalar values on the grid
    fillAllGhosts();
    BoxLoops::loop(Weyl4(m_p.extraction_params.extraction_center, m_dx),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    if (m_level == m_p.extraction_params.min_extraction_level())
    {
        // Now refresh the interpolator and do the interpolation
        m_gr_amr.m_interpolator->refresh();
        WeylExtraction my_extraction(m_p.extraction_params, m_dt, m_time);
        my_extraction.execute_query(m_gr_amr.m_interpolator);
    }
}
