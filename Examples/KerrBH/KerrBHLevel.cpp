/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "KerrBHLevel.hpp"
#include "BoxLoops.hpp"
#include "CCZ4.hpp"
#include "ChiTaggingCriterion.hpp"
#include "ComputePack.hpp"
#include "Constraints.hpp"
#include "KerrBHLevel.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SetValue.hpp"
#include "TraceARemoval.hpp"

// Initial data
#include "GammaCalculator.hpp"
#include "KerrBH.hpp"

void KerrBHLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, EXCLUDE_GHOST_CELLS,
                       disable_simd());
}

void KerrBHLevel::initialData()
{
    CH_TIME("KerrBHLevel::initialData");
    if (m_verbosity)
        pout() << "KerrBHLevel::initialData " << m_level << endl;

    // First set everything to zero (to avoid undefinded values on constraints)
    // then calculate initial data  Get the Kerr solution in the variables, then
    // calculate the \tilde\Gamma^i numerically as these  are non zero and not
    // calculated in the Kerr ICs
    BoxLoops::loop(
        make_compute_pack(SetValue(0.), KerrBH(m_p.kerr_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

void KerrBHLevel::preCheckpointLevel()
{
    fillAllGhosts();
    BoxLoops::loop(Constraints(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

void KerrBHLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                  const double a_time)
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side and set constraints to zero to avoid
    // undefined values
    BoxLoops::loop(make_compute_pack(CCZ4(m_p.ccz4_params, m_dx, m_p.sigma),
                                     SetValue(0, Interval(c_Ham, c_Mom3))),
                   a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

void KerrBHLevel::specificUpdateODE(GRLevelData &a_soln,
                                    const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

// Specify which variables to write at plot intervals
void KerrBHLevel::specificWritePlotHeader(std::vector<int> &plot_states) const
{
    // Specify the variables we want to output as plot
    plot_states = {c_chi, c_K, c_lapse, c_shift1};
}

void KerrBHLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                          const FArrayBox &current_state)
{
    BoxLoops::loop(ChiTaggingCriterion(m_dx), current_state, tagging_criterion);
}
