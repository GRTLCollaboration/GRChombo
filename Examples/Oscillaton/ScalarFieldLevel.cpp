/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4.hpp"

// For constraints calculation
#include "NewMatterConstraints.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"
#include "HamTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "OscillatonInitial.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
       pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // information about the csv file data
    const int lines = 100002;
    const double spacing = 0.01; // in r for the values

    std::array<double, 3> tmp = {0.0};
    std::vector<double> lapse_values, grr_values, Pi_values;

    std::string lapse_file(m_p.initial_data_prefix + "alpha001.csv");
    ifstream ifs0(lapse_file);

    std::string grr_file(m_p.initial_data_prefix + "grr001.csv");
    ifstream ifs1(grr_file);

    std::string Pi_file(m_p.initial_data_prefix + "Pi001.csv");
    ifstream ifs2(Pi_file);

    for (int i = 0; i < lines; ++i)
    {
        ifs0 >> tmp[0];
        ifs1 >> tmp[1];
        ifs2 >> tmp[2];

        lapse_values.push_back(tmp[0]);
        grr_values.push_back(tmp[1]);
        Pi_values.push_back(tmp[2]);
    }

    // Initial conditions for scalar field - Oscillaton
    BoxLoops::loop(SetValue(0.0),
                   m_state_new, m_state_new, FILL_GHOST_CELLS);
    BoxLoops::loop(SetValue(0.0),
                   m_state_diagnostics, m_state_diagnostics, FILL_GHOST_CELLS);
    BoxLoops::loop(OscillatonInitial(m_p.L, m_dx, m_p.center, spacing,
                                     lapse_values, grr_values, Pi_values),
                   m_state_new, m_state_new, FILL_GHOST_CELLS, disable_simd());

    // data is not conformally flat, so fix Gamma values
    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

// Things to do before outputting a checkpoint file
void ScalarFieldLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(
                       scalar_field, m_dx, m_p.G_Newton, c_Ham,
                       Interval(c_Mom, c_Mom), c_Ham_abs_sum,
                       Interval(c_Mom_abs_sum, c_Mom_abs_sum)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    MatterCCZ4<ScalarFieldWithPotential> my_ccz4_matter(
        scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
        m_p.G_Newton);
    BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void ScalarFieldLevel::preTagCells()
{
    // Pre tagging - fill ghost cells and calculate Ham terms
    fillAllEvolutionGhosts();
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(
                       scalar_field, m_dx, m_p.G_Newton, c_Ham,
                       Interval(c_Mom, c_Mom), c_Ham_abs_sum,
                       Interval(c_Mom_abs_sum, c_Mom_abs_sum)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    //fillAllGhosts();
}

void ScalarFieldLevel::computeDiagnosticsTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state_diagnostics)
{
    BoxLoops::loop(HamTaggingCriterion(m_dx), current_state_diagnostics,
                   tagging_criterion);
}

void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
//        Fixed grids tagging
//        BoxLoops::loop(
//            FixedGridsTaggingCriterion(m_dx, m_level, 2.0 * m_p.L,
//            m_p.center), current_state, tagging_criterion);
}
