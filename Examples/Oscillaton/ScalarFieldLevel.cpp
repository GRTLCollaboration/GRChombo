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
#include "AMRReductions.hpp"
#include "ComputePack.hpp"
#include "ExcisionDiagnostics.hpp"
#include "FluxExtraction.hpp"
#include "GammaCalculator.hpp"
#include "MatterEnergy.hpp"
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
    const int lines = 10000002;
    const double spacing = 0.0001; // in r for the values

    std::array<double, 3> tmp = {0.0};
    std::vector<double> lapse_values, grr_values, Pi_values;

    std::string lapse_file(m_p.initial_data_prefix + "alpha001_Thomas.csv");
    ifstream ifs0(lapse_file);

    std::string grr_file(m_p.initial_data_prefix + "grr001_Thomas.csv");
    ifstream ifs1(grr_file);

    std::string Pi_file(m_p.initial_data_prefix + "Pi001_Thomas.csv");
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
    BoxLoops::loop(SetValue(0.0), m_state_new, m_state_new, FILL_GHOST_CELLS);
    BoxLoops::loop(SetValue(0.0), m_state_diagnostics, m_state_diagnostics,
                   FILL_GHOST_CELLS);
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
    BoxLoops::loop(
        MatterEnergy<ScalarFieldWithPotential>(scalar_field, m_dx, m_p.center),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    // excise within horizon
    BoxLoops::loop(
        ExcisionDiagnostics(m_dx, m_p.center, m_p.inner_r, m_p.outer_r),
        m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
        disable_simd());
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
    MatterCCZ4RHS<ScalarFieldWithPotential> my_ccz4_matter(
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

void ScalarFieldLevel::specificPostTimeStep()
{
    // At any level, but after the coarsest timestep
    double coarsest_dt = m_p.coarsest_dx * m_p.dt_multiplier;
    const double remainder = fmod(m_time, coarsest_dt);
    if (min(abs(remainder), abs(remainder - coarsest_dt)) < 1.0e-8)
    {
        // calculate the density of the PF, but excise the BH region completely
        fillAllGhosts();
        Potential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);
        BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(
                           scalar_field, m_dx, m_p.G_Newton, c_Ham,
                           Interval(c_Mom, c_Mom), c_Ham_abs_sum,
                           Interval(c_Mom_abs_sum, c_Mom_abs_sum)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        BoxLoops::loop(MatterEnergy<ScalarFieldWithPotential>(scalar_field,
                                                              m_dx, m_p.center),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        // excise within horizon
        BoxLoops::loop(
            ExcisionDiagnostics(m_dx, m_p.center, m_p.inner_r, m_p.outer_r),
            m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
            disable_simd());
    }

    // write out the integral after each coarse timestep
    if (m_level == 0)
    {
        bool first_step = (m_time == 0.0);

        // integrate the densities and write to a file
        AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
        double rho1_sum = amr_reductions.sum(c_rho1);
        double rho2_sum = amr_reductions.sum(c_rho2);
        double source1_sum = amr_reductions.sum(c_source1);
        double source2_sum = amr_reductions.sum(c_source2);
	double Ham_sum = amr_reductions.norm(c_Ham);

        SmallDataIO integral_file("VolumeIntegrals", m_dt, m_time,
                                  m_restart_time, SmallDataIO::APPEND,
                                  first_step);
        // remove any duplicate data if this is post restart
        integral_file.remove_duplicate_time_data();
        std::vector<double> data_for_writing = {rho1_sum, rho2_sum, source1_sum, source2_sum, Ham_sum};
        // write data
        if (first_step)
        {
	  integral_file.write_header_line({"rho1", "rho2", "source1", "source2", "Ham"});
        }
        integral_file.write_time_data_line(data_for_writing);

        // Now refresh the interpolator and do the interpolation
        bool fill_ghosts = false;
        m_gr_amr.m_interpolator->refresh(fill_ghosts);
        m_gr_amr.fill_multilevel_ghosts(VariableType::diagnostic,
                                        Interval(c_flux1, c_flux2));
        FluxExtraction my_extraction(m_p.extraction_params, m_dt, m_time,
                                     m_restart_time);
        my_extraction.execute_query(m_gr_amr.m_interpolator);
    }
}

void ScalarFieldLevel::preTagCells()
{

    if (m_gr_amr.s_step == 0)
    {
        // Pre tagging - fill ghost cells and calculate Ham terms
        fillAllEvolutionGhosts();
    }
    else
    {
        preTagCellsTruncationTagging();
    }

    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(
                       scalar_field, m_dx, m_p.G_Newton, c_Ham,
                       Interval(c_Mom, c_Mom), c_Ham_abs_sum,
                       Interval(c_Mom_abs_sum, c_Mom_abs_sum)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);


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
