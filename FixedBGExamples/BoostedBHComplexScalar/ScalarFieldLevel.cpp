/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"

// For RHS update
#include "BoostedBHFixedBG.hpp"
#include "FixedBGEvolution.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "ComplexPotential.hpp"
#include "ExcisionDiagnostics.hpp"
#include "ExcisionEvolution.hpp"
#include "FixedBGComplexScalarField.hpp"
#include "FixedBGEnergyConservation.hpp"
#include "FixedBGLinMomConservation.hpp"
#include "FluxExtraction.hpp"
#include "InitialScalarData.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then initial conditions for scalar field
    SetValue set_zero(0.0);
    BoostedBHFixedBG boosted_bh(m_p.bg_params, m_dx); // just claculates chi
    InitialScalarData initial_sf(m_p.initial_params);
    auto compute_pack = make_compute_pack(set_zero, boosted_bh);

    BoxLoops::loop(compute_pack, m_state_diagnostics, m_state_diagnostics,
                   SKIP_GHOST_CELLS);
    BoxLoops::loop(initial_sf, m_state_new, m_state_new, FILL_GHOST_CELLS);

    // excise within horizon, no simd
    BoxLoops::loop(
        ExcisionEvolution<ScalarFieldWithPotential, BoostedBHFixedBG>(
            m_dx, m_p.center, boosted_bh),
        m_state_new, m_state_new, SKIP_GHOST_CELLS, disable_simd());
}

void ScalarFieldLevel::specificPostTimeStep()
{
    // At any level, but after the coarsest timestep
    int min_level = 0;
    bool calculate_quantities = at_level_timestep_multiple(min_level);

    if (calculate_quantities)
    {
        fillAllGhosts();
        ComplexPotential potential(m_p.initial_params);
        ScalarFieldWithPotential scalar_field(potential);
        BoostedBHFixedBG boosted_bh(m_p.bg_params, m_dx);
        FixedBGEnergyConservation<ScalarFieldWithPotential, BoostedBHFixedBG>
            Energies(scalar_field, boosted_bh, m_dx, m_p.center);
        FixedBGLinMomConservation<ScalarFieldWithPotential, BoostedBHFixedBG>
            LinMomenta(scalar_field, boosted_bh, m_dx, m_p.center);
        BoxLoops::loop(make_compute_pack(Energies, LinMomenta), m_state_new,
                       m_state_diagnostics, SKIP_GHOST_CELLS);

        // excise within/outside specified radii, no simd
        BoxLoops::loop(
            ExcisionDiagnostics<ScalarFieldWithPotential, BoostedBHFixedBG>(
                m_dx, m_p.center, boosted_bh, m_p.inner_r, m_p.outer_r),
            m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
            disable_simd());
    }

    // write out the integral after each coarse timestep
    if (m_p.activate_extraction == 1)
    {
        if (m_level == 0)
        {
            bool first_step = (m_time == m_dt);
            // integrate the densities and write to a file
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            double rhoEnergy_sum = amr_reductions.sum(c_rhoEnergy);
            double rhoLinMom_sum = amr_reductions.sum(c_rhoLinMom);
            double sourceLinMom_sum = amr_reductions.sum(c_sourceLinMom);

            SmallDataIO integral_file("EnMomSourceIntegrals", m_dt, m_time,
                                      m_restart_time, SmallDataIO::APPEND,
                                      first_step);
            // remove any duplicate data if this is post restart
            integral_file.remove_duplicate_time_data();

            std::vector<double> data_for_writing = {
                rhoEnergy_sum, rhoLinMom_sum, sourceLinMom_sum};

            // write data
            if (first_step)
            {
                integral_file.write_header_line({"Energy density.",
                                                 "Lin. Mom. density",
                                                 "Lin. Mom. source"});
            }
            integral_file.write_time_data_line(data_for_writing);

            // Now refresh the interpolator and do the interpolation
            bool fill_ghosts = false;
            m_gr_amr.m_interpolator->refresh(fill_ghosts);
            m_gr_amr.fill_multilevel_ghosts(
                VariableType::diagnostic, Interval(c_fluxEnergy, c_fluxLinMom));
            FluxExtraction my_extraction(m_p.extraction_params, m_dt, m_time,
                                         m_restart_time);
            my_extraction.execute_query(m_gr_amr.m_interpolator);
        }
    }
}

// Things to do before a plot level
void ScalarFieldLevel::prePlotLevel() {}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Calculate right hand side with matter_t = ScalarField
    // and background_t = BoostedBH
    // RHS for non evolution vars is zero, to prevent undefined values
    ComplexPotential potential(m_p.initial_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoostedBHFixedBG boosted_bh(m_p.bg_params, m_dx);
    FixedBGEvolution<ScalarFieldWithPotential, BoostedBHFixedBG> my_evolution(
        scalar_field, boosted_bh, m_p.sigma, m_dx, m_p.center);

    BoxLoops::loop(my_evolution, a_soln, a_rhs, SKIP_GHOST_CELLS);

    // Do excision within horizon
    BoxLoops::loop(
        ExcisionEvolution<ScalarFieldWithPotential, BoostedBHFixedBG>(
            m_dx, m_p.center, boosted_bh),
        a_soln, a_rhs, SKIP_GHOST_CELLS, disable_simd());
}

// Note that for the fixed grids this only happens on the initial timestep
// simd is disabled to allow simpler use of logical operators
void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.L, m_p.center),
                   current_state, tagging_criterion);
}
