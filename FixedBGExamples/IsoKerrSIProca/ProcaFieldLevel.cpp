/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ProcaFieldLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"

// For tag cells
#include "FixedGridsTaggingCriterionBH.hpp"

// Problem specific includes
#include "ExcisionProcaDiagnostics.hpp"
#include "ExcisionProcaEvolution.hpp"
#include "FixedBGDensityAndAngularMom.hpp"
#include "FixedBGEnergyAndAngularMomFlux.hpp"
#include "FixedBGEvolution.hpp"
#include "FixedBGProcaField.hpp"
#include "FluxExtraction.hpp"
#include "InitialConditions.hpp"
#include "IsotropicKerrFixedBG.hpp"
#include "Potential.hpp"
#include "XSquared.hpp"
//#include "ProcaConstraint.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ProcaFieldLevel::specificAdvance()
{
    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(m_dx), m_state_new, m_state_new,
                       SKIP_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void ProcaFieldLevel::initialData()
{
    CH_TIME("ProcaFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ProcaFieldLevel::initialData " << m_level << endl;

    // Set the ICs
    IsotropicKerrFixedBG kerr_bh(m_p.bg_params, m_dx); // just calculates chi
    InitialConditions set_field(m_p.field_amplitude, m_p.potential_params.mass,
                                m_p.center, m_p.bg_params, m_dx);
    BoxLoops::loop(set_field, m_state_new, m_state_new, FILL_GHOST_CELLS);
    BoxLoops::loop(kerr_bh, m_state_new, m_state_diagnostics, SKIP_GHOST_CELLS);

    // make excision data zero
    BoxLoops::loop(ExcisionProcaEvolution<ProcaField, IsotropicKerrFixedBG>(
                       m_dx, m_p.center, kerr_bh, 1.0),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS,
                   disable_simd());
}

// Things to do after each timestep
void ProcaFieldLevel::specificPostTimeStep()
{
    bool first_step = (m_time == m_dt);

    // At any level, but after the 2nd coarsest timestep
    int min_level = 1;
    bool calculate_fluxes = at_level_timestep_multiple(min_level);
    if (calculate_fluxes)
    {
        // calculate the density of the PF, but excise the BH region completely
        fillAllEvolutionGhosts();
        Potential potential(m_p.potential_params);
        ProcaField proca_field(potential, m_p.proca_damping);
        IsotropicKerrFixedBG kerr_bh(m_p.bg_params, m_dx);
        FixedBGDensityAndAngularMom<ProcaField, IsotropicKerrFixedBG> densities(
            proca_field, kerr_bh, m_dx, m_p.center);
        FixedBGEnergyAndAngularMomFlux<ProcaField, IsotropicKerrFixedBG> fluxes(
            proca_field, kerr_bh, m_dx, m_p.center,
            m_p.extraction_params.zaxis_over_xaxis);
        //        XSquared set_xsquared(m_p.potential_params, m_p.bg_params,
        //        m_p.center,
        //                              m_dx);
        BoxLoops::loop(make_compute_pack(densities, fluxes), //, set_xsquared),
                       m_state_new, m_state_diagnostics, SKIP_GHOST_CELLS);
        BoxLoops::loop(
            ExcisionProcaDiagnostics<ProcaField, IsotropicKerrFixedBG>(
                m_dx, m_p.center, kerr_bh, m_p.r_min, m_p.r_max),
            m_state_new, m_state_diagnostics, SKIP_GHOST_CELLS, disable_simd());
    }

    // write out the integral after each coarse but one timestep
    if (m_level == min_level)
    {
        // integrate the densities and write to a file
        AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
        double rho_sum = amr_reductions.sum(c_rho);
        double rho_J_sum = amr_reductions.sum(c_rhoJ);

        SmallDataIO integral_file(m_p.integral_filename, m_dt, m_time,
                                  m_restart_time, SmallDataIO::APPEND,
                                  first_step);
        // remove any duplicate data if this is post restart
        integral_file.remove_duplicate_time_data();
        std::vector<double> data_for_writing = {rho_sum, rho_J_sum};
        // write data
        if (first_step)
        {
            integral_file.write_header_line({"rho", "rhoJ"});
        }
        integral_file.write_time_data_line(data_for_writing);

        // Now refresh the interpolator and do the interpolation
        bool fill_ghosts = false;
        m_gr_amr.m_interpolator->refresh(fill_ghosts);
        m_gr_amr.fill_multilevel_ghosts(VariableType::diagnostic,
                                        Interval(c_Edot, c_Jdot), min_level);
        FluxExtraction my_extraction(m_p.extraction_params, m_dt, m_time,
                                     m_restart_time);
        my_extraction.execute_query(m_gr_amr.m_interpolator);
    }
}

// Things to do before outputting a plot file
void ProcaFieldLevel::prePlotLevel() {}

// Things to do in RHS update, at each RK4 step
void ProcaFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                      const double a_time)
{
    // Calculate MatterCCZ4 right hand side with matter_t = ProcaField
    Potential potential(m_p.potential_params);
    ProcaField proca_field(potential, m_p.proca_damping);
    IsotropicKerrFixedBG kerr_bh(m_p.bg_params, m_dx);
    FixedBGEvolution<ProcaField, IsotropicKerrFixedBG> my_matter(
        proca_field, kerr_bh, m_p.sigma, m_dx, m_p.center);
    BoxLoops::loop(my_matter, a_soln, a_rhs, SKIP_GHOST_CELLS);

    // excise within horizon for evolution vars
    BoxLoops::loop(ExcisionProcaEvolution<ProcaField, IsotropicKerrFixedBG>(
                       m_dx, m_p.center, kerr_bh, 0.9),
                   a_soln, a_rhs, SKIP_GHOST_CELLS, disable_simd());
}

void ProcaFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                              const FArrayBox &current_state)
{
    const double radius_bh = 0.5;
    BoxLoops::loop(FixedGridsTaggingCriterionBH(m_dx, m_level, m_p.max_level,
                                                m_p.L, m_p.center, radius_bh),
                   current_state, tagging_criterion, disable_simd());
}
