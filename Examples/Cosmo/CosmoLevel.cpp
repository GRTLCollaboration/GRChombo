/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "CosmoLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4RHS.hpp"

// For constraints calculation
#include "NewMatterConstraints.hpp"

// For tag cells
#include "CosmoHamTaggingCriterion.hpp"

// Problem specific includes
#include "AMRReductions.hpp"
#include "ComputePack.hpp"
#include "CosmoDiagnostics.hpp"
#include "CosmoMovingPunctureGauge.hpp"
#include "GammaCalculator.hpp"
#include "InitialK.hpp"
#include "InitialScalarData.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"

// For linout
#include "CustomExtraction.hpp"

// Things to do at each advance step, after the RK4 is calculated
void CosmoLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance"),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void CosmoLevel::initialData()
{
    CH_TIME("CosmoLevel::initialData");
    if (m_verbosity)
        pout() << "CosmoLevel::initialData " << m_level << endl;

    // First set everything to zero then initial conditions for scalar field -
    // Set initial condition of inflaton, see details in Potential.hpp and
    // InitialScalarData.hpp

    BoxLoops::loop(
        make_compute_pack(SetValue(0.),
                          InitialScalarData(m_p.initial_params, m_dx, m_p.L,
                                            m_p.scalar_field_mode)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);

    // Set initial K = -sqrt(24 pi <rho>)
    Potential potential(m_p.potential_params, m_p.L, m_p.scalar_field_mode);
    ScalarFieldWithPotential scalar_field(potential);
    InitialK<ScalarFieldWithPotential> my_initial_K(scalar_field, m_dx,
                                                    m_p.G_Newton);
    BoxLoops::loop(my_initial_K, m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
    // Calculate constraints and some diagnostics as we need it in tagging
    // criterion
    BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(
                       scalar_field, m_dx, m_p.G_Newton, c_Ham,
                       Interval(c_Mom, c_Mom), c_Ham_abs_sum,
                       Interval(c_Mom_abs_sum, c_Mom_abs_sum)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    CosmoDiagnostics<ScalarFieldWithPotential> cosmo_diagnostics(
        scalar_field, m_dx, m_p.G_Newton);
    BoxLoops::loop(cosmo_diagnostics, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
    // Assign initial rho_mean here
    // from rho = 1/2 m^2 phi^2;
    // phi0 = A sin(2 pi n x /L);
    // mean of phi0^2 = A^2 sin^2 = 0.5*A^2;
    // m = m_mode
    m_cosmo_amr.set_rho_mean(
        0.5 * m_p.scalar_field_mode * m_p.scalar_field_mode *
        (0.5 * m_p.initial_params.amplitude * m_p.initial_params.amplitude));
}

#ifdef CH_USE_HDF5
// Things to do before outputting a checkpoint file
void CosmoLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params, m_p.L, m_p.scalar_field_mode);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(
        MatterConstraints<ScalarFieldWithPotential>(
            scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom, c_Mom)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    CosmoDiagnostics<ScalarFieldWithPotential> cosmo_diagnostics(
        scalar_field, m_dx, m_p.G_Newton);
    BoxLoops::loop(cosmo_diagnostics, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
}
#endif

// Things to do in RHS update, at each RK4 step
void CosmoLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    Potential potential(m_p.potential_params, m_p.L, m_p.scalar_field_mode);
    ScalarFieldWithPotential scalar_field(potential);
    CosmoMovingPunctureGauge cosmo_moving_puncture_gauge(m_p.ccz4_params);
    cosmo_moving_puncture_gauge.set_K_mean(m_cosmo_amr.get_K_mean());
    if (m_p.max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, CosmoMovingPunctureGauge,
                      FourthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, CosmoMovingPunctureGauge,
                      SixthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// Things to do at ODE update, after soln + rhs
void CosmoLevel::specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void CosmoLevel::preTagCells()
{
    // Pre tagging - fill ghost cells and calculate Ham terms
    fillAllGhosts();
    Potential potential(m_p.potential_params, m_p.L, m_p.scalar_field_mode);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(
                       scalar_field, m_dx, m_p.G_Newton, c_Ham,
                       Interval(c_Mom, c_Mom), c_Ham_abs_sum,
                       Interval(c_Mom_abs_sum, c_Mom_abs_sum)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    CosmoDiagnostics<ScalarFieldWithPotential> cosmo_diagnostics(
        scalar_field, m_dx, m_p.G_Newton);
    BoxLoops::loop(cosmo_diagnostics, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
}

void CosmoLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    BoxLoops::loop(CosmoHamTaggingCriterion(m_dx, m_p.center_tag, m_p.rad,
                                            m_cosmo_amr.get_rho_mean()),
                   current_state_diagnostics, tagging_criterion);
}
void CosmoLevel::specificPostTimeStep()
{
    int min_level = 0;
    bool calculate_diagnostics = at_level_timestep_multiple(min_level);
    bool first_step = (m_time == 0.);

    // No need to evaluate the diagnostics more frequently than every coarse
    // timestep, but must happen on every level (not just level zero or data
    // will not be populated on finer levels)

    if (calculate_diagnostics)
    {
        fillAllGhosts();
        Potential potential(m_p.potential_params, m_p.L, m_p.scalar_field_mode);
        ScalarFieldWithPotential scalar_field(potential);
        BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(
                           scalar_field, m_dx, m_p.G_Newton, c_Ham,
                           Interval(c_Mom, c_Mom), c_Ham_abs_sum,
                           Interval(c_Mom_abs_sum, c_Mom_abs_sum)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        CosmoDiagnostics<ScalarFieldWithPotential> cosmo_diagnostics(
            scalar_field, m_dx, m_p.G_Newton);
        BoxLoops::loop(cosmo_diagnostics, m_state_new, m_state_diagnostics,
                       EXCLUDE_GHOST_CELLS);

        if (m_level == min_level)
        {
            // AMRReductions for diagnostic variables
            AMRReductions<VariableType::diagnostic> amr_reductions_diag(
                m_cosmo_amr);
            double phys_vol = amr_reductions_diag.sum(c_sqrt_gamma);
            double L2_Ham = amr_reductions_diag.norm(c_Ham);
            double L2_Mom = amr_reductions_diag.norm(c_Mom);
            double K_total = amr_reductions_diag.sum(c_K_scaled);
            m_cosmo_amr.set_rho_mean(amr_reductions_diag.sum(c_rho_scaled) /
                                     phys_vol);
            m_cosmo_amr.set_S_mean(amr_reductions_diag.sum(c_S_scaled) /
                                   phys_vol);
            m_cosmo_amr.set_K_mean(K_total / phys_vol);

            // AMRReductions for evolution variables
            AMRReductions<VariableType::evolution> amr_reductions_evo(
                m_cosmo_amr);

            double chi_mean = amr_reductions_evo.sum(c_chi) / phys_vol;

            // Write output file
            SmallDataIO constraints_file(m_p.data_path + "data_out", m_dt,
                                         m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            constraints_file.remove_duplicate_time_data();
            if (first_step)
            {
                constraints_file.write_header_line(
                    {"L^2_Ham", "L^2_Mom", "<chi>", "<rho>", "<K>"});
            }
            constraints_file.write_time_data_line({L2_Ham, L2_Mom, chi_mean,
                                                   m_cosmo_amr.get_rho_mean(),
                                                   m_cosmo_amr.get_K_mean()});

            // Use AMR Interpolator and do lineout data extraction
            // set up an interpolator
            // pass the boundary params so that we can use symmetries if
            // applicable
            AMRInterpolator<Lagrange<4>> interpolator(
                m_cosmo_amr, m_p.origin, m_p.dx, m_p.boundary_params,
                m_p.verbosity);

            // this should fill all ghosts including the boundary ones according
            // to the conditions set in params.txt
            interpolator.refresh();

            // set up the query and execute it
            std::array<double, CH_SPACEDIM> extr_point = {
                0., m_p.L / 2, m_p.L / 2}; // specified point {x \in [0,L],y \in
                                           // [0,L], z \in [0,L]}
            CustomExtraction extraction(c_rho, m_p.lineout_num_points, m_p.L,
                                        extr_point, m_dt, m_time);
            extraction.execute_query(&interpolator, m_p.data_path + "lineout");
        }
    }
}