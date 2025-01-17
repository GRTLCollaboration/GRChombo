/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldCosmoLevel.hpp"
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

// For lineout
#include "ConstraintsExtraction.hpp"
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
    double my_scalar_mass = m_p.potential_params.scalar_mass;
    BoxLoops::loop(make_compute_pack(SetValue(0.),
                                     InitialScalarData(m_p.initial_params, m_dx,
                                                       my_scalar_mass)),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();
    // Note that the GammaCaluculator is not necessary since the data is
    // conformally flat. It is left here for generality.
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);

    // Set initial K = -sqrt(24 pi <rho>)
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    InitialK<ScalarFieldWithPotential> my_initial_K(scalar_field, m_dx,
                                                    m_p.G_Newton);
    BoxLoops::loop(my_initial_K, m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);

    // Calculate constraints and some diagnostics as we need it in tagging
    // criterion
    BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(
                       scalar_field, m_dx, m_p.G_Newton, c_Ham,
                       Interval(c_Mom, c_Mom), m_p.min_chi, c_Ham_abs_sum,
                       Interval(c_Mom_abs_sum, c_Mom_abs_sum)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    CosmoDiagnostics<ScalarFieldWithPotential> cosmo_diagnostics(
        scalar_field, m_dx, m_p.G_Newton);
    BoxLoops::loop(cosmo_diagnostics, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
    // Assign initial rho_mean here
    InitialScalarData initial_scalar_data(m_p.initial_params, m_dx,
                                          my_scalar_mass);
    m_cosmo_amr.set_rho_mean(initial_scalar_data.compute_initial_rho_mean());
}

// Things to do when restarting from a checkpoint file
void CosmoLevel::postRestart()
{
    // only want to do this on the first restart and also every restart
    if (m_time == 0.0)
    {
        fillAllGhosts();
        Potential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);

        // Calculate constraints and some diagnostics as we need it in tagging
        // criterion
        BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(
                           scalar_field, m_dx, m_p.G_Newton, c_Ham,
                           Interval(c_Mom, c_Mom), m_p.min_chi, c_Ham_abs_sum,
                           Interval(c_Mom_abs_sum, c_Mom_abs_sum)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        CosmoDiagnostics<ScalarFieldWithPotential> cosmo_diagnostics(
            scalar_field, m_dx, m_p.G_Newton);
        BoxLoops::loop(cosmo_diagnostics, m_state_new, m_state_diagnostics,
                       EXCLUDE_GHOST_CELLS);

        pout() << "Setting K mean on restart at t = " << m_time << " on level "
               << m_level << endl;

        // AMRReductions for diagnostic variables
        AMRReductions<VariableType::diagnostic> amr_reductions_diagnostic(
            m_cosmo_amr);
        double phys_vol = amr_reductions_diagnostic.sum(c_sqrt_gamma);
        double K_total = amr_reductions_diagnostic.sum(c_K_scaled);

        // Set rho_mean
        m_cosmo_amr.set_rho_mean(amr_reductions_diagnostic.sum(c_rho_scaled) /
                                 phys_vol);
        pout() << "Set rho_mean = " << m_cosmo_amr.get_rho_mean()
               << " at t = " << m_time << " on initial data at level "
               << m_level << endl;
        m_cosmo_amr.set_S_mean(amr_reductions_diagnostic.sum(c_S_scaled) /
                               phys_vol);
        // Set K_mean
        m_cosmo_amr.set_K_mean(K_total / phys_vol);
        pout() << "Calculated K mean as " << m_cosmo_amr.get_K_mean()
               << " at t = " << m_time << " on restart at level " << m_level
               << endl;

        // Use AMR Interpolator and do lineout data extraction
        // pass the boundary params so that we can use symmetries
        AMRInterpolator<Lagrange<2>> interpolator(m_cosmo_amr, m_p.origin,
                                                  m_p.dx, m_p.boundary_params,
                                                  m_p.verbosity);

        // this should fill all ghosts including the boundary ones according
        // to the conditions set in params.txt
        interpolator.refresh();

        // restart works from level 0 to highest level, so want this to happen
        // last on finest level
        int write_out_level = m_p.max_level;
        if (m_level == write_out_level)
        {
            // AMRReductions for diagnostic variables
            AMRReductions<VariableType::diagnostic> amr_reductions_diagnostic(
                m_gr_amr);
            double L2_Ham = amr_reductions_diagnostic.norm(c_Ham);
            double L2_Mom = amr_reductions_diagnostic.norm(c_Mom);

            // only on rank zero write out the result
            if (procID() == 0)
            {
                pout()
                    << "The initial norm of the constraint vars on restart is "
                    << L2_Ham << " for the Hamiltonian constraint and "
                    << L2_Mom << " for the momentum constraints" << endl;
            }

            // set up the query and execute it
            int num_points = 3 * m_p.ivN[0];
            ConstraintsExtraction constraints_extraction(
                c_Ham, c_Mom, num_points, m_p.L, m_p.center, m_dt, m_time);
            constraints_extraction.execute_query(
                &interpolator, m_p.data_path + "constraints_lineout");
        }
    }
}

#ifdef CH_USE_HDF5
// Things to do before outputting a checkpoint file
void CosmoLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(
                       scalar_field, m_dx, m_p.G_Newton, c_Ham,
                       Interval(c_Mom, c_Mom), m_p.min_chi),
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
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    CosmoMovingPunctureGauge cosmo_moving_puncture_gauge(m_p.ccz4_params);
    cosmo_moving_puncture_gauge.set_K_mean(m_cosmo_amr.get_K_mean());
    MatterCCZ4RHS<ScalarFieldWithPotential, CosmoMovingPunctureGauge,
                  FourthOrderDerivatives>
        my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                       m_p.formulation, m_p.G_Newton);
    BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
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
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(
                       scalar_field, m_dx, m_p.G_Newton, c_Ham,
                       Interval(c_Mom, c_Mom), m_p.min_chi, c_Ham_abs_sum,
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
    double rho_mean = m_cosmo_amr.get_rho_mean();
    BoxLoops::loop(CosmoHamTaggingCriterion(m_dx, m_p.tagging_center,
                                            m_p.tagging_radius, rho_mean),
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
        Potential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);
        BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(
                           scalar_field, m_dx, m_p.G_Newton, c_Ham,
                           Interval(c_Mom, c_Mom), m_p.min_chi, c_Ham_abs_sum,
                           Interval(c_Mom_abs_sum, c_Mom_abs_sum)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        CosmoDiagnostics<ScalarFieldWithPotential> cosmo_diagnostics(
            scalar_field, m_dx, m_p.G_Newton);
        BoxLoops::loop(cosmo_diagnostics, m_state_new, m_state_diagnostics,
                       EXCLUDE_GHOST_CELLS);

        if (m_level == min_level)
        {
            // AMRReductions for diagnostic variables
            AMRReductions<VariableType::diagnostic> amr_reductions_diagnostic(
                m_cosmo_amr);
            double phys_vol = amr_reductions_diagnostic.sum(c_sqrt_gamma);
            double L2_Ham = amr_reductions_diagnostic.norm(c_Ham);
            double L2_Mom = amr_reductions_diagnostic.norm(c_Mom);
            double K_total = amr_reductions_diagnostic.sum(c_K_scaled);
            m_cosmo_amr.set_rho_mean(
                amr_reductions_diagnostic.sum(c_rho_scaled) / phys_vol);
            m_cosmo_amr.set_S_mean(amr_reductions_diagnostic.sum(c_S_scaled) /
                                   phys_vol);
            m_cosmo_amr.set_K_mean(K_total / phys_vol);

            // AMRReductions for evolution variables
            AMRReductions<VariableType::evolution> amr_reductions_evolution(
                m_cosmo_amr);

            double chi_mean = amr_reductions_evolution.sum(c_chi) / phys_vol;

            // Write output file
            SmallDataIO data_out_file(m_p.data_path + "data_out", m_dt, m_time,
                                      m_restart_time, SmallDataIO::APPEND,
                                      first_step);
            data_out_file.remove_duplicate_time_data();
            if (first_step)
            {
                data_out_file.write_header_line(
                    {"L^2_Ham", "L^2_Mom", "<chi>", "<rho>", "<K>"});
            }
            data_out_file.write_time_data_line({L2_Ham, L2_Mom, chi_mean,
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
            std::array<double, CH_SPACEDIM> extraction_origin = {
                0., m_p.L / 2, m_p.L / 2}; // specified point {x \in [0,L],y \in
                                           // [0,L], z \in [0,L]}
            // rho lineout
            CustomExtraction rho_extraction(c_rho, m_p.lineout_num_points,
                                            m_p.L, extraction_origin, m_dt,
                                            m_time);
            rho_extraction.execute_query(&interpolator,
                                         m_p.data_path + "rho_lineout");
        }
    }
}