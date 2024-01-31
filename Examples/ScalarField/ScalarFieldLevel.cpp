/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
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
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "InitialScalarData.hpp"
#include "KerrBH.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"

//For printing out mean diagnostics
#include "AMRReductions.hpp"
#include "MeansVars.hpp"
#include <cmath>
#include <string>
#include <sstream>
#include <typeinfo>
#include <vector>
#include <fstream>

// Start time
#include <ctime>
#include <typeinfo>

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
    time_t t;
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity) { pout() << "ScalarFieldLevel::initialData " << m_level << endl; }

    //Load in data from .dat files, for h and hdot initialisation

    ifstream gw_pos;
    ifstream gw_vel;
    gw_pos.open("/nfs/st01/hpc-gr-epss/eaf49/IC-files/32-and-64-kstar/gw-re-pos-64.dat", ios::in); //open the file with the waves in it
    gw_vel.open("/nfs/st01/hpc-gr-epss/eaf49/IC-files/32-and-64-kstar/gw-re-vel-64.dat", ios::in);

    if (!gw_pos)
    {
        MayDay::Error("GW position or velocity file failed to open.");
    }

    //int m,n = 0;
    int N = m_p.initial_params.N_init;

    std::string delim = " ";
    double A = 1e-6;

    std::string p_datline;
    std::string v_datline;
    std::stringstream p_number;
    std::stringstream v_number;

    std::vector<std::vector<double> > h(std::pow(N, 3.), std::vector<double>(6)); // input array memory allocation
    std::vector<std::vector<double> > hdot(std::pow(N, 3.), std::vector<double>(6));

    int n=0; //box position counter
    for (int i=0; i < std::pow(N, 3.); i++) //
    {
        p_datline = "";
        std::getline(gw_pos, p_datline);
        int m=0; //tensor index counter

        for(int j=0; j<p_datline.length(); j++)
        {
            if(p_datline[j] != delim[0])
            {
                p_number << p_datline[j];
            }
            else
            {
                p_number >> h[n][m];
                h[n][m] *= A;
                p_number.clear();
                m++;
            }
        }

        v_datline = "";
        std::getline(gw_vel, v_datline);
        m=0;

        for(int j=0; j<v_datline.length(); j++)
        {
            if(v_datline[j] != delim[0])
            {
                v_number << p_datline[j];
            }
            else
            {
                v_number >> hdot[n][m];
                hdot[n][m] *= A;
                v_number.clear();
                m++;
            }
        }

        n++;

        p_number.clear();
        v_number.clear();

        if(n > std::pow(N, 3.))
        {
            MayDay::Error("File length has exceeded N^3.");
        }
    }

    gw_pos.close();
    gw_vel.close();

    BoxLoops::loop(
    make_compute_pack(SetValue(0.),
                        InitialScalarData(m_p.initial_params, m_dx, h, hdot)),
    m_state_new, m_state_new, INCLUDE_GHOST_CELLS,disable_simd());

    h.clear();
    hdot.clear();
    
    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

#ifdef CH_USE_HDF5
// Things to do before outputting a checkpoint file
void ScalarFieldLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(
        MatterConstraints<ScalarFieldWithPotential>(
            scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom, c_Mom), m_p.min_chi),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}
#endif

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
    if (m_p.max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      SixthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
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
    // we don't need any ghosts filled for the fixed grids tagging criterion
    // used here so don't fill any
}

void ScalarFieldLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    BoxLoops::loop(
        FixedGridsTaggingCriterion(m_dx, m_level, 2.0 * m_p.L, m_p.center),
        current_state, tagging_criterion);
}

void ScalarFieldLevel::specificPostTimeStep()
{

    CH_TIME("ScalarFieldLevel::specificPostTimeStep");

    if (!FilesystemTools::directory_exists(m_p.data_path)) 
    {
        FilesystemTools::mkdir_recursive(m_p.data_path);
    }

    bool first_step = (m_time == 0.);
    fillAllGhosts();

    double mass = m_p.potential_params.scalar_mass;//0.01;

    AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
    AMRReductions<VariableType::evolution> amr_reductions_evo(m_gr_amr);
    double vol = amr_reductions.get_domain_volume();

    BoxLoops::loop(MeansVars(m_dx, m_p.grid_params, m_p.data_path), m_state_new, m_state_diagnostics, FILL_GHOST_CELLS);

    //Calculates means
    double phibar = amr_reductions.sum(c_sf)/vol;
    double pibar = amr_reductions.sum(c_sfd)/vol;

    double a = 1./sqrt(amr_reductions.sum(c_a)/vol);
    double H = -amr_reductions.sum(c_H)/vol/3.;

    double hambar = amr_reductions.sum(c_Ham)/vol;
    double mombar = amr_reductions.sum(c_Mom)/vol;

    //Calculates energy components and the slow-roll parameters
    double kinb = 0.5*pibar*pibar;
    double potb = 0.5*mass*mass*phibar*phibar;

    double epsilon = 3.*kinb/(kinb + potb);
    double delta = 3. + mass*mass*phibar/(pibar)/(H);

    //Calculates variances
    double phivar = amr_reductions.sum(c_sf2)/vol - phibar*phibar;
    double chivar = amr_reductions.sum(c_ch2)/vol - c_a*c_a;

    //Calculates gauge quantities
    double lapse = amr_reductions_evo.sum(c_lapse)/vol;

    //Prints all that out into the data/ directory
    SmallDataIO means_file(m_p.data_path+"means_file", m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step, ".dat");
    means_file.remove_duplicate_time_data(); // removes any duplicate data from previous run (for checkpointing)

    if(first_step) 
    {
        means_file.write_header_line({"Scalar field mean","Scalar field variance","Pi mean","Scale factor","Conformal factor variance","Hubble factor",
            "Kinetic ED","Potential ED","First SRP","Second SRP","Avg Ham constr","Avg Mom constr","Avg lapse"});
    }
    means_file.write_time_data_line({phibar, phivar, pibar, a, chivar, H, kinb, potb, epsilon, delta, hambar, mombar, lapse});
}
