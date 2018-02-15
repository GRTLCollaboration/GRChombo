/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETUP_FUNCTIONS_HPP_
#define SETUP_FUNCTIONS_HPP_
// This file incldues several functions that need to be called to
// set up the runs but aren't very interesting for the normal user.

#include "parstream.H" //Gives us pout()
#include <iostream>
using std::cerr;
using std::endl;
#include "AMRLevelFactory.H"
#include "GRAMR.hpp"
#include "ParmParse.H"

#ifdef EQUATION_DEBUG_MODE
#include "DebuggingTools.hpp"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/// This function calls MPI_Init, makes sure a parameter file is supplied etc...
void mainSetup(int argc, char *argv[]);

/// This function calls all finalisations
void mainFinalize();

/// Sets up the grid parameters, problem domain and AMR object
void setupAMRObject(AMR &gr_amr, AMRLevelFactory &a_factory);

void mainSetup(int argc, char *argv[])
{
#ifdef CH_MPI
    // Start MPI
    MPI_Init(&argc, &argv);
#ifdef CH_AIX
    H5dont_atexit();
#endif
// setChomboMPIErrorHandler();
#endif

    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank = 0;
    number_procs = 1;
#endif

#ifdef EQUATION_DEBUG_MODE
    EquationDebugging::check_no_omp();
    MayDay::Warning("GRChombo is running in equation debug mode. This mode is "
                    "intended only for debugging and leads to significantly "
                    "worse performance.");
#endif

    if (rank == 0)
    {
        pout() << " number_procs = " << number_procs << endl;
#ifdef _OPENMP
        pout() << " threads = " << omp_get_max_threads() << endl;
#endif
    }

    const int required_argc = 2;
    if (argc < required_argc)
    {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
    }
}

void mainFinalize()
{
#ifdef CH_MPI
    // Exit MPI
    dumpmemoryatexit();
    MPI_Finalize();
#endif
}

void setupAMRObject(GRAMR &gr_amr, AMRLevelFactory &a_factory)
{
    // Some hard-coded parameters:
    // The buffer is width of ghost cells + additional_grid_buffer
    // and defines the minimum number of level l cells there have to be
    // between level l+1 and level l-1
    const int additional_grid_buffer = 3;

    ParmParse pp;

    IntVect ivN = IntVect::Unit;
    // Setup the grid size
    for (int dir = 0; dir < SpaceDim; ++dir)
    {
        char dir_str[20];
        sprintf(dir_str, "N%d", dir + 1);
        int N;
        pp.get(dir_str, N);
        ivN[dir] = N - 1;
    }

    Box problem_domain(IntVect::Zero, ivN);
    ProblemDomain physdomain(problem_domain);

    // set periodicity
    std::vector<bool> isPeriodic;
    pp.getarr("isPeriodic", isPeriodic, 0, SpaceDim);
    for (int dir = 0; dir < SpaceDim; dir++)
    {
        physdomain.setPeriodic(dir, isPeriodic[dir]);
    }

    int max_level;
    pp.get("max_level", max_level);

    Vector<int> ref_ratios;
    pp.getarr("ref_ratio", ref_ratios, 0, max_level + 1);

    // Define the AMR object
    gr_amr.define(max_level, ref_ratios, physdomain, &a_factory);

    // To preserve proper nesting we need to know the maximum ref_ratio.
    int max_ref_ratio = ref_ratios[0];
    for (int i = 1; i < max_level + 1; ++i)
    {
        max_ref_ratio = std::max(max_ref_ratio, ref_ratios[i]);
    }

    int num_ghosts;
    pp.get("num_ghosts", num_ghosts);
    int grid_buffer_size =
        std::ceil(((double)num_ghosts) / (double)max_ref_ratio) +
        additional_grid_buffer;
    gr_amr.gridBufferSize(grid_buffer_size);

    int checkpoint_interval;
    pp.get("checkpoint_interval", checkpoint_interval);
    gr_amr.checkpointInterval(checkpoint_interval);

    // Number of coarse time steps from one regridding to the next
    Vector<int> regrid_intervals;
    pp.getarr("regrid_interval", regrid_intervals, 0, max_level + 1);
    gr_amr.regridIntervals(regrid_intervals);

    if (pp.contains("max_grid_size"))
    {
        int max_grid_size;
        pp.query("max_grid_size", max_grid_size);
        gr_amr.maxGridSize(max_grid_size);
    }

    if (pp.contains("block_factor"))
    {
        int block_factor;
        pp.query("block_factor", block_factor);
        gr_amr.blockFactor(block_factor);
    }

    if (pp.contains("fill_ratio"))
    {
        Real fill_ratio;
        pp.query("fill_ratio", fill_ratio);
        gr_amr.fillRatio(fill_ratio);
    }

    if (pp.contains("max_dt_grow"))
    {
        Real max_dt_grow;
        pp.query("max_dt_grow", max_dt_grow);
        gr_amr.maxDtGrow(max_dt_grow);
    }

    if (pp.contains("dt_tolerance_factor"))
    {
        Real dt_tolerance_factor;
        pp.query("dt_tolerance_factor", dt_tolerance_factor);
        gr_amr.dtToleranceFactor(dt_tolerance_factor);
    }

    if (pp.contains("chk_prefix"))
    {
        std::string prefix;
        pp.query("chk_prefix", prefix);
        gr_amr.checkpointPrefix(prefix);
    }

    if (pp.contains("plot_interval"))
    {
        int plot_interval;
        pp.get("plot_interval", plot_interval);
        gr_amr.plotInterval(plot_interval);
    }

    if (pp.contains("plot_prefix"))
    {
        std::string prefix;
        pp.query("plot_prefix", prefix);
        gr_amr.plotPrefix(prefix);
    }

    int verbosity;
    pp.get("verbosity", verbosity);
    gr_amr.verbosity(verbosity);

    // Set up input files
    if (!pp.contains("restart_file"))
    {
        gr_amr.setupForNewAMRRun();
    }
    else
    {
        std::string restart_file;
        pp.query("restart_file", restart_file);

#ifdef CH_USE_HDF5
        HDF5Handle handle(restart_file, HDF5Handle::OPEN_RDONLY);
        // read from checkpoint file
        gr_amr.setupForRestart(handle);
        handle.close();
#else
        MayDay::Error("GRChombo restart only defined with hdf5");
#endif
    }
}

#endif /* SETUP_FUNCTIONS_HPP_ */
