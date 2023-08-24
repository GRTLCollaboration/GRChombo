/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETUP_FUNCTIONS_HPP_
#define SETUP_FUNCTIONS_HPP_
// This file incldues several functions that need to be called to
// set up the runs but aren't very interesting for the normal user.

// Chombo includes
#include "AMRLevelFactory.H"
#include "parstream.H" //Gives us pout()

// Other includes
#include <iostream>
using std::cerr;
using std::endl;
#include "ChomboParameters.hpp"
#include "DerivativeSetup.hpp"
#include "FilesystemTools.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "IntegrationMethodSetup.hpp"

#include "simd.hpp"

#ifdef EQUATION_DEBUG_MODE
#include "DebuggingTools.hpp"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

// Chombo namespace
#include "UsingNamespace.H"

/// This function calls MPI_Init, makes sure a parameter file is supplied etc...
void mainSetup(int argc, char *argv[]);

/// This function calls all finalisations
void mainFinalize();

/// Sets up the grid parameters, problem domain and AMR object
void setupAMRObject(AMR &gr_amr, AMRLevelFactory &a_factory);

void mainSetup(int argc, char *argv[])
{
#ifdef CH_MPI
    int mpi_requested_thread_support = MPI_THREAD_SINGLE;
#ifdef _OPENMP
    if (omp_get_max_threads() > 1)
    {
        // MR: I don't think we need MPI_THREAD_MULTIPLE but it might be safer..
        mpi_requested_thread_support = MPI_THREAD_SERIALIZED;
    }
#endif
    int mpi_provided_thread_support;
    // Start MPI
    MPI_Init_thread(&argc, &argv, mpi_requested_thread_support,
                    &mpi_provided_thread_support);
    std::string mpi_provided_thread_support_str;
    switch (mpi_provided_thread_support)
    {
    case MPI_THREAD_SINGLE:
        mpi_provided_thread_support_str = "single";
        break;
    case MPI_THREAD_FUNNELED:
        mpi_provided_thread_support_str = "funneled";
        break;
    case MPI_THREAD_SERIALIZED:
        mpi_provided_thread_support_str = "serialized";
        break;
    case MPI_THREAD_MULTIPLE:
    default:
        mpi_provided_thread_support_str = "multiple";
        break;
    }
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
        std::cout << " number_procs = " << number_procs << endl;
#ifdef _OPENMP
        std::cout << " threads = " << omp_get_max_threads() << endl;
#ifdef CH_MPI
        std::cout << " mpi provided thread support = "
                  << mpi_provided_thread_support_str << endl;
#endif
#endif
        std::cout << " simd width (doubles) = " << simd_traits<double>::simd_len
                  << endl;
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
#ifdef CH_USE_MEMORY_TRACKING
    dumpmemoryatexit();
#endif
    MPI_Finalize();
#endif
}

void setupAMRObject(GRAMR &gr_amr, AMRLevelFactory &a_factory)
{
    // Reread the params - just the base ones
    // Note that we could have passed these through the function
    // but this way preserves backwards compatibility
    GRParmParse pp;
    ChomboParameters chombo_params(pp);

    // set size of box
    Box problem_domain(IntVect::Zero, chombo_params.ivN);
    ProblemDomain physdomain(problem_domain);

    // set periodicity
    for (int dir = 0; dir < SpaceDim; dir++)
    {
        physdomain.setPeriodic(dir,
                               chombo_params.boundary_params.is_periodic[dir]);
    }

    // Define the AMR object
    gr_amr.define(chombo_params.max_level, chombo_params.ref_ratios, physdomain,
                  &a_factory);

    // The buffer defines the minimum number of level l cells there have to be
    // between level l+1 and level l-1
    // It needs to be at least ceil(num_ghosts/max_ref_ratio) for proper nesting
    gr_amr.gridBufferSize(chombo_params.grid_buffer_size);

    // set checkpoint and plot intervals and prefixes
#ifdef CH_USE_HDF5
    gr_amr.checkpointInterval(chombo_params.checkpoint_interval);
    gr_amr.checkpointPrefix(chombo_params.hdf5_path +
                            chombo_params.checkpoint_prefix);
    if (chombo_params.plot_interval != 0)
    {
        gr_amr.plotInterval(chombo_params.plot_interval);
        gr_amr.plotPrefix(chombo_params.hdf5_path + chombo_params.plot_prefix);
    }
#endif

    // Number of coarse time steps from one regridding to the next
    gr_amr.regridIntervals(chombo_params.regrid_interval);

    // max and min box sizes, fill ratio determining accuracy of regrid
    gr_amr.maxGridSize(chombo_params.max_grid_size);
    gr_amr.blockFactor(chombo_params.block_factor);
    gr_amr.fillRatio(chombo_params.fill_ratio);

    // Set verbosity
    gr_amr.verbosity(chombo_params.verbosity);

    // Set timeEps to half of finest level dt
    // Chombo sets it to 1.e-6 by default (AMR::setDefaultValues in AMR.cpp)
    // This is only not enough for >~20 levels
    double eps = 1.;
    for (int ilevel = 0; ilevel < chombo_params.max_level; ++ilevel)
        eps /= chombo_params.ref_ratios[ilevel];
    gr_amr.timeEps(std::min(1.e-6, eps / 2.));

    // Set up input files
    if (!chombo_params.restart_from_checkpoint)
    {
#ifdef CH_USE_HDF5
        if (!FilesystemTools::directory_exists(chombo_params.hdf5_path))
            FilesystemTools::mkdir_recursive(chombo_params.hdf5_path);
#endif

        gr_amr.setupForNewAMRRun();
    }
    else
    {
#ifdef CH_USE_HDF5
        HDF5Handle handle(chombo_params.restart_file, HDF5Handle::OPEN_RDONLY);
        // read from checkpoint file
        gr_amr.setupForRestart(handle);
        handle.close();
#else
        MayDay::Error("GRChombo restart only defined with hdf5");
#endif
    }
}

#endif /* SETUP_FUNCTIONS_HPP_ */
