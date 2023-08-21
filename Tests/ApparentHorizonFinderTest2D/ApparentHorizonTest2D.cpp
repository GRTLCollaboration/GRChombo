/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __    _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to LICENSE, in Chombo's root directory.
 */
#endif
/*
The purpose of this test is to test:
    1) if the AHFinder code remains compatible with 2D
    2) if GRChombo is properly setup also to work with 2D (Chombo should
recompile again for 2D when compiling this example)
*/

// This test is valid for the below value of GR_SPACEDIM
#define GR_SPACEDIM 2

// General includes:
#include <iostream>

#include "GRAMR.hpp"
#include "parstream.H" //Gives us pout()

#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "AMRInterpolator.hpp"
#include "ApparentHorizonTest2DLevel.hpp"
#include "Lagrange.hpp"
#include "SmallDataIO.hpp"
#include "SmallDataIOReader.hpp"

#ifdef USE_AHFINDER
#include "AHFinder.hpp"
#endif

int runApparentHorizonTest2D(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameters class
    // To add more parameters edit the SimulationParameters file.
    std::string in_string = argv[argc - 1];
    pout() << in_string << std::endl;
    GRParmParse pp(0, argv + argc, NULL, in_string.c_str());
    SimulationParameters sim_params(pp);

    GRAMR gr_amr;
    DefaultLevelFactory<ApparentHorizonTest2DLevel> ah_test_level_fact(
        gr_amr, sim_params);
    setupAMRObject(gr_amr, ah_test_level_fact);

    int status = 0;

#ifdef USE_AHFINDER
    // Set up interpolator and PETSc subcommunicator when AH extraction is
    // active
    AMRInterpolator<Lagrange<4>> interpolator(
        gr_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    gr_amr.set_interpolator(&interpolator);

    AHFinder<AHSurfaceGeometry, AHFunction> ah_finder;
    AHStringGeometry sph(sim_params.L);

    // use 2. as initial guess, doesn't matter much (as long as it converges to
    // the correct solution)
    ah_finder.set_interpolator(&interpolator);
    ah_finder.add_ah(sph, 2., sim_params.AH_params);

    if (!ah_finder.get(0)->get_converged())
        status = 3;
    else
    {
        double area = 0;
        SmallDataIOReader file_reader;
        file_reader.open("stats_AH1.dat");
        // get area from file to determine status
        if (!file_reader.contains_data())
        {
            status = 2;
        }
        else
        {
            auto area_col = file_reader.get_column(2);
            area = area_col[0];
        }
        if (status == 0)
        {
            // Exact value from Mathematica
            /*
            Solve[y + A * Sin[2 Pi x / L] == Pi, y][[1]]
            NIntegrate[ Sqrt[ 1 + D[y/.%,x]^2 ], {x, 0, L} ]

            For A = 1 and L = 16
            */
            double area_exact = 16.600072;

            double error_perc = fabs(1. - area / area_exact) * 100;
            pout() << "error = " << error_perc << "%" << std::endl;
            status |= (error_perc > 0.1) ||
                      +std::isnan(
                          error_perc); // accept 0.1% error in area calculation
        }
    }
#endif

    return status;
}

int main(int argc, char *argv[])
{
#ifdef USE_AHFINDER
    mainSetup(argc, argv);

    int status = runApparentHorizonTest2D(argc, argv);

#else
    int status = -1;
#endif

    if (status == 0)
    {
        std::cout << "ApparentHorizon2D test passed." << endl;
        pout() << "ApparentHorizon2D test passed." << endl;
    }
    else if (status == -1)
    {
        pout() << "ApparentHorizon2D test skipped (USE_AHFINDER undefined)."
               << std::endl;
    }
    else
    {
        std::cout << "ApparentHorizon2D test FAILED with return code " << status
                  << endl;
        pout() << "ApparentHorizon2D test FAILED with return code " << status
               << endl;
    }
#ifdef USE_AHFINDER
    mainFinalize();
#endif
    return std::max(status, 0);
}