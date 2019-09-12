/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to LICENSE, in Chombo's root directory.
 */
#endif

// General includes:
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>

#include "parstream.H" //Gives us pout()
using std::endl;
#include "GRAMR.hpp"

#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "AMRInterpolator.hpp"
#include "DefaultLevelFactory.hpp"
#include "InterpolationQuery.hpp"
#include "InterpolatorTestLevel.hpp"
#include "Lagrange.hpp"
#include "UserVariables.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

int runInterpolatorTest(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    std::string in_string = argv[argc - 1];
    pout() << in_string << std::endl;
    char const *in_file = argv[argc - 1];
    GRParmParse pp(0, argv + argc, NULL, in_file);
    SimulationParameters sim_params(pp);

    GRAMR gr_amr;
    DefaultLevelFactory<InterpolatorTestLevel> interpolator_test_level_fact(
        gr_amr, sim_params);
    setupAMRObject(gr_amr, interpolator_test_level_fact);

    // Setup the AMRInterpolator
    int num_points = 2;
    pp.query("num_points", num_points);

    std::unique_ptr<double[]> chi_ptr{new double[num_points]};
    std::unique_ptr<double[]> phi_ptr{new double[num_points]};
    std::unique_ptr<double[]> interp_x{new double[num_points]};
    std::unique_ptr<double[]> interp_y{new double[num_points]};
    std::unique_ptr<double[]> interp_z{new double[num_points]};

    double L;
    pp.get("L", L);
    double extract_center[3] = {L / 2, L / 2, L / 2};
    double extract_radius = L / 4;

    for (int iPhi = 0; iPhi < num_points; ++iPhi)
    {
        double extract_angle = iPhi * M_PI / num_points;
        interp_x[iPhi] =
            extract_center[0] + extract_radius * cos(extract_angle);
        interp_y[iPhi] =
            extract_center[1] + extract_radius * sin(extract_angle);
        interp_z[iPhi] = extract_center[2];
    }

    InterpolationQuery query(num_points);
    query.setCoords(0, interp_x.get())
        .setCoords(1, interp_y.get())
        .setCoords(2, interp_z.get())
        .addComp(c_chi, chi_ptr.get())
        .addComp(c_phi, phi_ptr.get());

    auto dx_scalar = GRAMRLevel::gr_cast(gr_amr.getAMRLevels()[0])
                         ->get_dx(); // coarsest grid spacing
    std::array<double, 3> dx;
    dx.fill(dx_scalar);
    std::array<double, CH_SPACEDIM> origin;
    origin.fill(dx_scalar / 2);

    AMRInterpolator<Lagrange<4>> interpolator(gr_amr, origin, dx, 0);
    interpolator.interp(query);

    int status = 0;

    for (int ipoint = 0; ipoint < num_points; ++ipoint)
    {
        status |= (abs(chi_ptr[ipoint] - 42.) > 1e-10);
        status |= (abs(phi_ptr[ipoint] - 42.) > 1e-10);
    }

    return status;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runInterpolatorTest(argc, argv);

    if (status == 0)
        pout() << "BasicAMRInterpolator test passed." << endl;
    else
        pout() << "BasicAMRInterpolator test failed with return code " << status
               << endl;

    mainFinalize();
    return status;
}
