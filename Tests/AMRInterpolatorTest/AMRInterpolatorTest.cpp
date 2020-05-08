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
    const int num_points = sim_params.num_points;

    std::vector<double> A_ptr(num_points);
    std::vector<double> B_ptr(num_points);
    std::vector<double> B_dx_ptr(num_points);
    std::vector<double> interp_x(num_points);
    std::vector<double> interp_y(num_points);
    std::vector<double> interp_z(num_points);

    double extract_radius = sim_params.L / 4;

    for (int iPhi = 0; iPhi < num_points; ++iPhi)
    {
        double extract_angle = iPhi * M_PI / num_points;
        interp_x[iPhi] =
            sim_params.center[0] + extract_radius * cos(extract_angle);
        interp_y[iPhi] =
            sim_params.center[1] + extract_radius * sin(extract_angle);
        interp_z[iPhi] = sim_params.center[2];
    }

    InterpolationQuery query(num_points);
    query.setCoords(0, interp_x.data())
        .setCoords(1, interp_y.data())
        .setCoords(2, interp_z.data())
        .addComp(c_A, A_ptr.data())
        .addComp(c_B, B_ptr.data())
        .addComp(c_B, B_dx_ptr.data(), Derivative::dx);

    AMRInterpolator<Lagrange<4>> interpolator(gr_amr, sim_params.origin,
                                              sim_params.dx,
                                              sim_params.boundary_params, 0);
    interpolator.interp(query);

    int status = 0;

    for (int ipoint = 0; ipoint < num_points; ++ipoint)
    {
        double x = interp_x[ipoint] - sim_params.center[0];
        double y = interp_y[ipoint] - sim_params.center[1];
        double z = interp_z[ipoint] - sim_params.center[2];

        double value_A =
            42. + x * x +
            (y - sim_params.center[1]) * (y - sim_params.center[1]) * z * z;
        double value_B = pow(x, 3);
        double value_B_dx = 3. * pow(x, 2);

        status |= (abs(A_ptr[ipoint] - value_A) > 1e-10);
        status |= (abs(B_ptr[ipoint] - value_B) > 1e-10);
        status |= (abs(B_dx_ptr[ipoint] - value_B_dx) > 1e-10);
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
