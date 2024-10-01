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

// Chombo includes
#include "parstream.H" //Gives us pout()

#include "FileInterpolator1D.hpp"
#include "SetupFunctions.hpp"
#include <vector>

// Chombo namespace
#include "UsingNamespace.H"

int runCubicSplineInterpolatorTest(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    std::string in_string = argv[argc - 1];
    pout() << in_string << std::endl;
    char const *in_file = argv[argc - 1];
    GRParmParse pp(0, argv + argc, NULL, in_file);

    FileInterpolator1D file("N", "x", "y");
    file.allow_extrapolation(false);

    // using sin(2 * pi * x) as a test function
    // notice it already satisfies the default of: f''(0) = 0 = f''(1)
    // file.set_boundary_conditions(CubicSplineInterpolator::second_deriv, 0.,
    //                              CubicSplineInterpolator::second_deriv, 0.);

    double point = 1. / M_PI; // some irrational number
    double exact = sin(2.);
    double result = file.interpolate(point);

    pout() << result << std::endl;
    pout() << exact << std::endl;

    double err = abs(result / exact - 1.) * 100.;

    pout() << std::setprecision(9);
    pout() << "Exact value is: " << exact << std::endl;
    pout() << "Cubic Spline Interpolator gave: " << result << " (" << err
           << "% error)" << std::endl;

    bool wrong = (err > 1e-3);

    return wrong;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runCubicSplineInterpolatorTest(argc, argv);

    if (status == 0)
        std::cout << "CubicSplineInterpolator test passed." << endl;
    else
        std::cout << "CubicSplineInterpolator test failed with return code "
                  << status << endl;

    mainFinalize();
    return status;
}
