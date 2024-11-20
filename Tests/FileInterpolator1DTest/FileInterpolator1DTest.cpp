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

// Chombo includes:
#include "parstream.H" //Gives us pout()

// Other includes
#include "FileInterpolator1D.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"

#include <cmath>

// Chombo namespace
#include "UsingNamespace.H"

double test_function(double x) { return sin(x * M_PI); }
double test_function_1st_der(double x) { return M_PI * cos(x * M_PI); }

int runFileInterpolator1D(int argc, char *argv[])
{
    std::string in_string = argv[argc - 1];
    pout() << in_string << std::endl;
    char const *in_file = argv[argc - 1];
    GRParmParse pp(0, argv + argc, NULL, in_file);

    FileInterpolator1D file("N", "x", "y");

    double test_point = 1. / M_PI; // value surely not in the points in the file
    double exact = test_function(test_point);
    double val = file.interpolate(test_point);

    double exact_der = test_function_1st_der(test_point);
    double val_der = file.interpolate(test_point, 1);

    double err = abs(val / exact - 1.) * 100.;
    double err_der = abs(val_der / exact_der - 1.) * 100.;

    pout() << std::setprecision(9);
    pout() << "Exact value is: " << exact << std::endl;
    pout() << "FileInterpolator1D gave: " << val << " (" << err << "% error)"
           << std::endl;
    pout() << "Exact 1st derivative value is: " << exact_der << std::endl;
    pout() << "FileInterpolator1D gave: " << val_der << " (" << err_der
           << "% error)" << std::endl;

    bool wrong = (err + err_der > 1e-2);

    return wrong;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runFileInterpolator1D(argc, argv);

    if (status == 0)
        std::cout << "FileInterpolator1D test passed." << endl;
    else
        std::cout << "FileInterpolator1D test failed with return code "
                  << status << endl;

    mainFinalize();
    return status;
}
