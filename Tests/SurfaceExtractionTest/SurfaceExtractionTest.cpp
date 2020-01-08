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
#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"

#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "AMRInterpolator.hpp"
#include "Lagrange.hpp"
#include "SphericalExtraction.hpp"
#include "SurfaceExtractionTestLevel.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

int runSurfaceExtractionTest(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    std::string in_string = argv[argc - 1];
    pout() << in_string << std::endl;
    char const *in_file = argv[argc - 1];
    GRParmParse pp(0, argv + argc, NULL, in_file);
    SimulationParameters sim_params(pp);

    GRAMR gr_amr;
    DefaultLevelFactory<SurfaceExtractionTestLevel>
        surface_extraction_test_level_fact(gr_amr, sim_params);
    setupAMRObject(gr_amr, surface_extraction_test_level_fact);

    AMRInterpolator<Lagrange<4>> interpolator(gr_amr, sim_params.origin,
                                              sim_params.dx, 0);

    // extract chi and dphi/dx on spheres of radii specified by the extraction
    // parameters
    SphericalExtraction spherical_extraction(
        sim_params.extraction_params,
        sim_params.coarsest_dx * sim_params.dt_multiplier, 0.0, true, 0.0);
    spherical_extraction.add_var(c_chi);
    spherical_extraction.add_var(c_phi, Derivative::dx);
    spherical_extraction.extract(&interpolator);
    spherical_extraction.write_extraction("ExtractionOut_");

    // chi is the first variable so data[0] is chi
    auto integrand_chi = [](std::vector<double> data, double, double, double) {
        return data[0]; // chi
    };

    // dphi/dx is the second variable so data[1] is dphi/dx
    auto integrand_dphi_dx = [](std::vector<double> data, double, double,
                                double) {
        return data[1]; // dphi/dx
    };

    // integrate chi and dphi/dx over the extraction spheres
    auto integral_chi = spherical_extraction.integrate(
        integrand_chi, IntegrationMethod::simpson);
    auto integral_dphi_dx = spherical_extraction.integrate(
        integrand_dphi_dx, IntegrationMethod::simpson);

    int status = 0;

    for (int iradius = 0;
         iradius < sim_params.extraction_params.num_extraction_radii; ++iradius)
    {
        double r = sim_params.extraction_params.extraction_radii[iradius];
        double result_chi = integral_chi[iradius];
        double analytic_result_chi = 4.0 * M_PI * r * r * 42.0;
        double relative_error_chi = abs(result_chi / analytic_result_chi - 1.0);
        double result_dphi_dx = integral_dphi_dx[iradius];
        // analytic_result_dphi_dx = 0.0
        pout() << "At r = " << r << ", integral_chi = " << result_chi
               << ", and integral_dphi_dx = " << result_dphi_dx << "\n";
        pout() << "analytic_integral_chi = " << analytic_result_chi << "\n";
        pout() << "relative_error in integral_chi = " << relative_error_chi
               << endl;
        // want max 0.01% relative error in integral_chi
        // and 1.0e-10 absolute error in integral_dphi_dx
        // (since the latter is just a sum of zeros)
        status |= (relative_error_chi > 1.0e-4 );
        status |= (abs(result_dphi_dx) > 1.0e-10);
    }

    return status;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runSurfaceExtractionTest(argc, argv);

    if (status == 0)
        pout() << "SurfaceExtractionTest test passed." << endl;
    else
        pout() << "SurfaceExtractionTest test failed with return code "
               << status << endl;

    mainFinalize();
    return status;
}
