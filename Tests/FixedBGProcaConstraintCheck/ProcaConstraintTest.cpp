/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef _OPENMP
#include <omp.h>
#endif

// General Chombo includes
#include "BoxIterator.H"
#include "BoxLoops.hpp"
#include "Cell.hpp"
#include "ComputePack.hpp"
#include "DebuggingTools.hpp"
#include "FArrayBox.H"
#include "SetValue.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sys/time.h>

// Problem specific for tests
#include "ExcisionProcaTest.hpp"
#include "FixedBGEvolution.hpp"
#include "FixedBGProcaConstraintTest.hpp"
#include "FixedBGProcaFieldTest.hpp"
#include "InitialConditions.hpp"
#include "KerrSchildFixedBG.hpp"
#include "Potential.hpp"
#include "UserVariables.hpp"

/*  This is a test to make sure that the Proca field equations are working
 * properly. It takes rnd-ish intial data and runs one time-step and makes sure
 * that constraints are still fulfilled. This is a necessary but not sufficient
 * test!
 */

int main()
{
#ifdef _OPENMP
    std::cout << "#threads = " << omp_get_max_threads() << std::endl;
#endif

    // setup grid and outputs etc
    const int failed = 0;
    const int success = 1;
    const bool debug_plots_on = true;
    const int num_resolutions = 3;
    const int N_MIN = 24;
    const double length = 12.0;
    const double center = length / 2.0;
    const std::array<double, CH_SPACEDIM> center_vector = {center, center,
                                                           center};
    const double dt_multiplier = 0.2;
    // setup a vector of norms for checking convergence
    std::array<double, num_resolutions> error_norms;
    error_norms.fill(0.0);

    // Just make it easier to write the Proca Field class
    // by using an alias
    typedef FixedBGProcaFieldTest<Potential> ProcaField;
    // Proca field potential params
    const double proca_mass = 0.5;
    const double proca_self_interaction = 0.5;
    const double proca_damping = 0.1;
    Potential::params_t potential_params;
    potential_params.proca_mass = proca_mass;
    potential_params.proca_self_interaction = proca_self_interaction;
    Potential potential(potential_params);
    ProcaField proca_field(proca_mass, proca_damping, potential);
    const double sigma = 0.0; // kreiss oliger dissipation off for test

    // metric background
    KerrSchildFixedBG::params_t kerr_params;
    kerr_params.mass = 0.5;
    kerr_params.spin = 0.25;
    kerr_params.center = center_vector;

    // loops over resolutions
    for (int ires = 0; ires < num_resolutions; ires++)
    {
        // set up the array boxes for the vars inputs/outputs
        const int N_GRID = N_MIN * pow(2, ires);
        std::cout << " Calculating error norm with " << N_GRID << "^3 points "
                  << std::endl;
        const double dx = length / (N_GRID);
        const double dt = dt_multiplier * dx;
        Box box(IntVect(0, 0, 0), IntVect(N_GRID - 1, N_GRID - 1, N_GRID - 1));
        Box ghosted_box(IntVect(-3 - 3 * pow(2, ires), -3 - 3 * pow(2, ires),
                                -3 - 3 * pow(2, ires)),
                        IntVect(N_GRID + 2 + 3 * pow(2, ires),
                                N_GRID + 2 + 3 * pow(2, ires),
                                N_GRID + 2 + 3 * pow(2, ires)));
        Box double_ghosted_box(IntVect(-6 - 6 * pow(2, ires),
                                       -6 - 6 * pow(2, ires),
                                       -6 - 6 * pow(2, ires)),
                               IntVect(N_GRID + 5 + 3 * pow(2, ires),
                                       N_GRID + 5 + 3 * pow(2, ires),
                                       N_GRID + 5 + 3 * pow(2, ires)));
        FArrayBox fixedbg_data(double_ghosted_box, NUM_VARS);
        FArrayBox fixedbg_rhs_data(ghosted_box, NUM_VARS);
        FArrayBox constraint_data(box, NUM_VARS);

        // Initalising stuff with zeros
        BoxLoops::loop(make_compute_pack(SetValue(0.0)), fixedbg_rhs_data,
                       fixedbg_rhs_data);
        BoxLoops::loop(make_compute_pack(SetValue(0.0)), constraint_data,
                       constraint_data);
        BoxLoops::loop(make_compute_pack(SetValue(0.0)), fixedbg_data,
                       fixedbg_data);

        // Setting up inital data for the field - for now only A_i is non zero
        // because this satisfies the constraints trivially
        KerrSchildFixedBG kerr_metric(kerr_params, dx);
        BoxLoops::loop(
            InitialConditions<KerrSchildFixedBG>(
                kerr_metric, proca_self_interaction, length, dx, center_vector),
            fixedbg_data, fixedbg_data); //, disable_simd());

        // for each resolution the timestep is halved
        for (int itime = 0; itime < pow(2, ires); itime++)
        {
            // Calculating the RHS for the scalar
            FixedBGEvolution<ProcaField, KerrSchildFixedBG> my_matter(
                proca_field, kerr_metric, sigma, dx, center_vector);
            // set other non evolution vars to zero
            SetValue set_constraints_zero(0.0, Interval(c_chi, NUM_VARS - 1));
            auto compute_pack =
                make_compute_pack(my_matter, set_constraints_zero);
            BoxLoops::loop(compute_pack, fixedbg_data, fixedbg_rhs_data);
            // Make a single timestep var(t+dt) = var + dvardt * dt
            fixedbg_rhs_data *= dt;
            fixedbg_data += fixedbg_rhs_data;
        }

        // Calculating gauss constraint after the timestep
        // excise region in BH as usually big errors there
        FixedBGProcaConstraintTest<Potential, KerrSchildFixedBG> my_constraint(
            kerr_metric, dx, proca_mass, proca_damping, potential);
        ExcisionProcaTest<ProcaField, KerrSchildFixedBG> excision(
            dx, center_vector, kerr_metric);
        BoxLoops::loop(my_constraint, fixedbg_data, constraint_data);
        BoxLoops::loop(excision, fixedbg_data, constraint_data, disable_simd());

        // Checking the results
        const int max_norm = 0;
        const int L1_norm = 1;
        const int num_comps = 1;

        // save the L1 norm for the convergence check - gives sum of abs
        // values if L1_norm is chosen, so divide by number of cells for average
        error_norms[ires] =
            constraint_data.norm(L1_norm, c_gauss, num_comps) * pow(N_GRID, -3);
        std::cout << "the error in the Proca constraint is : ";
        std::cout << error_norms[ires] << std::endl;

        // Output slice of data on top res, useful for debugging
        // activate by setting debug_plots to true above
        if (ires == num_resolutions - 1 && debug_plots_on)
        {
            std::string filename = "output.txt";
            std::ofstream outfile;
            outfile.clear();
            outfile.open(filename);
            outfile << std::setw(20) << "# x, y, z, phi, E^x, A^x, Constraint"
                    << "\n";
            BoxIterator bit(box);
            for (bit.begin(); bit.ok(); ++bit)
            {
                // work out location on the grid, taking slice through center
                IntVect iv = bit();
                if (iv[1] == N_GRID / 2)
                {
                    double x = dx * (iv[0] + 0.5) - center;
                    double y = dx * (iv[1] + 0.5) - center;
                    double z = dx * (iv[2] + 0.5) - center;
                    double out1 = fixedbg_data(iv, c_phi);
                    double out2 = fixedbg_data(iv, c_Evec1);
                    double out3 = fixedbg_data(iv, c_Avec1);
                    double out4 = constraint_data(iv, c_gauss);

                    outfile << std::setw(20) << x << std::setw(20) << y;
                    outfile << std::setw(20) << z;
                    outfile << std::setw(20) << out1 << std::setw(20) << out2;
                    outfile << std::setw(20) << out3 << std::setw(20) << out4;
                    outfile << "\n";
                }
            }
            outfile.close();
        }
    }

    // check convergence
    for (int ires = 0; ires < num_resolutions - 1; ires++)
    {
        if ((abs(error_norms[ires + 1]) - abs(error_norms[ires])) > 0)
        {
            std::cout << "ERROR DOESN'T REDUCE (NB if very small could be "
                         "numerical noise)"
                      << std::endl;
            return failed;
        }
    }
    std::cout << "OK, ERROR REDUCES, TEST PASSED" << std::endl;
    return success;
}
