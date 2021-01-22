/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "BoxIterator.H"
#include "FArrayBox.H"

// Other includes
#include <iostream>

// Our includes
#include "BoxLoops.hpp"
#include "DerivativeTestsCompute.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SixthOrderDerivatives.hpp"
#include "UserVariables.hpp"

// Chombo namespace
#include "UsingNamespace.H"

bool is_wrong(double value, double correct_value, std::string deriv_type)
{
    if (abs(value - correct_value) > 1e-10)
    {
        std::cout.precision(17);
        std::cout << "Test of " << deriv_type << " failed "
                  << " with value " << value << " instad of " << correct_value
                  << ".\n";
        return true;
    }
    else
    {
        return false;
    }
}

int main()
{
    const int num_cells = 512;
    // box is flat in y direction to make test cheaper
    IntVect domain_hi_vect(num_cells - 1, 0, num_cells - 1);
    Box box(IntVect(0, 0, 0), domain_hi_vect);
    Box ghosted_box(IntVect(-4, -4, -4),
                    IntVect(num_cells + 3, 4, num_cells + 3));

    FArrayBox in_fab(ghosted_box, NUM_VARS);
    FArrayBox out_fab(box, NUM_VARS);

    const double dx = 1.0 / num_cells;

    BoxIterator bit_ghost(ghosted_box);
    for (bit_ghost.begin(); bit_ghost.ok(); ++bit_ghost)
    {
        const double x = (0.5 + bit_ghost()[0]) * dx;
        const double z = (0.5 + bit_ghost()[2]) * dx;
        for (int i = 0; i < NUM_VARS; ++i)
        {
            in_fab(bit_ghost(), i) = x * z * (z - 1);
        }
        // The dissipation component is special:
        in_fab(bit_ghost(), c_diss) = (pow(z - 0.5, 6) - 0.015625) / 720. +
                                      (z - 1) * z * pow(x, 6) / 720.;
    }

    // Fourth order derivatives
    BoxLoops::loop(DerivativeTestsCompute<FourthOrderDerivatives>(dx), in_fab,
                   out_fab);

    BoxIterator bit(box);
    for (bit.begin(); bit.ok(); ++bit)
    {
        const double x = (0.5 + bit()[0]) * dx;
        const double z = (0.5 + bit()[2]) * dx;

        bool error = false;
        error |= is_wrong(out_fab(bit(), c_d1), 2 * x * (z - 0.5),
                          "diff1 (fourth order)");
        error |= is_wrong(out_fab(bit(), c_d2), 2 * x, "diff2 (fourth order)");
        error |= is_wrong(out_fab(bit(), c_d2_mixed), 2 * (z - 0.5),
                          "mixed diff2 (fourth order)");

        double correct_dissipation = (1. + z * (z - 1)) * pow(dx, 5) / 64;
        error |= is_wrong(out_fab(bit(), c_diss), correct_dissipation,
                          "dissipation (fourth order)");

        double correct_advec_down = -2 * z * (z - 1) - 3 * x * (2 * z - 1);
        error |= is_wrong(out_fab(bit(), c_advec_down), correct_advec_down,
                          "advection down (fourth order)");

        double correct_advec_up = 2 * z * (z - 1) + 3 * x * (2 * z - 1);
        error |= is_wrong(out_fab(bit(), c_advec_up), correct_advec_up,
                          "advection up (fourth order)");

        if (error)
        {
            std::cout << "Derivative unit tests NOT passed.\n";
            return error;
        }
    }

    // Sixth order derivatives
    BoxLoops::loop(DerivativeTestsCompute<SixthOrderDerivatives>(dx), in_fab,
                   out_fab);

    for (bit.begin(); bit.ok(); ++bit)
    {
        const double x = (0.5 + bit()[0]) * dx;
        const double z = (0.5 + bit()[2]) * dx;

        bool error = false;
        error |= is_wrong(out_fab(bit(), c_d1), 2 * x * (z - 0.5),
                          "diff1 (sixth order)");
        error |= is_wrong(out_fab(bit(), c_d2), 2 * x, "diff2 (sixth order)");
        error |= is_wrong(out_fab(bit(), c_d2_mixed), 2 * (z - 0.5),
                          "mixed diff2 (sixth order)");

        double correct_dissipation = (1. + z * (z - 1)) * pow(dx, 5) / 64;
        error |= is_wrong(out_fab(bit(), c_diss), correct_dissipation,
                          "dissipation (sixth order)");

        double correct_advec_down = -2 * z * (z - 1) - 3 * x * (2 * z - 1);
        error |= is_wrong(out_fab(bit(), c_advec_down), correct_advec_down,
                          "advection down (sixth order)");

        double correct_advec_up = 2 * z * (z - 1) + 3 * x * (2 * z - 1);
        error |= is_wrong(out_fab(bit(), c_advec_up), correct_advec_up,
                          "advection up (sixth order)");

        if (error)
        {
            std::cout << "Derivative unit tests NOT passed.\n";
            return error;
        }
    }

    std::cout << "Derivative unit tests passed.\n";
    return 0;
}
