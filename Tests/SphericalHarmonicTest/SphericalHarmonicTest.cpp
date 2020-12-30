/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "FArrayBox.H"

// Other includes
#ifdef _OPENMP
#include <omp.h>
#endif

#include "BoxLoops.hpp"
#include "Cell.hpp"
#include "ComputePack.hpp"
#include "DebuggingTools.hpp"
#include "HarmonicTest.hpp"
#include "SetValue.hpp"
#include "UserVariables.hpp"
#include <sys/time.h>

// Chombo namespace
#include "UsingNamespace.H"

int main()
{
#ifdef _OPENMP
    std::cout << "#threads = " << omp_get_max_threads() << std::endl;
#endif

    const int N_GRID = 64;
    Box box(IntVect(0, 0, 0), IntVect(N_GRID - 1, N_GRID - 1, N_GRID - 1));
    FArrayBox in_fab(box, NUM_VARS);
    BoxLoops::loop(make_compute_pack(SetValue(0.0)), in_fab, in_fab);
    FArrayBox out_fab(box, NUM_VARS);
    BoxLoops::loop(make_compute_pack(SetValue(0.0)), out_fab, out_fab);
    double length = 64.0;

    const double dx = length / (N_GRID);
    const double center = length / 2.0;

    for (int iz = 0; iz < N_GRID; ++iz)
    {
        const double z = (iz + 0.5) * dx - center;
        for (int iy = 0; iy < N_GRID; ++iy)
        {
            const double y = (iy + 0.5) * dx - center;
            for (int ix = 0; ix < N_GRID; ++ix)
            {
                const double x = (ix + 0.5) * dx - center;
                double r = sqrt(x * x + y * y + z * z);
                double rho = sqrt(x * x + y * y);
                const IntVect iv(ix, iy, iz);
                if (r < 1e-6)
                {
                    r = 1e-6;
                }
                if (rho < 1e-6)
                {
                    rho = 1e-6;
                }

                // here testing the es = -1, el = 2, em = -1 case
                // and also the calculation of r in coords
                double harmonic;
                harmonic = sqrt(5.0 / 16.0 / M_PI) * x *
                           (2 * z * z - z * r - r * r) / rho / r / r;
                in_fab(iv, c_phi) = harmonic / r / r;
            }
        }
    }

    std::array<double, CH_SPACEDIM> center_vector = {center, center, center};

    // Test the spherical harmonics across grid
    BoxLoops::loop(HarmonicTest(center_vector, dx), in_fab,
                   out_fab); // disable_simd());
    out_fab -= in_fab;

    int failed = 0;

    for (int i = 0; i < NUM_VARS; ++i)
    {
        double max_err = out_fab.norm(0, i, 1);
        double max_act = in_fab.norm(0, i, 1);
        if (max_err / max_act > 1e-10)
        {
            std::cout << "COMPONENT " << UserVariables::variable_names[i]
                      << " DOES NOT AGREE: MAX ERROR = "
                      << out_fab.norm(0, i, 1) << std::endl;
            std::cout << "COMPONENT " << UserVariables::variable_names[i]
                      << " DOES NOT AGREE: MAX Actual Value = " << max_act
                      << std::endl;
            failed = -1;
        }
    }

    if (failed == 0)
        std::cout << "Spherical Harmonic test passed..." << std::endl;
    else
        std::cout << "Spherical Harmonic test failed..." << std::endl;

    return failed;
}
