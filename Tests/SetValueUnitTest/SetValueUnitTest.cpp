/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "FArrayBox.H"

// Other includes
#include "BoxLoops.hpp"
#include "SetValue.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp"
#include <iostream>

// Chombo namespace
#include "UsingNamespace.H"

int main()
{
    int failed = 0;

    const int N_GRID = 8;
    Box box(IntVect(0, 0, 0), IntVect(N_GRID - 1, N_GRID - 1, N_GRID - 1));
    FArrayBox in_fab(box, NUM_VARS);

    double value = 42.;
    BoxLoops::loop(SetValue(42.), in_fab, in_fab);

    for (int iz = 0; iz < N_GRID; ++iz)
    {
        for (int iy = 0; iy < N_GRID; ++iy)
        {
            for (int ix = 0; ix < N_GRID; ++ix)
            {
                const IntVect iv(ix, iy, iz);
                if (in_fab(iv, c_chi) != value)
                {
                    failed = -1;
                    pout() << iv << std::endl;
                }
            }
        }
    }

    if (failed == 0)
        std::cout << "SetValue test passed" << std::endl;
    else
        std::cout << "SetValue test NOT passed" << std::endl;

    return failed;
}
