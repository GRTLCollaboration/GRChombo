/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#define COMPARE_WITH_CHF
#define COVARIANTZ4

// #include <omp.h>

// Chombo includes
#include "BoxIterator.H"
#include "FArrayBox.H"

// Other includes
#include <iomanip>
#include <iostream>
#include <sys/time.h>

// Our includes
#include "BoxLoops.hpp"
#include "ConstraintTestF_F.H"
#include "NewConstraints.hpp"

// Chombo namespace
#include "UsingNamespace.H"

#define CHF_FRAn(a, n, c)                                                      \
    a.dataPtr(n),                                                              \
        D_DECL6(&a.loVect()[0], &a.loVect()[1], &a.loVect()[2],                \
                &a.loVect()[3], &a.loVect()[4], &a.loVect()[5]),               \
        D_DECL6(&a.hiVect()[0], &a.hiVect()[1], &a.hiVect()[2],                \
                &a.hiVect()[3], &a.hiVect()[4], &a.hiVect()[5]),               \
        &c

#define CHF_CONST_FRAn(a, n, c) CHF_FRAn(a, n, c)

int main()
{
    //    std::cout << "#threads = " << omp_get_max_threads() << std::endl;

    const int N_GRID = 32;
    Box box(IntVect(0, 0, 0), IntVect(N_GRID - 1, N_GRID - 1, N_GRID - 1));
    Box ghosted_box(IntVect(-3, -3, -3),
                    IntVect(N_GRID + 2, N_GRID + 2, N_GRID + 2));
    FArrayBox in_fab(ghosted_box, NUM_VARS);
    FArrayBox in_fab_cpp_result(ghosted_box, NUM_VARS);

    const double dx = 1.0 / (N_GRID - 1);

    for (int iz = -3; iz < N_GRID + 3; ++iz)
    {
        const double z = iz * dx;
        for (int iy = -3; iy < N_GRID + 3; ++iy)
        {
            const double y = iy * dx;
            for (int ix = -3; ix < N_GRID + 3; ++ix)
            {
                const double x = ix * dx;
                const IntVect iv(ix, iy, iz);

                double g[3][3];
                double g_UU[3][3];
                double chi;
                {
                    g[0][0] = 1.36778 + 2.39731 * x + 4.53541 * x * x +
                              19.9771 * x * y * y * y + 6.13801 * y * z +
                              5.65185 * z * z + 9.35842 * z * z * z * z;
                    g[0][1] = -0.07646 - 0.48786 * x - 0.75098 * x * x -
                              1.73683 * x * y * y * y + 1.71676 * y * z +
                              1.03662 * z * z + 0.35630 * z * z * z * z;
                    g[0][2] = -0.10083 + 0.12445 * x - 1.26649 * x * x -
                              1.95052 * x * y * y * y + 0.73091 * y * z -
                              1.49835 * z * z - 2.39024 * z * z * z * z;
                    g[1][1] = 0.84072 + 2.31163 * x + 3.32275 * x * x +
                              15.1662 * x * y * y * y + 8.48730 * y * z +
                              3.05098 * z * z + 17.8448 * z * z * z * z;
                    g[1][2] = -0.42495 - 0.33464 * x - 0.47012 * x * x -
                              7.38477 * x * y * y * y + 0.41896 * y * z -
                              1.36394 * z * z + 5.25894 * z * z * z * z;
                    g[2][2] = 0.60995 + 1.30428 * x + 3.86237 * x * x +
                              22.7614 * x * y * y * y + 6.93818 * y * z +
                              4.39250 * z * z + 19.0244 * z * z * z * z;
                    g[1][0] = g[0][1];
                    g[2][0] = g[0][2];
                    g[2][1] = g[1][2];

                    const double detg = g[0][0] * g[1][1] * g[2][2] +
                                        2 * g[0][1] * g[0][2] * g[1][2] -
                                        g[0][0] * g[1][2] * g[1][2] -
                                        g[1][1] * g[0][2] * g[0][2] -
                                        g[2][2] * g[0][1] * g[0][1];
                    g_UU[0][0] = (g[1][1] * g[2][2] - g[1][2] * g[1][2]) / detg;
                    g_UU[0][1] = (g[0][2] * g[1][2] - g[0][1] * g[2][2]) / detg;
                    g_UU[0][2] = (g[0][1] * g[1][2] - g[0][2] * g[1][1]) / detg;
                    g_UU[1][1] = (g[0][0] * g[2][2] - g[0][2] * g[0][2]) / detg;
                    g_UU[1][2] = (g[0][1] * g[0][2] - g[0][0] * g[1][2]) / detg;
                    g_UU[2][2] = (g[0][0] * g[1][1] - g[0][1] * g[0][1]) / detg;
                    g_UU[1][0] = g_UU[0][1];
                    g_UU[2][0] = g_UU[0][2];
                    g_UU[2][1] = g_UU[1][2];

                    chi = pow(fabs(detg), -1.0 / GR_SPACEDIM);
                    in_fab(iv, c_chi) = chi;
                    in_fab(iv, c_h11) = chi * g[0][0];
                    in_fab(iv, c_h12) = chi * g[0][1];
                    in_fab(iv, c_h13) = chi * g[0][2];
                    in_fab(iv, c_h22) = chi * g[1][1];
                    in_fab(iv, c_h23) = chi * g[1][2];
                    in_fab(iv, c_h33) = chi * g[2][2];
                }

                {
                    double K[3][3];
                    K[0][0] = -0.16238 - 0.74295 * x + 0.51595 * x * x -
                              6.60239 * x * y * y * y - 0.76401 * y * z -
                              1.81131 * z * z - 3.88228 * z * z * z * z;
                    K[0][1] = 0.15054 - 0.60088 * x - 0.15428 * x * x +
                              3.16779 * x * y * y * y - 2.00687 * y * z -
                              1.35442 * z * z - 0.67601 * z * z * z * z;
                    K[0][2] = -0.02174 - 0.36243 * x + 0.81531 * x * x +
                              4.34918 * x * y * y * y + 0.90419 * y * z -
                              0.85088 * z * z - 6.45097 * z * z * z * z;
                    K[1][1] = -0.47653 - 0.43889 * x + 0.87342 * x * x +
                              4.24684 * x * y * y * y + 0.26290 * y * z +
                              1.90095 * z * z + 3.69515 * z * z * z * z;
                    K[1][2] = 0.37472 + 0.03657 * x - 0.10327 * x * x -
                              0.95744 * x * y * y * y - 1.20800 * y * z -
                              0.43064 * z * z - 0.25419 * z * z * z * z;
                    K[2][2] = 0.34184 + 0.21495 * x - 0.73195 * x * x +
                              7.81626 * x * y * y * y + 2.48359 * y * z +
                              1.89657 * z * z - 4.10980 * z * z * z * z;
                    K[1][0] = K[0][1];
                    K[2][0] = K[0][2];
                    K[2][1] = K[1][2];

                    double trK = 0;
                    FOR(a, b) { trK += g_UU[a][b] * K[a][b]; }

                    in_fab(iv, c_K) = trK;
                    in_fab(iv, c_A11) =
                        chi * (K[0][0] - trK * g[0][0] / GR_SPACEDIM);
                    in_fab(iv, c_A12) =
                        chi * (K[0][1] - trK * g[0][1] / GR_SPACEDIM);
                    in_fab(iv, c_A13) =
                        chi * (K[0][2] - trK * g[0][2] / GR_SPACEDIM);
                    in_fab(iv, c_A22) =
                        chi * (K[1][1] - trK * g[1][1] / GR_SPACEDIM);
                    in_fab(iv, c_A23) =
                        chi * (K[1][2] - trK * g[1][2] / GR_SPACEDIM);
                    in_fab(iv, c_A33) =
                        chi * (K[2][2] - trK * g[2][2] / GR_SPACEDIM);
                }

                in_fab(iv, c_Theta) = 0.27579 + 0.25791 * x + 1.40488 * x * x +
                                      5.68276 * x * y * y * y +
                                      3.04325 * y * z + 1.81250 * z * z +
                                      1.01832 * z * z * z * z;
                in_fab(iv, c_Gamma1) =
                    -0.49482 + 0.89227 * x + 0.05571 * x * x -
                    5.38570 * x * y * y * y + 0.13979 * y * z -
                    0.68588 * z * z - 4.39964 * z * z * z * z;
                in_fab(iv, c_Gamma2) =
                    -0.09082 - 0.31017 * x + 1.06980 * x * x +
                    7.81524 * x * y * y * y - 1.65016 * y * z -
                    0.53352 * z * z - 3.20997 * z * z * z * z;
                in_fab(iv, c_Gamma3) =
                    -0.42367 + 0.03891 * x - 0.87898 * x * x +
                    6.67657 * x * y * y * y - 3.44662 * y * z -
                    0.19655 * z * z + 2.97524 * z * z * z * z;

                in_fab(iv, c_lapse) = 0.73578 + 0.36898 * x + 0.64348 * x * x +
                                      9.33487 * x * y * y * y +
                                      0.99469 * y * z + 0.20515 * z * z +
                                      8.88385 * z * z * z * z;
                in_fab(iv, c_shift1) = 0.00000 + 0.18795 * x - 0.52389 * x * x -
                                       4.14079 * x * y * y * y +
                                       0.73135 * y * z - 0.27057 * z * z +
                                       3.24187 * z * z * z * z;
                in_fab(iv, c_shift2) = 0.00000 - 0.30316 * x - 0.15184 * x * x -
                                       0.48815 * x * y * y * y +
                                       2.45991 * y * z - 0.79248 * z * z +
                                       7.14007 * z * z * z * z;
                in_fab(iv, c_shift3) = 0.00000 + 0.68835 * x - 0.52219 * x * x -
                                       7.50449 * x * y * y * y -
                                       2.35372 * y * z - 0.21476 * z * z +
                                       4.36363 * z * z * z * z;
                in_fab(iv, c_B1) = -0.26928 + 0.35045 * x - 0.48884 * x * x +
                                   2.72465 * x * y * y * y - 2.59022 * y * z -
                                   0.27384 * z * z + 0.38748 * z * z * z * z;
                in_fab(iv, c_B2) = 0.40234 + 0.26741 * x + 1.94822 * x * x -
                                   0.78276 * x * y * y * y + 2.12346 * y * z +
                                   0.69086 * z * z - 4.47639 * z * z * z * z;
                in_fab(iv, c_B3) = 0.40313 + 0.00569 * x - 1.12452 * x * x -
                                   5.49255 * x * y * y * y - 2.21932 * y * z +
                                   0.49523 * z * z + 1.29460 * z * z * z * z;
            }
        }
    }

    struct timeval begin, end;

    // Make sure instructions are hot in cache
    BoxLoops::loop(Constraints(dx, c_Ham, Interval(c_Mom1, c_Mom3)), in_fab,
                   in_fab_cpp_result, box);

    gettimeofday(&begin, NULL);

    BoxLoops::loop(Constraints(dx, c_Ham, Interval(c_Mom1, c_Mom3)), in_fab,
                   in_fab_cpp_result, box);

    gettimeofday(&end, NULL);

    int cxx_time = end.tv_sec * 1000 + end.tv_usec / 1000 -
                   begin.tv_sec * 1000 - begin.tv_usec / 1000;
    std::cout << "C++ version took " << cxx_time << "ms" << std::endl;

#ifdef COMPARE_WITH_CHF

    int SIX = 6;
    int THREE = 3;

    FORT_GETBSSNCONSTRF(
        CHF_FRA1(in_fab, c_Ham), CHF_FRAn(in_fab, c_Mom, THREE),
        CHF_CONST_FRA1(in_fab, c_chi), CHF_CONST_FRAn(in_fab, c_h, SIX),
        CHF_CONST_FRA1(in_fab, c_K), CHF_CONST_FRAn(in_fab, c_A, SIX),
        CHF_CONST_FRAn(in_fab, c_Gamma, THREE), CHF_CONST_REAL(dx),
        CHF_BOX(box));

    gettimeofday(&begin, NULL);

    FORT_GETBSSNCONSTRF(
        CHF_FRA1(in_fab, c_Ham), CHF_FRAn(in_fab, c_Mom, THREE),
        CHF_CONST_FRA1(in_fab, c_chi), CHF_CONST_FRAn(in_fab, c_h, SIX),
        CHF_CONST_FRA1(in_fab, c_K), CHF_CONST_FRAn(in_fab, c_A, SIX),
        CHF_CONST_FRAn(in_fab, c_Gamma, THREE), CHF_CONST_REAL(dx),
        CHF_BOX(box));

    gettimeofday(&end, NULL);

    int fort_time = end.tv_sec * 1000 + end.tv_usec / 1000 -
                    begin.tv_sec * 1000 - begin.tv_usec / 1000;
    std::cout << "ChF version took " << fort_time << "ms" << std::endl;

    std::cout << "C++ speedup = " << setprecision(2)
              << (double)fort_time / cxx_time << "x" << std::endl;

    in_fab_cpp_result -= in_fab;

    int error = 0;

    for (int i = c_Ham; i < NUM_VARS; ++i)
    {
        BoxIterator bit(box);
        double max_err = 0;
        IntVect location = IntVect::Zero;
        for (bit.begin(); bit.ok(); ++bit)
        {
            max_err = max(max_err, abs(in_fab_cpp_result(bit(), i)));
            if (max_err == in_fab_cpp_result(bit(), i))
                location = bit();
        }
        // double max_err = in_fab_cpp_result.norm(0, i, 1);
        if (max_err > 1e-10)
        {
            std::cout << "COMPONENT " << i
                      << " DOES NOT AGREE: MAX ERROR = " << max_err
                      << " AT INTVECT " << location << std::endl;
            error = -1;
        }
    }

    if (error == 0)
        std::cout << "Constraint cpp ChF comparison test passed" << std::endl;
    else
        std::cout << "Constraint cpp ChF comparison test NOT passed"
                  << std::endl;

#endif
    return error;
}
