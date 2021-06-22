#define COMPARE_WITH_CHF
#define COVARIANTZ4

// Chombo includes
#include "FArrayBox.H"

// Other includes
#ifdef _OPENMP
#include <omp.h>
#endif

#include <iomanip>
#include <iostream>
#include <sys/time.h>

#include "BoxLoops.hpp"
#include "CCZ4RHS.hpp"
#include "EMTensor.hpp"
#include "GravWavDecF_F.H"
#include "MatterWeyl4.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"
#include "UserVariables.hpp"

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
#ifdef _OPENMP
    std::cout << "#threads = " << omp_get_max_threads() << std::endl;
#endif

    const int N_GRID = 32;
    Box box(IntVect(0, 0, 0), IntVect(N_GRID - 1, N_GRID - 1, N_GRID - 1));
    Box ghosted_box(IntVect(-3, -3, -3),
                    IntVect(N_GRID + 2, N_GRID + 2, N_GRID + 2));
    FArrayBox in_fab(ghosted_box, NUM_VARS);
    FArrayBox out_fab(box, NUM_VARS);
    FArrayBox out_fab_chf(box, NUM_VARS);

    const double dx = 0.25 / (N_GRID - 1);

    for (int zz = -3; zz < N_GRID + 3; ++zz)
    {
        const double z = zz * dx;
        for (int yy = -3; yy < N_GRID + 3; ++yy)
        {
            const double y = yy * dx;
            for (int xx = -3; xx < N_GRID + 3; ++xx)
            {
                const double x = xx * dx;
                const IntVect iv(xx, yy, zz);

                double g[3][3];
                double g_UU[3][3];
                double chi;
                {
                    g[0][0] = 1. + 0.22 * sin(x * 2 * 3.14) *
                                       cos(y * 2 * 3.15) * cos(z * 2 * 3.14);
                    g[0][1] = 0.4 * sin(x * 2 * 3.14) * cos(y * 3 * 3.15) *
                              cos(z * 2 * 3.14);
                    g[0][2] = 0.1 * sin(x * 0.2 * 3.14) * cos(y * 2 * 3.15) *
                              cos(z * 2 * 3.14);
                    g[1][1] = 1. + 0.1 * sin(x * 2.1 * 3.14) *
                                       cos(y * 2 * 3.15) * cos(z * 2 * 3.14);
                    g[1][2] = 0.4 * sin(x * 1 * 3.14) * cos(y * 2 * 3.15) *
                              cos(z * 2 * 3.14);
                    g[2][2] = 1. + 0.4 * sin(x * 1.02 * 3.14) *
                                       cos(y * 3.112 * 3.15) *
                                       cos(z * 2.22 * 3.14);
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
                    in_fab(iv, c_chi2) = sqrt(chi);
                    in_fab(iv, c_h11) = chi * g[0][0];
                    in_fab(iv, c_h12) = chi * g[0][1];
                    in_fab(iv, c_h13) = chi * g[0][2];
                    in_fab(iv, c_h22) = chi * g[1][1];
                    in_fab(iv, c_h23) = chi * g[1][2];
                    in_fab(iv, c_h33) = chi * g[2][2];

                    in_fab(iv, c_Weyl4_Re) = 0.0;
                    in_fab(iv, c_Weyl4_Im) = 0.0;
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

                in_fab(iv, c_Theta) = 0;

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

                in_fab(iv, c_phi) = 0.21232 * sin(x * 2.1232 * 3.14) *
                                    cos(y * 2.5123 * 3.15) *
                                    cos(z * 2.1232 * 3.14);
                in_fab(iv, c_Pi) = 0.4112 * sin(x * 4.123 * 3.14) *
                                   cos(y * 2.2312 * 3.15) *
                                   cos(z * 2.5123 * 3.14);
                // compute the rho from other variables as this is what
                // MatterWeyl4 will do
                // in_fab(iv, c_Rho) = 0;
            }
        }
    }
    using DefaultScalarField = ScalarField<DefaultPotential>;

    // compute rho and store in FAB
    BoxLoops::loop(EMTensor<DefaultScalarField>(
                       DefaultScalarField(DefaultPotential()), dx, c_Rho),
                   in_fab, in_fab, /* can't compute in ghost cells */ box);

    Real null = 0;

    std::array<double, CH_SPACEDIM> centerGW = {0, 0, 0};

    struct timeval begin, end;
    gettimeofday(&begin, NULL);

    BoxLoops::loop(
        MatterWeyl4<DefaultScalarField>(DefaultScalarField(DefaultPotential()),
                                        centerGW, dx, CCZ4RHS<>::USE_BSSN),
        in_fab, out_fab);

    gettimeofday(&end, NULL);

    int cxx_time = end.tv_sec * 1000 + end.tv_usec / 1000 -
                   begin.tv_sec * 1000 - begin.tv_usec / 1000;
    std::cout << "C++ version took " << cxx_time << "ms" << std::endl;

#ifdef COMPARE_WITH_CHF

    int SIX = 6;
    int THREE = 3;

    gettimeofday(&begin, NULL);

    // NB I had to fudge the chi derivatives in the ChF code to get the same
    // result because of the use of chi^2 rather than chi as conformal factor.
    // Also in the calculation of R(d0,d1) the conformal Ricci tensor, we used
    // to use the calculated Gamma^i, not the evolved one. In theory they would
    // be the same, but in this test The Gamma^i data is random, so it won't
    // work, so I changed it to match the new form
    // (MR): I changed it back to reflect the changes I have made to the
    // calculation of the Ricci tensor

    FORT_GRAVWAVDEC(
        CHF_FRA1(out_fab_chf, c_Weyl4_Re), CHF_FRA1(out_fab_chf, c_Weyl4_Im),
        CHF_CONST_FRAn(in_fab, c_h, SIX),
        CHF_CONST_FRAn(in_fab, c_Gamma, THREE), CHF_CONST_FRA1(in_fab, c_chi2),
        CHF_CONST_FRA1(in_fab, c_K), CHF_CONST_FRAn(in_fab, c_A, SIX),
        CHF_CONST_FRAn(in_fab, c_shift, THREE), CHF_CONST_FRA1(in_fab, c_phi),
        CHF_CONST_FRA1(in_fab, c_Pi), CHF_CONST_FRA1(in_fab, c_lapse),
        CHF_CONST_FRA1(in_fab, c_Rho), CHF_CONST_REAL(null),
        CHF_CONST_REAL(null), CHF_CONST_REAL(dx), CHF_CONST_REAL(centerGW[0]),
        CHF_CONST_REAL(centerGW[1]), CHF_CONST_REAL(centerGW[2]), CHF_BOX(box));

    gettimeofday(&end, NULL);

    int fort_time = end.tv_sec * 1000 + end.tv_usec / 1000 -
                    begin.tv_sec * 1000 - begin.tv_usec / 1000;
    std::cout << "ChF version took " << fort_time << "ms" << std::endl;

    std::cout << "C++ speedup = " << setprecision(2)
              << (double)fort_time / cxx_time << "x" << std::endl;

    out_fab -= out_fab_chf;
    int failed = 0;
    for (int i = 0; i < NUM_VARS; ++i)
    {
        double max_err = out_fab.norm(0, i, 1);
        double max_chf = out_fab_chf.norm(0, i, 1);

        // Note that due to the change from chi' to chi = chi'^2, errors in chi
        // may be larger than other variables.
        if (max_err / max_chf > 0.00001)
        {
            std::cout << std::endl
                      << " -------------------------- Remark "
                         "---------------------------- "
                      << std::endl;
            std::cout << " --  " << UserVariables::variable_names[i]
                      << " Value = " << max_chf << std::endl;
            std::cout << " --  " << UserVariables::variable_names[i]
                      << " relative error = " << max_err / max_chf * 100 << " %"
                      << std::endl;
            failed = 1;
            std::cout << " ----------------------------------------------------"
                         "-----------------------------------------"
                      << std::endl;
        }
    }

    if (failed == 0)
        std::cout << "MatterWeyl4 test passed..." << std::endl;
    else
        std::cout << "MatterWeyl4 test failed..." << std::endl;

    return failed;

#endif
}
