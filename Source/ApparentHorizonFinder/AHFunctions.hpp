/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHFUNCTIONS_HPP_
#define _AHFUNCTIONS_HPP_

// Chombo includes
#include "CH_Timer.H"

// Other includes
#include "DimensionDefinitions.hpp" // make sure GR_SPACEDIM exists
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include <cmath> // sqrt

#include "AHFunctionDefault.hpp"

// Chombo namespace
#include "UsingNamespace.H"

/////////////////////////////////////////////////////////
// Predefined optimization functions
/////////////////////////////////////////////////////////

struct ExpansionFunction : AHFunctionDefault
{
    //////////////////////////////////
    // needed to calculate expansion:

    // induced metric
    Tensor<2, double> g;
    Tensor<2, double> g_UU;
    Tensor<3, double> dg;

    // extrinsic curvature
    Tensor<2, double> K;
    double trK;

    // target coordinate (e.g. radius)
    double f;

    // case of reduced Cartoon methods that have an extra metric component and
    // require the coordinates to calculate the expansion
#if GR_SPACEDIM != CH_SPACEDIM
    // hd - higher dimensions
    Tensor<1, double> coords; // cartesian coordinates
    double g_hd;
    Tensor<1, double> dg_hd;

    double Kww; // not needed to calculate expansion, but may be needed for ut
    ALWAYS_INLINE double get_metric_hd() const { return g_hd; }
#endif
    //////////////////////////////////

    // only require variables up to Aij (chi, hij, K, Aij)
    // this assumes c_Theta comes after the last of Aij
    // (and like this the code is generic for 2D and 3D)
    static ALWAYS_INLINE int vars_min() { return c_chi; }
    static ALWAYS_INLINE int vars_max() { return c_Theta - 1; }
    // Derivatives required only for the metric component (chi and hij)
    // this assumes c_K comes after the last hij
    // (and like this the code is generic for 2D and 3D)
    static ALWAYS_INLINE int d1_vars_min() { return c_chi; }
    static ALWAYS_INLINE int d1_vars_max() { return c_K - 1; }

    ALWAYS_INLINE const Tensor<2, double> get_metric() const { return g; }
    ALWAYS_INLINE const Tensor<2, double> get_inverse_metric() const
    {
        return g_UU;
    }
    ALWAYS_INLINE const Tensor<2, double> get_extrinsic_curvature() const
    {
        return K;
    }

    static int num_write_vars()
    {
        int num_components = CH_SPACEDIM * (CH_SPACEDIM + 1) / 2 * 2;
#if GR_SPACEDIM != CH_SPACEDIM
        num_components += 2; // hww and Kww
#endif
        return num_components;
    }
    static void write_headers(std::string *vec)
    {
        int el = 0;
        for (int i = 0; i < CH_SPACEDIM; ++i)
            for (int j = i; j < CH_SPACEDIM; ++j)
                vec[el++] = "g" + std::to_string(i + 1) + std::to_string(j + 1);
#if GR_SPACEDIM != CH_SPACEDIM
        vec[el++] = "gww";
#endif

        for (int i = 0; i < CH_SPACEDIM; ++i)
            for (int j = i; j < CH_SPACEDIM; ++j)
                vec[el++] = "K" + std::to_string(i + 1) + std::to_string(j + 1);
#if GR_SPACEDIM != CH_SPACEDIM
        vec[el++] = "Kww";
#endif
    }
    void write_vars(double *vec)
    {
        int el = 0;
        for (int i = 0; i < CH_SPACEDIM; ++i)
            for (int j = i; j < CH_SPACEDIM; ++j)
                vec[el++] = g[i][j];
#if GR_SPACEDIM != CH_SPACEDIM
        vec[el++] = g_hd;
#endif
        for (int i = 0; i < CH_SPACEDIM; ++i)
            for (int j = i; j < CH_SPACEDIM; ++j)
                vec[el++] = K[i][j];
#if GR_SPACEDIM != CH_SPACEDIM
        vec[el++] = Kww;
#endif
    }

    ExpansionFunction(const AHVarsData<int, double> &a_data,
                      const Tensor<1, double> &a_coords,
                      const Tensor<1, double> &a_coords_cartesian)
    {
        CH_TIME("ExpansionFunction::calculate_data");

        f = a_coords[CH_SPACEDIM - 1];

        // * ---------------------------
        // *       GR-RELATED DATA
        // * ---------------------------

        int h[CH_SPACEDIM][CH_SPACEDIM];
        int A[CH_SPACEDIM][CH_SPACEDIM];

        int comp_h = c_h11;
        int comp_A = c_A11;
        for (int i = 0; i < CH_SPACEDIM; ++i)
        {
            for (int j = i; j < CH_SPACEDIM; ++j)
            {
                h[i][j] = comp_h;
                A[i][j] = comp_A;
                if (i != j)
                {
                    h[j][i] = comp_h;
                    A[j][i] = comp_A;
                }
                ++comp_h;
                ++comp_A;
            }
        }

        const double chi = a_data.vars.at(c_chi);
        trK = a_data.vars.at(c_K);

        // INVERSE METRIC
        Tensor<2, double, CH_SPACEDIM> h_DD;
        FOR(i, j) { h_DD[i][j] = a_data.vars.at(h[i][j]); }

        Tensor<2, double, CH_SPACEDIM> h_UU =
            TensorAlgebra::compute_inverse_sym(h_DD);

        // Reconstructing ADM variables
        Tensor<1, double, CH_SPACEDIM> dchi;

        FOR(i) { dchi[i] = a_data.d1.at(c_chi)[i]; }

        for (int i = 0; i < CH_SPACEDIM; ++i)
        {
            for (int j = i; j < CH_SPACEDIM; ++j)
            {
                {
                    const double gij = h_DD[i][j] / chi;
                    g[i][j] = gij;
                    g[j][i] = gij;
                }
                {
                    const double g_UUij = h_UU[i][j] * chi;
                    g_UU[i][j] = g_UUij;
                    g_UU[j][i] = g_UUij;
                }

                {
                    const double Aij = a_data.vars.at(A[i][j]);
                    const double Kij = Aij / chi + trK * g[i][j] / GR_SPACEDIM;
                    K[i][j] = Kij;
                    K[j][i] = Kij;
                }

                for (int k = 0; k < CH_SPACEDIM; ++k)
                {
                    {
                        const double dhij = a_data.d1.at(h[i][j])[k];
                        const double dgijk =
                            (dhij - (h_DD[i][j] * dchi[k]) / chi) / chi;
                        dg[i][j][k] = dgijk;
                        dg[j][i][k] = dgijk;
                    }
                }
            }
        }

        // part for higher dimensions that use Cartoon Method
        // Uli's paper does it for 3+1D (from 5D) - arxiv 1808.05834
#if GR_SPACEDIM != CH_SPACEDIM
        int comp_hww =
            c_K - 1; // c_K-1 expected to be c_hww (is c_hww direct better?)
        int comp_Aww = c_Aww; // is (c_Theta-1) more general?
        Tensor<1, double, CH_SPACEDIM> dhww;

        coords = a_coords_cartesian;

        double hww = a_data.vars.at(comp_hww);
        FOR(a) { dhww[a] = a_data.d1.at(comp_hww)[a]; }
        double Aww = a_data.vars.at(comp_Aww);

        g_hd = hww / chi;
        FOR(a) { dg_hd[a] = (dhww[a] - (hww * dchi[a]) / chi) / chi; }

        Kww = Aww / chi + trK * g_hd / GR_SPACEDIM;
#endif
    }

    struct params
    {
        double expansion_radius_power;

        void read_params(GRParmParse &pp)
        {
            pp.load("AH_expansion_radius_power", expansion_radius_power, 1.);
        }
    };

    double get(const AHGeometryData &geo_data, const AHDerivData &deriv,
               const params &a_params) const
    {
        Tensor<1, double> s_L = get_level_function_derivative(geo_data, deriv);
        Tensor<2, double> Ds =
            get_level_function_2nd_covariant_derivative(geo_data, deriv, s_L);
        double norm_s;
        Tensor<1, double> S_U;
        get_spatial_normal_U_and_norm(S_U, norm_s, s_L);

        // calculate D_i S^i and S^i S^j K_ij
        double DiSi = 0.;
        double Kij_dot_Si_Sj = 0.;
        FOR(a, b)
        {
            DiSi += (g_UU[a][b] - S_U[a] * S_U[b]) * Ds[a][b] / norm_s;
            Kij_dot_Si_Sj += S_U[a] * S_U[b] * K[a][b];
        }

        // Calculation of expansion - as in (6.7.9 / 6.7.13) of Alcubierre
        double expansion = DiSi - trK + Kij_dot_Si_Sj;

        // part from extra dimensions in the case of Cartoon methods
#if GR_SPACEDIM != CH_SPACEDIM
        expansion += (GR_SPACEDIM - CH_SPACEDIM) * S_U[CH_SPACEDIM - 1] /
                     coords[CH_SPACEDIM - 1];
        FOR(a)
        {
            expansion +=
                (GR_SPACEDIM - CH_SPACEDIM) * 0.5 * dg_hd[a] / g_hd * S_U[a];
        }
#endif

        // using "r * Expansion" significantly improves the convergence
        // (making a Schw. BH converge for any radius >~ 0.5*r_AH instead of
        // only up to ~ 3 * r_AH as it happens just with the expansion)
        // (see arXiv:gr-qc/0512169, around Fig 3.4)
        return expansion * pow(f, a_params.expansion_radius_power);
    }

    // extra stuff:
    Tensor<1, double>
    get_level_function_derivative(const AHGeometryData &geo_data,
                                  const AHDerivData &deriv) const
    {
        // calculate D_a L of 6.7.12 of Alcubierre for some level function
        // L picking L = f - F(u,v)
        // D_a L = d_a f - dF/du * du/dx^a - dF/dv * dv/dx^a

        Tensor<1, double> s_L = {0.}; // not normalized, just D_a L
        FOR(a)
        {
            s_L[a] = D_TERM(geo_data.df[a], -(deriv.duF * geo_data.du[a]),
                            -(deriv.dvF * geo_data.dv[a]));
        }

        return s_L;
    }
    Tensor<1, double> get_spatial_normal_U(const Tensor<1, double> &s_L) const
    {
        Tensor<1, double> S_U;
        double norm_s;
        get_spatial_normal_U_and_norm(S_U, norm_s, s_L);
        return S_U;
    }
    void get_spatial_normal_U_and_norm(Tensor<1, double> &S_U, double &norm_s,
                                       const Tensor<1, double> &s_L) const
    {
        // calculate S_U = the real 's' of Alcubierre

        // norm of s_L = | D_a L| (the 'u' in 6.7.12 of Alcubierre)
        norm_s = 0.0;
        FOR(a, b) { norm_s += g_UU[a][b] * s_L[a] * s_L[b]; }
        norm_s = sqrt(norm_s);

        // raise s_L
        Tensor<1, double> s_U = {0.};
        FOR(a, b) { s_U[a] += g_UU[a][b] * s_L[b]; }

        FOR(a) { S_U[a] = s_U[a] / norm_s; }
    }

    Tensor<2, double> get_level_function_2nd_covariant_derivative(
        const AHGeometryData &geo_data, const AHDerivData &deriv,
        const Tensor<1, double> &s_L) const
    {
        // calculates D_a D_b L, required for 6.7.13 of Alcubierre
        // for the level function L = f - F(u,v)

        Tensor<2, double> ds = {0.};
        FOR(a, b)
        {
            ds[a][b] =
                D_TERM(geo_data.ddf[a][b],
                       -(deriv.duF * geo_data.ddu[a][b]) -
                           (deriv.duduF * geo_data.du[a] * geo_data.du[b]),
                       -(deriv.dvF * geo_data.ddv[a][b]) -
                           (deriv.dvdvF * geo_data.dv[a] * geo_data.dv[b]) -
                           (deriv.dudvF * (geo_data.du[a] * geo_data.dv[b] +
                                           geo_data.du[b] * geo_data.dv[a])));
        }

        // calculate Christoffels on this point (u,v,f)
        Tensor<3, double> chris = {0.};
        FOR(a, b, c, d)
        {
            chris[a][b][c] +=
                0.5 * g_UU[a][d] * (dg[b][d][c] + dg[c][d][b] - dg[b][c][d]);
        }

        // covariant derivatrive of s_a to use for DS
        Tensor<2, double> Ds = {0.};
        FOR(a, b)
        {
            Ds[a][b] = ds[a][b];
            FOR(c) { Ds[a][b] -= chris[c][a][b] * s_L[c]; }
        }

        return Ds;
    }
};

// to look for chi contours
struct ChiContourFunction : AHFunctionDefault
{
    double chi;

    static int vars_min() { return c_chi; }
    static int vars_max() { return c_chi; }

    ChiContourFunction(const AHVarsData<int, double> &a_data,
                       const Tensor<1, double> &a_coords,
                       const Tensor<1, double> &a_coords_cartesian)
    {
        chi = a_data.vars.at(c_chi);
    }

    struct params
    {
        double look_for_chi_contour;

        void read_params(GRParmParse &pp)
        {
            pp.load("AH_look_for_chi_contour", look_for_chi_contour);
            CH_assert(look_for_chi_contour > 0.);
        }
    };

    double get(const AHGeometryData &geo_data, const AHDerivData &deriv,
               const params &a_params) const
    {
        return chi - a_params.look_for_chi_contour;
    }
};

#endif /* _AHFUNCTIONS_HPP_ */
