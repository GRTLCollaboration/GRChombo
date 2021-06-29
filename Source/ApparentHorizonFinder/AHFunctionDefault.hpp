/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHFUNCTIONDEFAULT_HPP_
#define _AHFUNCTIONDEFAULT_HPP_

#include "AHDerivData.hpp"
#include "AHGeometryData.hpp"
#include "AHVarsData.hpp"
#include "AlwaysInline.hpp"
#include "GRParmParse.hpp"
#include "Tensor.hpp"

/////////////////////////////////////////////////////////
// Predefined optimization functions
/////////////////////////////////////////////////////////

struct AHFunctionDefault
{
    // set the minimum/maximum variable needed (and same for 1st and 2nd
    // derivatives)
    static ALWAYS_INLINE int vars_min() { return 0; }
    static ALWAYS_INLINE int vars_max() { return -1; }
    static ALWAYS_INLINE int d1_vars_min() { return 0; }
    static ALWAYS_INLINE int d1_vars_max() { return -1; }
    static ALWAYS_INLINE int d2_vars_min() { return 0; }
    static ALWAYS_INLINE int d2_vars_max() { return -1; }

    // set how many and what variables to print for the geometry when
    // 'AH_print_geometry_data' is on
    // (useful to print more complex variables as the non-conformal metric or
    // the extrinsic curvature, which are not elementary variables in
    // BSSN/CCZ4)
    static ALWAYS_INLINE int num_write_vars() { return 0; }
    static ALWAYS_INLINE void write_headers(std::string *vec) {}
    ALWAYS_INLINE void write_vars(double *vec) {}

    // metric to use when calculating the spin and area of the AH
    // (override with spacetime metric)
    ALWAYS_INLINE const Tensor<2, double> get_metric() const
    {
        // cartesian flat metric
        Tensor<2, double> g = {0.};
        FOR(i) { g[i][i] = 1.; }
        return g;
    }
    ALWAYS_INLINE const Tensor<2, double> get_inverse_metric() const
    {
        // cartesian flat metric
        return get_metric();
    }
    ALWAYS_INLINE const Tensor<2, double> get_extrinsic_curvature() const
    {
        return {0.};
    }

    // case of reduced Cartoon methods that have an extra metric component
#if GR_SPACEDIM != CH_SPACEDIM // hd - higher dimensions
    ALWAYS_INLINE double get_metric_hd() const { return 1.; }
#endif

    struct params // no params needed
    {
        void read_params(GRParmParse &pp) {}
    };

    // not defined by default
    Tensor<1, double>
    get_level_function_derivative(const AHGeometryData &geo_data,
                                  const AHDerivData &deriv) const
    {
        return {0.};
    }
    Tensor<1, double> get_spatial_normal_U(const Tensor<1, double> &s_L) const
    {
        return {0.};
    }
    // not defined by default
    Tensor<2, double> get_level_function_2nd_covariant_derivative(
        const AHGeometryData &geo_data, const AHDerivData &deriv,
        const Tensor<1, double> &s_L) const
    {
        return {0.};
    }

    // WHAT TO ADD TO YOUR OWN FUNCTIONS:
    // some constructor with these arguments:
    // AHFunctionDefault(const AHVarsData<int, double> &a_data,
    // const Tensor<1, double> &a_coords,
    // const Tensor<1, double> &a_coords_cartesian);

    // and some 'get'
    // double get(const AHGeometryData &geo_data, const AHDerivData &deriv,
    // const params &a_params) const;
};

#endif /* _AHFUNCTIONDEFAULT_HPP_ */
