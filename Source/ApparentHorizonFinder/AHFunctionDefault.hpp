/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHFUNCTIONDEFAULT_HPP_
#define _AHFUNCTIONDEFAULT_HPP_

#include "AlwaysInline.hpp"
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
        Tensor<2, double> g;
        FOR1(i) { g[i][i] = 1.; }
        return g;
    }

    // case of reduced Cartoon methods that have an extra metric component
#if GR_SPACEDIM != CH_SPACEDIM // hd - higher dimensions
    ALWAYS_INLINE double get_metric_hd() const { return 1.; }
#endif

    struct params // no params needed
    {
    };

    // some constructor with these arguments:
    // AHFunctionDefault(const AHData<int, double> &a_data,
    // const Tensor<1, double> &a_coords,
    // const Tensor<1, double> &a_coords_cartesian);

    // and some 'get'
    // double get(const AHGeometryData &geo_data, const AHDeriv &deriv,
    // const params &a_params) const;
};

#endif /* _AHFUNCTIONDEFAULT_HPP_ */
