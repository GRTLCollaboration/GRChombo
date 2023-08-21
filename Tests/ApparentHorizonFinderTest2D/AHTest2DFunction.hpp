/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHTEST2DFUNCTION_HPP_
#define _AHTEST2DFUNCTION_HPP_

#include "AHFunctionDefault.hpp"

#include "AHDerivData.hpp"
#include "AHGeometryData.hpp"
#include "AHVarsData.hpp"

#include "UserVariables.hpp"

// to look for x == y
struct AHTest2DFunction : AHFunctionDefault
{
    double v;

    static int vars_min() { return c_V; }
    static int vars_max() { return c_V; }

    AHTest2DFunction(const AHVarsData<int, double> &a_data,
                     const Tensor<1, double> &a_coords,
                     const Tensor<1, double> &a_cart_coords)
    {
        v = a_data.vars.at(c_V);
    }

    double get(const AHGeometryData &geo_data, const AHDerivData &deriv,
               const params &a_params) const
    {
        return v;
    }
};

#endif /* _AHTEST2DFUNCTION_HPP_ */
