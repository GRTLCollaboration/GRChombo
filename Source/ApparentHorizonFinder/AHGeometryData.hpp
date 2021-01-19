/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHGEOMETRYDATA_HPP_
#define _AHGEOMETRYDATA_HPP_

#include "DimensionDefinitions.hpp" // make sure GR_SPACEDIM exists
#include "Tensor.hpp"

#include <cstring> // memset

// The d prefix refers to partial derivatives wrt Cartesian coordinates

struct AHGeometryData
{

    // jacobian
    Tensor<1, double> du;
#if CH_SPACEDIM == 3
    Tensor<1, double> dv;
#endif
    Tensor<1, double> df;

    // hessian
    Tensor<2, double> ddu;
#if CH_SPACEDIM == 3
    Tensor<2, double> ddv;
#endif
    Tensor<2, double> ddf;

    // inverse jacobian
    Tensor<1, double> dxdu;
#if CH_SPACEDIM == 3
    Tensor<1, double> dxdv;
#endif
    Tensor<1, double> dxdf;

    // inverse hessian
    // never needed

    AHGeometryData()
    {
        // TF: no need to set to 0
        // force all (double) elements of AHGeometryData to be 0
        // memset(this, 0, sizeof(AHGeometryData));
    }
};

#endif /* _AHGEOMETRYDATA_HPP_ */
