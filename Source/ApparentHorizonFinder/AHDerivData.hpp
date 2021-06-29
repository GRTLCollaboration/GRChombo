/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHDERIV_HPP_
#define _AHDERIV_HPP_

#include "DimensionDefinitions.hpp" // make sure GR_SPACEDIM exists

#include <cstring> // memset

#define DWIDTH 5
#define DDWIDTH 6

//! class to store 1st and 2nd derivatives of 'F'
struct AHDerivData
{
    double duF;
    double duduF;
#if CH_SPACEDIM == 3
    double dvF;
    double dvdvF;
    double dudvF;
#endif

    int du_stencil_start;
    int dudu_stencil_start;
#if CH_SPACEDIM == 3
    int dv_stencil_start;
    int dvdv_stencil_start;
#endif

    double du_weights[DWIDTH];
    double dudu_weights[DDWIDTH];
#if CH_SPACEDIM == 3
    double dv_weights[DWIDTH];
    double dvdv_weights[DDWIDTH];
#endif

    AHDerivData()
    {
        // force all (double) elements of AHDerivData to be 0
        memset(this, 0, sizeof(AHDerivData));
    }
};

#endif /* _AHDERIV_HPP_ */
