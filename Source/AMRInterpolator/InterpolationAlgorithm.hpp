/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTERPOLATIONALGORITHM_HPP_
#define INTERPOLATIONALGORITHM_HPP_

// Chombo includes
#include "FArrayBox.H"

// Other includes
#include <array>

// Chombo namespace
#include "UsingNamespace.H"

class InterpolationAlgorithm
{
  public:
    virtual ~InterpolationAlgorithm() = 0;
};

class NearestNeighbour : public InterpolationAlgorithm
{
  public:
    static inline double
    interpPoint(const std::array<double, CH_SPACEDIM> &gridCoord,
                const FArrayBox &fab, int comps, const IntVect &nearest)
    {
        return fab.get(nearest, comps);
    }
};

#endif /* INTERPOLATIONALGORITHM_HPP_ */
