/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTERPOLATIONLAYOUT_HPP_
#define INTERPOLATIONLAYOUT_HPP_

#include <vector>

class InterpolationLayout
{
  private:
    template <typename InterpAlgo> friend class AMRInterpolator;

    std::vector<int> rank;
    std::vector<int> level_idx;
    std::vector<int> box_idx;

    InterpolationLayout(int num_points)
        : rank(num_points, -1), level_idx(num_points, -1),
          box_idx(num_points, -1)
    {
    }
};

#endif /* INTERPOLATIONLAYOUT_HPP_ */
