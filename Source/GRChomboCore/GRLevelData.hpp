/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRLEVELDATA_HPP_
#define GRLEVELDATA_HPP_

// Chombo includes
#include "FArrayBox.H"
#include "LevelData.H"

// Chombo namespace
#include "UsingNamespace.H"

class GRLevelData : public LevelData<FArrayBox>
{
  public:
    GRLevelData();

    void setVal(const double a_val);

    void setVal(const double a_val, const int a_comp);

    void setVal(const double a_val, const Interval a_comps);

    // loop only goes over a_disjoint_box_layout
    void plus(const GRLevelData &a_src, const double a_scale,
              const DisjointBoxLayout &a_disjoint_box_layout);
};

#endif /* GRLEVELDATA_HPP_ */
