/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTERPSOURCE_H_
#define INTERPSOURCE_H_

#include <array>
#include "LevelData.H"

// Abstrace base class to get the FABs out of an AMRLevel
class InterpSource
{
  public:
    virtual const LevelData<FArrayBox> &getLevelData() const = 0;
    virtual bool
    contains(const std::array<double, CH_SPACEDIM> &point) const = 0;
    virtual void fillAllGhosts() = 0;
};

#endif /* INTERPSOURCE_H_ */
