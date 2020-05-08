/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BHAMR_HPP_
#define BHAMR_HPP_

#include "GRAMR.hpp"
#include "PunctureTracker.hpp"

/// A child of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy, and those of GRAMR
/**
 * This object inherits from GRAMR and adds tools required for BH spacetimes
 */
class BHAMR : public GRAMR
{
  public:
    PunctureTracker puncture_tracker;

    BHAMR() {}
};

#endif /* BHAMR_HPP_ */
