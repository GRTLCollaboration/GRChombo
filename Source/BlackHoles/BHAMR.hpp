/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BHAMR_HPP_
#define BHAMR_HPP_

#include "GRAMR.hpp"
#if CH_SPACEDIM == 3
#include "PunctureTracker.hpp"
#endif

#ifdef USE_AHFINDER
#include "AHFinder.hpp"
#endif

/// A child of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy, and those of GRAMR
/**
 * This object inherits from GRAMR and adds tools required for BH spacetimes
 */
class BHAMR : public GRAMR
{
  public:
#if CH_SPACEDIM == 3
    PunctureTracker m_puncture_tracker;
#endif

#ifdef USE_AHFINDER
    AHFinder<> m_ah_finder;
#endif

    BHAMR() {}

    void set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator) override
    {
        GRAMR::set_interpolator(a_interpolator);
#if CH_SPACEDIM == 3
        m_puncture_tracker.set_interpolator(a_interpolator);
#endif
#ifdef USE_AHFINDER
        m_ah_finder.set_interpolator(a_interpolator);
#endif
    }
};

#endif /* BHAMR_HPP_ */
