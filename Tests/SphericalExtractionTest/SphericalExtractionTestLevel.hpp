/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPHERICALEXTRACTIONTESTLEVEL_HPP_
#define SPHERICALEXTRACTIONTESTLEVEL_HPP_

#include "GRAMRLevel.hpp"

class SphericalExtractionTestLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<SphericalExtractionTestLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // initialize data
    virtual void initialData() { m_state_new.setVal(42.); }

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time)
    {
    }

    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state){};
};

#endif /* SURFACEEXTRACTIONTESTLEVEL_HPP_ */
