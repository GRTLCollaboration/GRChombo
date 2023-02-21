/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INSITUTESTLEVEL_HPP_
#define INSITUTESTLEVEL_HPP_

#include "BoxLoops.hpp"
#include "GRAMRLevel.hpp"
#include "Polynomial.hpp"
#include "SetValue.hpp"

class InsituTestLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<InsituTestLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // initialize data
    virtual void initialData()
    {
        BoxLoops::loop(Polynomial(m_p.center, m_dx), m_state_new, m_state_new,
                       FILL_GHOST_CELLS);
    }

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time)
    {
    }

    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state){};
};

#endif /* INSITUTESTLEVEL_HPP_ */
