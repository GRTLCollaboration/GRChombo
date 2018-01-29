/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BINARYBHLEVEL_HPP_
#define BINARYBHLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"

class BinaryBHLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<BinaryBHLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    virtual void specificAdvance() override;

    /// Initial data calculation
    virtual void initialData() override;

    /// Any actions that should happen just before checkpointing
    virtual void preCheckpointLevel() override;

    /// Calculation of the right hand side for the time stepping
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) override;

    /// Things to do after dt*rhs has been added to the solution
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   Real a_dt) override;

    /// Identify and tag the cells that need higher resolution
    virtual void
    computeTaggingCriterion(FArrayBox &tagging_criterion,
                            const FArrayBox &current_state) override;
};

#endif /* BINARYBHLEVEL_HPP_ */
