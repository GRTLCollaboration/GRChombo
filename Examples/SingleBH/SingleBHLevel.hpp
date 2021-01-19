/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SINGLEBHLEVEL_HPP_
#define SINGLEBHLEVEL_HPP_

#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"

class SingleBHLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<SingleBHLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    BHAMR &m_bh_amr = dynamic_cast<BHAMR &>(m_gr_amr);

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    virtual void specificAdvance() override;

    /// Initial data calculation
    virtual void initialData() override;

    /// Any actions that should happen just before checkpointing
    virtual void prePlotLevel() override;

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

    // to do post each time step on every level
    virtual void specificPostTimeStep() override;
};

#endif /* SINGLEBHLEVEL_HPP_ */
