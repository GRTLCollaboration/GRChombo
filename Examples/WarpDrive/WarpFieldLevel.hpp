/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WARPFIELDLEVEL_HPP_
#define WARPFIELDLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "WarpField.hpp"

//!  A class for the evolution of a Warp field, minimally coupled to gravity
/*!
     The class takes some initial data for a Warp field
     and evolves it using the CCZ4 equations, with an ansatz for the matter.
     ConstraintsMatter(), WarpField()
*/
class WarpFieldLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<WarpFieldLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    double m_vs;
    double m_xs;

    //! Things to do at the end of the advance step, after RK4 calculation
    virtual void specificAdvance();

    //! Initialize data for the field and metric variables
    virtual void initialData();

    //! routines to do before outputing plot file
    virtual void prePlotLevel();

    //! RHS routines used at each RK4 step
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time);

    //! Things to do in UpdateODE step, after soln + rhs update
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs, Real a_dt);

    /// Things to do before tagging cells (i.e. filling ghosts)
    virtual void preTagCells() override;

    //! Tell Chombo how to tag cells for regridding
    virtual void computeDiagnosticsTaggingCriterion(
        FArrayBox &tagging_criterion,
        const FArrayBox &current_state_diagnostics);

    //! Tell Chombo how to tag cells for regridding
    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state);

    //! Things to do after a timestep eg extract GWs
    virtual void specificPostTimeStep();
};

#endif /* WARPFIELDLEVEL_HPP_ */
