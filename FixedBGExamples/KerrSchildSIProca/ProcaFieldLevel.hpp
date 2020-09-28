/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PROCAFIELDLEVEL_HPP_
#define PROCAFIELDLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "FixedBGProcaField.hpp"
#include "Potential.hpp"

//!  A class for the evolution of a proca field, minimally coupled to gravity
/*!
     The class takes some initial data for a proca field (variables phi and Pi)
     and evolves it on a fixed BG
*/
class ProcaFieldLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<ProcaFieldLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // Typedef for proca field
    typedef FixedBGProcaField<Potential> ProcaField;

    //! Things to do at the end of the advance step, after RK4 calculation
    virtual void specificAdvance();

    //! Initialize data for the field and metric variables
    virtual void initialData();

    //! routines to do before outputing plot file
    virtual void prePlotLevel();

    //! RHS routines used at each RK4 step
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time);

    //! To do after each timestep
    virtual void specificPostTimeStep();

    //! Things to do in UpdateODE step, after soln + rhs update
    //    virtual void specificUpdateODE(GRLevelData &a_soln,
    //                                   const GRLevelData &a_rhs, Real a_dt);

    //! Tell Chombo how to tag cells for regridding
    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state);
};

#endif /* PROCAFIELDLEVEL_HPP_ */
