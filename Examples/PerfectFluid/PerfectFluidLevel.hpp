/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PERFECTFLUIDLEVEL_HPP_
#define PERFECTFLUIDLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"

// Problem specific includes
#include "EquationOfState.hpp"
#include "PerfectFluid.hpp"

//!  A class for the evolution of a scalar field, minimally coupled to gravity
/*! FIXME: Update for Fluid
     The class takes some initial data for a scalar field (variables phi and Pi)
     and evolves it using the CCZ4 equations. It is possible to specify an
   initial period of relaxation for the conformal factor chi, for non analytic
   initial conditions (for example, a general field configuration at a moment of
   time symmetry assuming conformal flatness). \sa MatterCCZ4(),
   ConstraintsMatter(), ScalarField(), RelaxationChi()
*/

class PerfectFluidLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<PerfectFluidLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // FIXME: Remove?
    // // Typedef for scalar field
    typedef PerfectFluid<EquationOfState> PerfectFluidWithEOS;

    //! Things to do at the end of the advance step, after RK4 calculation
    virtual void specificAdvance();

    //! Initialize data for the field and metric variables
    virtual void initialData();

    //! routines to do before outputing checkpoint file
    virtual void preCheckpointLevel();

    //! RHS routines used at each RK4 step
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time);

    //! Things to do in UpdateODE step, after soln + rhs update
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs, Real a_dt);

    //! Specify which variables to write at plot intervals
    virtual void specificWritePlotHeader(std::vector<int> &plot_states) const;

    //! Tell Chombo how to tag cells for regridding
    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state);
};

#endif /* PERFECTFLUIDLEVEL_HPP_ */
