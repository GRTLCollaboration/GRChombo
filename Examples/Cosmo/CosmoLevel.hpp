/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COSMOLEVEL_HPP_
#define COSMOLEVEL_HPP_

#include "CosmoAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "Potential.hpp"
#include "ScalarField.hpp"

//!  A class for the evolution of a scalar field in cosmological spacetime
/*!
     The class takes initial data for a scalar field (variables phi and Pi)
     and evolves it using the CCZ4 equations. An example of specifying initial
   scalar field data is provided as an example. The given initial data example
   is obtained by analytically solving the constraint using a sinusoidal scalar
   field profile with a quadratic potential. See initial scalar field details in
   InitialData(), InitialScalarData.hpp and Potential.hpp
*/
class CosmoLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<CosmoLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    CosmoAMR &m_cosmo_amr = dynamic_cast<CosmoAMR &>(m_gr_amr);

    // Typedef for scalar field
    typedef ScalarField<Potential> ScalarFieldWithPotential;

    //! Things to do at the end of the advance step, after RK4 calculation
    virtual void specificAdvance() override;

    //! Initialize data for the field and metric variables
    virtual void initialData() override;

    //! Recalculate and set K after restart
    virtual void postRestart() override;

#ifdef CH_USE_HDF5
    //! routines to do before outputting plot file
    virtual void prePlotLevel() override;
#endif

    //! RHS routines used at each RK4 step
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) override;

    //! Things to do in UpdateODE step, after soln + rhs update
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   Real a_dt) override;

    /// Things to do before tagging cells (i.e. filling ghosts)
    virtual void preTagCells() override;

    //! Tell Chombo how to tag cells for regridding
    virtual void computeTaggingCriterion(
        FArrayBox &tagging_criterion, const FArrayBox &current_state,
        const FArrayBox &current_state_diagnostics) override;

    //! to do post each time step on every level
    virtual void specificPostTimeStep() override;
};

#endif /* COSMOLEVEL_HPP_ */
