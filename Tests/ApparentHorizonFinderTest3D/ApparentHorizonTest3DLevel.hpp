/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef APPARENT_HORIZON_TEST3DLEVEL_HPP_
#define APPARENT_HORIZON_TEST3DLEVEL_HPP_

#include "AMRLevel.H"
#include "BoxLoops.hpp"
#include "CoarseAverage.H"
#include "DefaultLevelFactory.hpp"
#include "FourthOrderFillPatch.H"
#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"
#include "GRLevelData.hpp"
#include "GammaCalculator.hpp"
#include "InterpSource.hpp"
#include "KerrBH.hpp"
#include "LevelFluxRegister.H" //We don't actually use flux conservation but Chombo assumes we do
#include "LevelRK4.H"
#include "SetValue.hpp"

class ApparentHorizonTest3DLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<ApparentHorizonTest3DLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // initialize data
    virtual void initialData()
    {
        /*
        m_state_new.setVal(0);
        m_state_new.setVal(1, c_h11);
        m_state_new.setVal(1, c_h22);
        m_state_new.setVal(1, c_h33);
        const DisjointBoxLayout &level_domain = m_state_new.disjointBoxLayout();

        DataIterator dit = level_domain.dataIterator();
        double mass = m_p.mass;

        for (dit.begin(); dit.ok(); ++dit)
        {
            FArrayBox &fab = m_state_new[dit];
            const Box &b = level_domain[dit];

            const IntVect &smallEnd = b.smallEnd();
            const IntVect &bigEnd = b.bigEnd();

            const int xmin = smallEnd[0];
            const int ymin = smallEnd[1];
            const int zmin = smallEnd[2];

            const int xmax = bigEnd[0];
            const int ymax = bigEnd[1];
            const int zmax = bigEnd[2];

            // assigning values of 'chi' to a Schwarzschild BH in isotropic
            // coordinates
            for (int iz = zmin - 3; iz <= zmax + 3; ++iz)
                for (int iy = ymin - 3; iy <= ymax + 3; ++iy)
                    for (int ix = xmin - 3; ix <= xmax + 3; ++ix)
                    {
                        const double x = (ix + 0.5) * m_dx - m_p.center[0];
                        const double y = (iy + 0.5) * m_dx - m_p.center[1];
                        const double z = (iz + 0.5) * m_dx - m_p.center[2];
                        const double r = sqrt(x * x + y * y + z * z);
                        const IntVect iv(ix, iy, iz);
                        fab(iv, c_chi) = pow((1 + mass / (2 * r)), -4.0);
                    }
        }
        */

        // First set everything to zero then calculate initial data
        // Get the Kerr solution in the variables, then no need to calculate the
        // \tilde\Gamma^i numerically (not calculated in the Kerr ICs) as not
        // needed in AHFinder
        BoxLoops::loop(
            make_compute_pack(SetValue(0.), KerrBH(m_p.kerr_params, m_dx)),
            m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

        // Gamma's are not needed
        // fillAllGhosts();
        // BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
        // EXCLUDE_GHOST_CELLS);
    }

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time)
    {
    }

    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state)
    {
        tagging_criterion.setVal(0.);
    };
};

#endif /* APPARENT_HORIZON_TEST3DLEVEL_HPP_ */
