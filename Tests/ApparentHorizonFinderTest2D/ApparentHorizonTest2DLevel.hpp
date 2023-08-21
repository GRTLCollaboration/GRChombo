/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef APPARENT_HORIZON_TEST2DLEVEL_HPP_
#define APPARENT_HORIZON_TEST2DLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "GRLevelData.hpp"

class ApparentHorizonTest2DLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<ApparentHorizonTest2DLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // initialize data
    virtual void initialData()
    {

        const DisjointBoxLayout &level_domain = m_state_new.disjointBoxLayout();

        DataIterator dit = level_domain.dataIterator();

        for (dit.begin(); dit.ok(); ++dit)
        {
            FArrayBox &fab = m_state_new[dit];
            const Box &b = level_domain[dit];

            const IntVect &smallEnd = b.smallEnd();
            const IntVect &bigEnd = b.bigEnd();

            const int xmin = smallEnd[0];
            const int ymin = smallEnd[1];

            const int xmax = bigEnd[0];
            const int ymax = bigEnd[1];

            // assigning values of 'chi' to a Schwarzschild BH in isotropic
            // coordinates
            for (int iy = ymin - 3; iy <= ymax + 3; ++iy)
                for (int ix = xmin - 3; ix <= xmax + 3; ++ix)
                {
                    const double x = (ix + 0.5) * m_dx - m_p.center[0];
                    const double y = (iy + 0.5) * m_dx - m_p.center[1];
                    const IntVect iv(ix, iy);
                    static const double perturb = 1;
                    const double y_perturbed =
                        abs(y) + perturb * sin(2. * M_PI * x / m_p.L);
                    // this is basically a decaying wave (sin(y)/y) slightly
                    // deformed in the 'x' direction
                    fab(iv, c_V) = sin(y_perturbed) / y_perturbed;
                }
        }
    }

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time)
    {
    }

    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state){};
};

#endif /* APPARENT_HORIZON_TEST2DLEVEL_HPP_ */
