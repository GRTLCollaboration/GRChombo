/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "FArrayBox.H"

// Our includes
#include "GRLevelData.hpp"

// Chombo namespace
#include "UsingNamespace.H"

GRLevelData::GRLevelData() : LevelData<FArrayBox>() {}

void GRLevelData::setVal(const double a_val)
{
    DataIterator dit = m_disjointBoxLayout.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &fab = (*this)[dit];
        fab.setVal(a_val);
    }
}

void GRLevelData::setVal(const double a_val, const int a_comp)
{
    DataIterator dit = m_disjointBoxLayout.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &fab = (*this)[dit];
        fab.setVal(a_val, a_comp);
    }
}

void GRLevelData::setVal(const double a_val, const Interval a_comps)
{
    DataIterator dit = m_disjointBoxLayout.dataIterator();
    // Want component loop inside so unfortunately we have to duplicate the
    // outer loop
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &fab = (*this)[dit];
        for (int i = a_comps.begin(); i <= a_comps.end(); ++i)
        {
            fab.setVal(a_val, i);
        }
    }
}

void GRLevelData::plus(const GRLevelData &a_src, const double a_scale,
                       const DisjointBoxLayout &a_disjoint_box_layout)
{
    CH_TIME("GRLevelData::plus");
    DataIterator dit = a_disjoint_box_layout.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &fab = (*this)[dit];
        const FArrayBox &src_fab = a_src[dit];
        Real *fab_data = fab.dataPtr();
        const Real *const src_fab_data = src_fab.dataPtr();

        const Box &this_box = fab.box();
        const Box &src_box = src_fab.box();
        const Box loop_box = a_disjoint_box_layout[dit];
        const int *this_lo = this_box.loVect();
        const int *src_lo = src_box.loVect();
        const int *loop_lo = loop_box.loVect();
        const int *loop_hi = loop_box.hiVect();
        const IntVect &this_size = this_box.size();
        const IntVect &src_size = src_box.size();
        const int num_comps = fab.nComp();

        CH_assert(num_comps == src_fab.nComp());
#ifdef _OPENMP
#pragma omp parallel for default(shared) collapse(CH_SPACEDIM)
#endif
        for (int icomp = 0; icomp < num_comps; ++icomp)
#if CH_SPACEDIM == 3
            for (int iz = loop_lo[2]; iz <= loop_hi[2]; ++iz)
#endif
                for (int iy = loop_lo[1]; iy <= loop_hi[1]; ++iy)
                {
#pragma omp simd
                    for (int ix = loop_lo[0]; ix <= loop_hi[0]; ++ix)
                    {
                        const int this_index =
                            ix - this_lo[0] + (iy - this_lo[1]) * this_size[0] +
#if CH_SPACEDIM == 3
                            (iz - this_lo[2]) * this_size[0] * this_size[1] +
#endif
                            icomp * this_size[0] * this_size[1]
#if CH_SPACEDIM == 3
                                * this_size[2]
#endif
                            ;
                        const int src_index =
                            ix - src_lo[0] + (iy - src_lo[1]) * src_size[0] +
#if CH_SPACEDIM == 3
                            (iz - src_lo[2]) * src_size[0] * src_size[1] +
#endif
                            icomp * src_size[0] * src_size[1]
#if CH_SPACEDIM == 3
                                * src_size[2]
#endif
                            ;
                        fab_data[this_index] +=
                            a_scale * src_fab_data[src_index];
                    }
                }
    }
}

// old plus function
// void GRLevelData::plus(const GRLevelData &a_src, const double a_scale,
//                        const DisjointBoxLayout &a_disjoint_box_layout)
// {
//    CH_TIME("GRLevelData::plus");
//    DataIterator dit = a_disjoint_box_layout.dataIterator();
//     for (dit.begin(); dit.ok(); ++dit)
//     {
//         FArrayBox &fab = (*this)[dit];
//         const FArrayBox &src_fab = a_src[dit];
//         fab.plus(src_fab, a_scale);
//     }
// }
