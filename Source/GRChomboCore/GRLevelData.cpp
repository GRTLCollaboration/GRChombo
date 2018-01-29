/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "GRLevelData.hpp"
#include "FArrayBox.H"

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

// a_src and this must have the same box layout
void GRLevelData::plus(const GRLevelData &a_src, const double a_scale)
{
    DataIterator dit = m_disjointBoxLayout.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &fab = (*this)[dit];
        const FArrayBox &src_fab = a_src[dit];
        fab.plus(src_fab, a_scale);
    }
}
