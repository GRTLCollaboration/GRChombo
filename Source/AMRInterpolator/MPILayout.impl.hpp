/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MPILAYOUT_IMPL_HPP_
#define MPILAYOUT_IMPL_HPP_

inline MPILayout::MPILayout(int num_process)
    : m_num_process(num_process), m_counts(m_num_process, 0),
      m_displs(m_num_process, 0), m_dirty(false)
{
}

inline int MPILayout::count(int rank) const { return m_counts[rank]; }

inline int MPILayout::totalCount() const
{
    if (m_dirty)
        updateDirty();
    return m_total_count;
}

inline int MPILayout::displ(int rank) const
{
    if (m_dirty)
        updateDirty();
    return m_displs[rank];
}

inline void MPILayout::setCount(int rank, int count)
{
    CH_assert(rank < m_num_process && count >= 0);
    m_counts[rank] = count;
    m_dirty = true;
}

inline void MPILayout::incrementCount(int rank)
{
    CH_assert(rank < m_num_process);
    ++m_counts[rank];
    m_dirty = true;
}

inline void MPILayout::clearCounts()
{
    m_counts.assign(m_num_process, 0);
    m_dirty = true;
}

inline void MPILayout::updateDirty() const
{
    m_total_count = m_counts[0];
    for (int i = 1; i < m_num_process; ++i)
    {
        m_total_count += m_counts[i];
        m_displs[i] = m_displs[i - 1] + m_counts[i - 1];
    }
    m_dirty = false;
}

inline int *MPILayout::countsPtr() { return &m_counts[0]; }

inline int *MPILayout::displsPtr()
{
    if (m_dirty)
        updateDirty();
    return &m_displs[0];
}

#endif /* MPILAYOUT_IMPL_HPP_ */
