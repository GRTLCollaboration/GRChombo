/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOXPOINTERS_HPP_
#define BOXPOINTERS_HPP_

// Chombo includes
#include "FArrayBox.H"

// Other includes
#include "CellIndex.hpp"
#include "UserVariables.hpp"
#include <array>

// Chombo namespace
#include "UsingNamespace.H"

/// This class provides information about where a Chombo box lies in memory and
/// how it is laid out.
class BoxPointers
{
  private:
    const int *m_in_lo;
    const int *m_in_hi;

    const int *m_out_lo;
    const int *m_out_hi;

  public:
    std::vector<const double *> m_in_ptr;
    std::array<int, CH_SPACEDIM> m_in_stride; //!< Distance in memory between
                                              //! two values corresponding to
                                              //! adjacent coordinates in all
                                              //! directions.
    const int m_in_num_comps;

    std::vector<double *> m_out_ptr;
    std::array<int, CH_SPACEDIM> m_out_stride; //!< Distance in memory between
                                               //! two values corresponding to
                                               //! adjacent coordinates in all
                                               //! directions.
    const int m_out_num_comps;

    BoxPointers(const FArrayBox &in, FArrayBox &out)
        : m_in_lo(in.loVect()), m_in_hi(in.hiVect()), m_out_lo(out.loVect()),
          m_out_hi(out.hiVect()), m_in_num_comps(in.nComp()),
          m_out_num_comps(out.nComp())
    {
        m_in_ptr.resize(m_in_num_comps);
        m_out_ptr.resize(m_out_num_comps);
        // dataPtr in Chombo does CH_assert bound check
        // which we don't want to do in a loop
        for (int i = 0; i < m_in_num_comps; ++i)
            m_in_ptr[i] = in.dataPtr(i);
        // If the output FArrayBox doesn't have all components, just don't set
        // the pointers
        for (int i = 0; i < m_out_num_comps; ++i)
            m_out_ptr[i] = out.dataPtr(i);

        m_in_stride[0] = 1;
        m_in_stride[1] = m_in_hi[0] - m_in_lo[0] + 1;
#if CH_SPACEDIM >= 3
        m_in_stride[2] = (m_in_hi[1] - m_in_lo[1] + 1) * m_in_stride[1];
#endif

        m_out_stride[0] = 1;
        m_out_stride[1] = m_out_hi[0] - m_out_lo[0] + 1;
#if CH_SPACEDIM >= 3
        m_out_stride[2] = (m_out_hi[1] - m_out_lo[1] + 1) * m_out_stride[1];
#endif
    }

    CellIndexIn get_in_index(IntVect integer_coords) const
    {
        return (m_in_stride[2] * (integer_coords[2] - m_in_lo[2]) +
                m_in_stride[1] * (integer_coords[1] - m_in_lo[1]) +
                (integer_coords[0] - m_in_lo[0]));
    }

    CellIndexOut get_out_index(IntVect integer_coords) const
    {
        return (m_out_stride[2] * (integer_coords[2] - m_out_lo[2]) +
                m_out_stride[1] * (integer_coords[1] - m_out_lo[1]) +
                (integer_coords[0] - m_out_lo[0]));
    }
};

#endif /* BOXPOINTERS_HPP_ */
