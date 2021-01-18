/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETVALUE_HPP
#define SETVALUE_HPP

// Chombo includes
#include "Interval.H"

// Other includes
#include "Cell.hpp"
#include <limits>

// Chombo namespace
#include "UsingNamespace.H"

/// This compute class can be used to set the value of a cell
/** This compute class together with BoxLoops::loop(...) can be used to set the
 * values in a Box or in a LevelData to a constant value. Note: The
 * functionality is the same as GRLevelData.setVal but this compute class makes
 * it possible to bundle calculations up in a ComputeClassPack to avoid doing
 * several loops.
 */

class SetValue
{
    double m_value;
    Interval m_interval;

  public:
    SetValue(double a_value,
             Interval a_interval = Interval(0, std::numeric_limits<int>::max()))
        : m_value(a_value), m_interval(a_interval)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        int start_var = m_interval.begin();
        int end_var =
            std::min(m_interval.end(), current_cell.get_num_out_vars() - 1);
        for (int i = start_var; i <= end_var; ++i)
        {
            current_cell.store_vars(m_value, i);
        }
    }
};

#endif /* SETVALUE_HPP */
