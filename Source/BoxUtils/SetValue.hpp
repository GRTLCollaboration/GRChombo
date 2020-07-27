/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETVALUE_HPP
#define SETVALUE_HPP

#include "Cell.hpp"
#include "Interval.H"

/// This compute class can be used to set the value of a cell
/** This compute class together with BoxLoops::loop(...) can be used to set the
 * values in a Box or in a LevelData to a constant value. Note: The
 * functionality is the same as GRLevelData.setVal but this compute class makes
 * it possible to bundle calculations up in a ComputeClassPack to avoid doing
 * several loops.
 */


template <typename...> class SetValueTemplate;

template <typename T> class SetValueTemplate<T>
{
    double m_value;
    T m_interval;

  public:
    SetValueTemplate(double a_value, T a_interval)
        : m_value(a_value), m_interval(a_interval)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        int start_var = m_interval.begin();
        int end_var = m_interval.end();
        for (int i = start_var; i <= end_var; ++i)
        {
            current_cell.store_vars(m_value, i);
        }
    }
};

template <> class SetValueTemplate<>
{
    double m_value;

  public:
    SetValueTemplate(double a_value)
        : m_value(a_value)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        for (int i = 0; i < current_cell.get_num_out_vars(); ++i)
        {
            current_cell.store_vars(m_value, i);
        }
    }
};

using SetValue = SetValueTemplate<>;
#endif /* SETVALUE_HPP */
