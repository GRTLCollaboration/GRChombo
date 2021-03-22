/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TRUNCATIONERRORTAGGING_HPP_
#define TRUNCATIONERRORTAGGING_HPP_

#include "Cell.hpp"

class TruncationErrorTagging
{
  protected:
    const int m_num_truncation_error_vars;
    // There are twice as many in vars as there are truncation error vars
    // first half: truncation error vars from this level
    // second half: truncation error vars interpolated from coarser level
    const int m_num_in_vars;
    const int m_level;
    const int m_store_var;

  public:
    TruncationErrorTagging(const int a_num_truncation_error_vars,
                           const int a_level, const int a_store_var = 0)
        : m_num_truncation_error_vars(a_num_truncation_error_vars),
          m_num_in_vars(2 * m_num_truncation_error_vars), m_level(a_level),
          m_store_var(a_store_var)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        data_t criterion;
        if (m_level == 0)
        {
            // refine everywhere on the coarsest level
            criterion = 100;
        }
        else
        {
            data_t in_vars[m_num_in_vars];
            for (int ivar = 0; ivar < m_num_in_vars; ++ivar)
            {
                in_vars[ivar] = current_cell.load_vars(ivar);
            }

            data_t truncation_errors[m_num_truncation_error_vars];
            for (int ivar = 0; ivar < m_num_truncation_error_vars; ++ivar)
            {
                // truncation error is difference between interpolated coarser
                // level data and this level
                truncation_errors[ivar] =
                    in_vars[ivar] - in_vars[m_num_truncation_error_vars + ivar];
            }

            data_t criterion_sq;
            for (int ivar = 0; ivar < m_num_truncation_error_vars; ++ivar)
            {
                criterion_sq +=
                    truncation_errors[ivar] * truncation_errors[ivar];
            }
            criterion = sqrt(criterion_sq);
        }
        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, m_store_var);
    }
};

#endif /* TRUNCATIONERRORTAGGING_HPP_ */
