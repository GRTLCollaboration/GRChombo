/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef NANCHECK_HPP_
#define NANCHECK_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "UserVariables.hpp"

/// This compute class checks for nans or very large values and aborts the
/// execution when they appear displying a custom or default error message.
class NanCheck
{
  protected:
    const double m_dx;
    const std::array<double, CH_SPACEDIM> &m_center;
    const std::string m_error_info;
    const double m_max_abs;

  public:
    NanCheck(const std::string a_error_info = "NanCheck",
             const double a_max_abs = 1e20)
        : NanCheck(-1, {0.}, a_error_info, a_max_abs)
    {
    }

    NanCheck(const double a_dx, const std::array<double, CH_SPACEDIM> &a_center,
             const std::string a_error_info = "NanCheck",
             const double a_max_abs = 1e20)
        : m_dx(a_dx), m_center(a_center), m_error_info(a_error_info),
          m_max_abs(a_max_abs)
    {
    }

    void compute(Cell<double> current_cell) const
    {
        bool stop = false;
        int num_vars = current_cell.get_num_in_vars();
        for (int ivar = 0; ivar < num_vars; ++ivar)
        {
            double val;
            current_cell.load_vars(val, ivar);
            if (std::isnan(val) || abs(val) > m_max_abs)
                stop = true;
        }
        if (stop)
        {
#pragma omp single
            {
                pout() << m_error_info
                       << "::Values have become nan. The current state is: \n";
                for (int ivar = 0; ivar < num_vars; ++ivar)
                {
                    pout() << UserVariables::variable_names[ivar] << ": "
                           << current_cell.load_vars(ivar) << std::endl;
                }
                pout() << "Integer coordinates: " << current_cell.get_int_vect()
                       << std::endl;
                if (m_dx > 0)
                {
                    Coordinates<double> coords{current_cell, m_dx, m_center};
                    pout() << coords << std::endl;
                }
            }
            MayDay::Error("Values have become nan.");
        }
    }
};

#endif /* NANCHECK_HPP_ */
