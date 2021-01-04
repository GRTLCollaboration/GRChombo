/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This compute class enforces the positive chi and alpha condition
#ifndef POSITIVECHIANDALPHA_HPP_
#define POSITIVECHIANDALPHA_HPP_

#include "Cell.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"

class PositiveChiAndAlpha
{
  private:
    const double m_min_chi;
    const double m_min_lapse;

  public:
    //! Constructor for class
    PositiveChiAndAlpha(const double a_min_chi = 1e-4,
                        const double a_min_lapse = 1e-4)
        : m_min_chi(a_min_chi), m_min_lapse(a_min_lapse)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto chi = current_cell.load_vars(c_chi);
        auto lapse = current_cell.load_vars(c_lapse);

        chi = simd_max(chi, m_min_chi);
        lapse = simd_max(lapse, m_min_lapse);

        current_cell.store_vars(chi, c_chi);
        current_cell.store_vars(lapse, c_lapse);
    }
};

#endif /* POSITIVECHIANDALPHA_HPP_ */
