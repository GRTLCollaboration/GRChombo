/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef HAMTAGGINGCRITERION_HPP_
#define HAMTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

class HamTaggingCriterion
{
  protected:
    const double m_dx;

  public:
    HamTaggingCriterion(double dx) : m_dx(dx){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto Ham_abs_sum = current_cell.load_vars(c_Ham_abs_sum);
        data_t criterion = sqrt(Ham_abs_sum) * m_dx;

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* HAMTAGGINGCRITERION_HPP_ */
