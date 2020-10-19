/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDGRIDSTAGGINGCRITERION_HPP_
#define FIXEDGRIDSTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

class FixedGridsTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_L;
    const FourthOrderDerivatives m_deriv;
    const int m_level;
    const std::array<double, CH_SPACEDIM> m_center;

  public:
    FixedGridsTaggingCriterion(const double dx, const int a_level,
                               const double a_L,
                               const std::array<double, CH_SPACEDIM> a_center)
        : m_dx(dx), m_deriv(dx), m_level(a_level), m_L(a_L),
          m_center(a_center){};

    void compute(Cell<double> current_cell) const
    {
        double criterion = 0.0;

        // make sure the inner part is regridded around the horizon
        // take L as the length of full grid, so tag inner 1/2
        // of it, which means inner \pm L/4
        double ratio = pow(2.0, -(m_level + 2.0));
        const Coordinates<double> coords(current_cell, m_dx, m_center);
        if (abs(coords.x) < m_L * ratio && abs(coords.y) < m_L * ratio &&
            abs(coords.z) < m_L * ratio)
        {
            criterion = 100;
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* FIXEDGRIDSTAGGINGCRITERION_HPP_ */
