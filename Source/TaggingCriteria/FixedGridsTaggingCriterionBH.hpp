/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDGRIDSTAGGINGCRITERIONBH_HPP_
#define FIXEDGRIDSTAGGINGCRITERIONBH_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

class FixedGridsTaggingCriterionBH
{
  protected:
    const double m_dx;
    const double m_L;
    const FourthOrderDerivatives m_deriv;
    const int m_level;
    const int m_max_level;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_radius;

  public:
    FixedGridsTaggingCriterionBH(const double dx, const int a_level,
                                 const int a_max_level, const double a_L,
                                 const std::array<double, CH_SPACEDIM> a_center,
                                 const double a_radius = 1.75)
        : m_dx(dx), m_deriv(dx), m_level(a_level), m_max_level(a_max_level),
          m_L(a_L), m_center(a_center), m_radius(a_radius){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        double criterion = 0.0;

        // make sure the inner part is regridded around the horizon
        // take L as the length of full grid, so tag inner 1/2
        // of it, which means inner \pm L/4
        double ratio = pow(2.0, -(m_level + 2.0));
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        if (abs(coords.x) < m_L * ratio && abs(coords.y) < m_L * ratio &&
            abs(coords.z) < m_L * ratio)
        {
            criterion = 100;
        }

        // make sure horizon always covered
        if (m_level == m_max_level - 1)
        {
            if (abs(coords.x) < m_radius && abs(coords.y) < m_radius &&
                abs(coords.z) < m_radius)
            {
                criterion = 100.0;
            }
        }
        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* FIXEDGRIDSTAGGINGCRITERIONBH_HPP_ */
