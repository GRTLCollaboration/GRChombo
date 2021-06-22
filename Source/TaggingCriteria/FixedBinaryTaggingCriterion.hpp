/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBINARYTAGGINGCRITERION_HPP_
#define FIXEDBINARYTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

class FixedBinaryTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_L;
    const FourthOrderDerivatives m_deriv;
    const int m_level;
    const std::array<double, CH_SPACEDIM> m_center;
    const std::array<double, CH_SPACEDIM> m_centerA;
    const std::array<double, CH_SPACEDIM> m_centerB;

  public:
    FixedBinaryTaggingCriterion(const double dx, const int a_level,
                                const double a_L,
                                const std::array<double, CH_SPACEDIM> a_center,
                                const std::array<double, CH_SPACEDIM> a_centerA,
                                const std::array<double, CH_SPACEDIM> a_centerB)
        : m_dx(dx), m_deriv(dx), m_level(a_level), m_L(a_L), m_center(a_center),
          m_centerA(a_centerA), m_centerB(a_centerB){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        double criterion = 0.0;

        // make sure the inner part is regridded around the horizon
        double ratio = pow(2.0, -(m_level + 3.0));
        const Coordinates<data_t> coordsA(current_cell, m_dx, m_centerA);
        if (abs(coordsA.x) < m_L * ratio && abs(coordsA.y) < m_L * ratio &&
            abs(coordsA.z) < m_L * ratio)
        {
            criterion = 100;
        }
        const Coordinates<data_t> coordsB(current_cell, m_dx, m_centerB);
        if (abs(coordsB.x) < m_L * ratio && abs(coordsB.y) < m_L * ratio &&
            abs(coordsB.z) < m_L * ratio)
        {
            criterion = 100;
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* FIXEDBINARYTAGGINGCRITERION_HPP_ */
