/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXTRACTIONTAGGINGCRITERION_HPP_
#define EXTRACTIONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

class ExtractionTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_L;
    const double m_extract_radius;
    int m_level;
    int m_extraction_level;

  public:
    ExtractionTaggingCriterion(double dx, double L, double extract_radius,
                               int a_level, int a_extraction_level)
        : m_dx(dx), m_L(L), m_extract_radius(extract_radius), m_level(a_level),
          m_extraction_level(a_extraction_level){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        data_t criterion = 0.0; // current_cell.load_vars(0);
        if (m_level < m_extraction_level)
        {
            std::array<double, CH_SPACEDIM> center;
            center.fill(0.5 * m_L);
            const Coordinates<data_t> coords(current_cell, m_dx, center);
            const data_t r = coords.get_radius();
            data_t regrid = simd_compare_lt(r, m_extract_radius);
            criterion = simd_conditional(regrid, 1.0, criterion);
        }
        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* EXTRACTIONTAGGINGCRITERION_HPP_ */
