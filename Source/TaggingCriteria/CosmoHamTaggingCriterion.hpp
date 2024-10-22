/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COSMOHAMTAGGINGCRITERION_HPP_
#define COSMOHAMTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "CosmoAMR.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

class CosmoHamTaggingCriterion
{
  protected:
    const double m_dx;
    std::array<double, CH_SPACEDIM> m_tagging_center;
    double m_tagging_radius;
    double m_rho_mean;

  public:
    CosmoHamTaggingCriterion(double dx,
                             std::array<double, CH_SPACEDIM> tagging_center,
                             double tagging_radius, double rho_mean)
        : m_dx(dx), m_tagging_center(tagging_center),
          m_tagging_radius(tagging_radius), m_rho_mean(rho_mean){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto Ham_abs_sum = current_cell.load_vars(c_Ham_abs_sum);
        auto sqrt_gamma = current_cell.load_vars(c_sqrt_gamma);

        // Divide Ham_abs_sum by rho_mean
        data_t criterion = Ham_abs_sum / m_rho_mean * sqrt_gamma * m_dx;

        const Coordinates<data_t> coords(current_cell, m_dx, m_tagging_center);
        const data_t r = coords.get_radius();

        // Regrid where r is within tagging radius
        auto regrid = simd_compare_gt(r, m_tagging_radius);

        criterion = simd_conditional(regrid, 0.0, criterion);

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* COSMOHAMTAGGINGCRITERION_HPP_ */
