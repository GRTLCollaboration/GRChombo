/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PHIANDCHITAGGINGCRITERION_HPP_
#define PHIANDCHITAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

class PhiAndChiTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_L;
    const double m_extract_radius;
    const int m_level;
    const int m_extraction_level;
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_phi;
    const double m_threshold_chi;

  public:
    PhiAndChiTaggingCriterion(const double dx, const double threshold_phi,
                              const double threshold_chi, const double L,
                              const double extract_radius, const int a_level,
                              const int a_extraction_level)
        : m_dx(dx), m_deriv(dx), m_threshold_phi(threshold_phi),
          m_threshold_chi(threshold_chi), m_L(L),
          m_extract_radius(extract_radius), m_level(a_level),
          m_extraction_level(a_extraction_level){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Tensor<1, data_t> d1_phi;
        FOR1(idir) m_deriv.diff1(d1_phi, current_cell, idir, c_phi);

        auto chi = current_cell.load_vars(c_chi);
        Tensor<1, data_t> d1_chi;
        FOR1(idir) m_deriv.diff1(d1_chi, current_cell, idir, c_chi);

        data_t mod_d1_phi = 0;
        data_t mod_d1_chi = 0;
        FOR1(idir)
        {
            mod_d1_phi += d1_phi[idir] * d1_phi[idir];
            mod_d1_chi += d1_chi[idir] * d1_chi[idir];
        }

        data_t criterion1 =
            m_dx * (sqrt(mod_d1_phi) / m_threshold_phi +
                    sqrt(mod_d1_chi) / m_threshold_chi * pow(chi, -1.0));

        if (m_level <= m_extraction_level)
        {
            std::array<double, CH_SPACEDIM> center;
            center.fill(0.5 * m_L);
            const Coordinates<data_t> coords(current_cell, m_dx, center);
            const data_t r = coords.get_radius();
            const data_t regrid = simd_compare_lt(r, m_extract_radius);
            data_t criterion2 = simd_conditional(regrid, 1.0, criterion1);

            // Write back into the flattened Chombo box
            current_cell.store_vars(criterion2, 0);
        }
        else
        {
            current_cell.store_vars(criterion1, 0);
        }
    }
};

#endif /* PHIANDCHITAGGINGCRITERION_HPP_ */
