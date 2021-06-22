/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef OSCILLOTONTAGGINGCRITERION_HPP_
#define OSCILLOTONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtractionComplex.hpp"
#include "Tensor.hpp"

class OscillotonTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_chi;
    const double m_threshold_Pi;
    const int m_level;
    const SphericalExtractionComplex::params_t m_params;

  public:
    OscillotonTaggingCriterion(
        const double dx, const double threshold_chi, const double threshold_Pi,
        const int a_level, const SphericalExtractionComplex::params_t a_params)
        : m_dx(dx), m_deriv(dx), m_threshold_chi(threshold_chi),
          m_threshold_Pi(threshold_Pi), m_params(a_params), m_level(a_level){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto chi = current_cell.load_vars(c_chi);
        Tensor<1, data_t> d1_chi;
        FOR1(idir) m_deriv.diff1(d1_chi, current_cell, idir, c_chi);

        auto Pi = current_cell.load_vars(c_Pi);
        Tensor<1, data_t> d1_Pi;
        FOR1(idir) m_deriv.diff1(d1_Pi, current_cell, idir, c_Pi);

        // to stop zooming out when Pi goes flat but phi large we apply
        // condition to both values
        auto phi = current_cell.load_vars(c_phi);
        Tensor<1, data_t> d1_phi;
        FOR1(idir) m_deriv.diff1(d1_phi, current_cell, idir, c_phi);

        data_t mod_d1_chi = 0;
        data_t mod_d1_Pi = 0;
        data_t mod_d1_phi = 0;
        FOR1(idir)
        {
            mod_d1_chi += d1_chi[idir] * d1_chi[idir];
            mod_d1_Pi += d1_Pi[idir] * d1_Pi[idir];
            mod_d1_phi += d1_phi[idir] * d1_phi[idir];
        }

        data_t criterion =
            m_dx *
            (1.0 / m_threshold_chi * (sqrt(mod_d1_chi) / pow(chi, 2.0)) +
             1.0 / m_threshold_Pi * (sqrt(mod_d1_phi) + sqrt(mod_d1_Pi)));

        // regrid if within extraction level and not at required refinement
        if (m_level < m_params.extraction_level)
        {
            const Coordinates<data_t> coords(current_cell, m_dx,
                                             m_params.extraction_center);
            /*          // For fixed boxes not radius
                        data_t regridx =
                            simd_compare_lt(abs(coords.x -
               m_params.extraction_center[0]), m_params.extraction_radius*1.1);
                        data_t regridy =
                            simd_compare_lt(abs(coords.y -
               m_params.extraction_center[1]), m_params.extraction_radius*1.1);
                        data_t regridz =
                            simd_compare_lt(abs(coords.z -
               m_params.extraction_center[2]), m_params.extraction_radius*1.1);
                        criterion = simd_conditional(regridx, 1.0, criterion);
                        criterion = simd_conditional(regridy, 1.0, criterion);
                        criterion = simd_conditional(regridz, 1.0, criterion);
            */
            const data_t r = coords.get_radius();
            // add a 20% buffer to extraction zone so not too near to boundary
            data_t regrid =
                simd_compare_lt(r, 1.2 * m_params.extraction_radius);
            criterion = simd_conditional(regrid, 1.0, criterion);
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* OSCILLOTONTAGGINGCRITERION_HPP_ */
