/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TRUNCATIONERROREXTRACTIONTAGGING_HPP_
#define TRUNCATIONERROREXTRACTIONTAGGING_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "SphericalExtraction.hpp"

class TruncationErrorExtractionTagging
{
  protected:
    const double m_dx;
    const int m_num_truncation_error_vars;
    // There are twice as many in vars as there are truncation error vars
    // first half: truncation error vars from this level
    // second half: truncation error vars interpolated from coarser level
    const int m_num_in_vars;
    const int m_level;
    const bool m_activate_extraction;
    const SphericalExtraction::params_t m_extraction_params;
    const int m_store_var;

  public:
    TruncationErrorExtractionTagging(
        const double a_dx, const int a_num_truncation_error_vars,
        const int a_level, const bool a_activate_extraction,
        const SphericalExtraction::params_t &a_extraction_params,
        const int a_store_var = 0)
        : m_dx(a_dx), m_num_truncation_error_vars(a_num_truncation_error_vars),
          m_num_in_vars(2 * m_num_truncation_error_vars), m_level(a_level),
          m_activate_extraction(a_activate_extraction),
          m_extraction_params(a_extraction_params), m_store_var(a_store_var)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        data_t criterion;
        if (m_level == 0)
        {
            // refine everywhere on the coarsest level
            criterion = 100;
        }
        else
        {
            data_t in_vars[m_num_in_vars];
            for (int ivar = 0; ivar < m_num_in_vars; ++ivar)
            {
                in_vars[ivar] = current_cell.load_vars(ivar);
            }

            data_t truncation_errors[m_num_truncation_error_vars];
            for (int ivar = 0; ivar < m_num_truncation_error_vars; ++ivar)
            {
                // truncation error is difference between interpolated
                // coarser level data and this level
                truncation_errors[ivar] =
                    in_vars[ivar] - in_vars[m_num_truncation_error_vars + ivar];
            }

            data_t criterion_sq;
            for (int ivar = 0; ivar < m_num_truncation_error_vars; ++ivar)
            {
                criterion_sq +=
                    truncation_errors[ivar] * truncation_errors[ivar];
            }
            criterion = sqrt(criterion_sq);

            // if extracting on spheres at a given radius, enforce a given
            // resolution there
            if (m_activate_extraction)
            {
                for (int iradius = 0;
                     iradius < m_extraction_params.num_extraction_radii;
                     ++iradius)
                {
                    // regrid if within extraction level and not at required
                    // refinement
                    if (m_level <
                        m_extraction_params.extraction_levels[iradius])
                    {
                        const Coordinates<data_t> coords(
                            current_cell, m_dx, m_extraction_params.center);
                        const data_t r = coords.get_radius();
                        // add a 20% buffer to extraction zone so not too
                        // near to boundary
                        auto regrid = simd_compare_lt(
                            r,
                            1.2 *
                                m_extraction_params.extraction_radii[iradius]);
                        criterion = simd_conditional(regrid, 100.0, criterion);
                    }
                }
            }
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, m_store_var);
    }
};

#endif /* TRUNCATIONERROREXTRACTIONTAGGING_HPP_ */
