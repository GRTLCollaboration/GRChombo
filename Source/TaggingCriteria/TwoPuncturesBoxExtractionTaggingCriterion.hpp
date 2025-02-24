/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TWOPUNCTURESBOXEXTRACTIONTAGGINGCRITERION_HPP_
#define TWOPUNCTURESBOXEXTRACTIONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"

class TwoPuncturesBoxExtractionTaggingCriterion
{
  protected:
    const double m_dx;
    const int m_level;
    const int m_max_level;
    const spherical_extraction_params_t m_params;
    const bool m_activate_extraction;
    const std::vector<double> m_puncture_masses;
    const std::vector<std::array<double, CH_SPACEDIM>> &m_puncture_coords;

  public:
    TwoPuncturesBoxExtractionTaggingCriterion(
        const double dx, const int a_level, const int a_max_level,
        const spherical_extraction_params_t a_params,
        const std::vector<std::array<double, CH_SPACEDIM>> &a_puncture_coords,
        const bool activate_extraction = false,
        const std::vector<double> a_puncture_masses = {0.5, 0.5})
        : m_dx(dx), m_level(a_level), m_max_level(a_max_level),
          m_params(a_params), m_activate_extraction(activate_extraction),
          m_puncture_masses(a_puncture_masses),
          m_puncture_coords(a_puncture_coords){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        data_t criterion = 0.0;

        double puncture_separation_squared = 0.0;
        FOR(i)
        {
            double displacement =
                m_puncture_coords[0][i] - m_puncture_coords[1][i];
            puncture_separation_squared += displacement * displacement;
        }
        // make sure the inner part is regridded around the horizon
        // take L as the length of full grid, so tag inner 1/2
        // of it, which means inner \pm L/4
        // we want each level to be double the innermost one in size
        const double factor = pow(2.0, m_max_level - m_level - 1);
        double sum_masses = m_puncture_masses[0] + m_puncture_masses[1];

        if (puncture_separation_squared > 2.25 * sum_masses * sum_masses)
        {
            // loop over puncture masses
            for (int ipuncture = 0; ipuncture < m_puncture_masses.size();
                 ++ipuncture)
            {
                // where am i?
                const Coordinates<data_t> coords(current_cell, m_dx,
                                                 m_puncture_coords[ipuncture]);
                const data_t max_abs_xy =
                    simd_max(abs(coords.x), abs(coords.y));
                const data_t max_abs_xyz = simd_max(max_abs_xy, abs(coords.z));
                auto regrid = simd_compare_lt(
                    max_abs_xyz, 2.5 * factor * m_puncture_masses[ipuncture]);
                criterion = simd_conditional(regrid, 100.0, criterion);
            }
        }
        else if (puncture_separation_squared < 2.25 * sum_masses * sum_masses &&
                 m_level == m_max_level - 1)
        {
            criterion = 0;
        }
        else if (puncture_separation_squared < 2.25 * sum_masses * sum_masses &&
                 m_level == m_max_level - 2)
        {
            // where am i?
            std::array<double, CH_SPACEDIM> puncture_centre;
            FOR(i)
            puncture_centre[i] =
                0.5 * (m_puncture_coords[0][i] + m_puncture_coords[1][i]);
            const Coordinates<data_t> coords(current_cell, m_dx,
                                             puncture_centre);
            const data_t max_abs_xy = simd_max(abs(coords.x), abs(coords.y));
            const data_t max_abs_xyz = simd_max(max_abs_xy, abs(coords.z));
            auto regrid =
                simd_compare_lt(max_abs_xyz, 2.5 * factor * sum_masses);
            criterion = simd_conditional(regrid, 100.0, criterion);
        }

        // if extracting weyl data at a given radius, enforce a given resolution
        // there
        if (m_activate_extraction)
        {
            for (int iradius = 0; iradius < m_params.num_extraction_radii;
                 ++iradius)
            {
                // regrid if within extraction level and not at required
                // refinement
                if (m_level < m_params.extraction_levels[iradius])
                {
                    const Coordinates<data_t> coords(current_cell, m_dx,
                                                     m_params.center);
                    const data_t r = coords.get_radius();
                    // add a 20% buffer to extraction zone so not too near to
                    // boundary
                    auto regrid = simd_compare_lt(
                        r, 1.2 * m_params.extraction_radii[iradius]);
                    criterion = simd_conditional(regrid, 100.0, criterion);
                }
            }
        }
        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* TWOPUNCTURESBOXEXTRACTIONTAGGINGCRITERION_HPP_ */
