/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef OSCILLOTONTAGGINGCRITERION_HPP_
#define OSCILLOTONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "ScalarField.hpp"
#include "SimulationParametersBase.hpp"
#include "Tensor.hpp"

class OscillotonTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_L;
    const FourthOrderDerivatives m_deriv;
    const int m_level;
    const int m_min_level;
    const std::array<double, CH_SPACEDIM> m_center;

    /// Vars object for chi
    template <class data_t> struct Vars
    {
        data_t chi; //!< Conformal factor

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
        }
    };

  public:
    OscillotonTaggingCriterion(const double dx, const int a_min_level,
                               const int a_level, const double a_L,
                               const std::array<double, CH_SPACEDIM> a_center)
        : m_dx(dx), m_deriv(dx), m_min_level(a_min_level), m_L(a_L),
          m_level(a_level), m_center(a_center){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto d2chi = m_deriv.template diff2<Vars>(current_cell);

        data_t mod_d2_chi = 0;

        FOR2(idir, jdir)
        {
           mod_d2_chi += d2chi.chi[idir][jdir] * d2chi.chi[idir][jdir]; 
        }

        // criterion is primarily based on whether chi gradients high
        data_t criterion = m_dx * sqrt(mod_d2_chi);

        // make sure the inner part is regridded up to some level
        // take L as the length of full grid, so tag inner 1/2
        // of it, which means inner \pm L/4 (on level zero)
        const double ratio = pow(2.0, -(m_level + 2.0));
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        if (abs(coords.y) < m_L * ratio && abs(coords.z) < m_L * ratio
                                        && m_level <= m_min_level)
        {
            auto should_refine = simd_compare_lt(abs(coords.x), m_L * ratio);
            // just make the criterion big so Chombo definitely tags it
            criterion = simd_conditional(should_refine, 100.0, criterion);
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* OSCILLOTONTAGGINGCRITERION_HPP_ */
