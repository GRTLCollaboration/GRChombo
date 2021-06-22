/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISIONPROCATEST_HPP_
#define EXCISIONPROCATEST_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Does excision for fixed BG BH solutions
//! Note that it is does not using simd so one must set disable_simd()
template <class matter_t, class background_t> class ExcisionProcaTest
{
    // Use matter_t class
    using Vars = typename matter_t::template Vars<double>;

  protected:
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The BH center
    const background_t m_background;

  public:
    ExcisionProcaTest(const double a_dx,
                      const std::array<double, CH_SPACEDIM> a_center,
                      background_t a_background)
        : m_dx(a_dx), m_center(a_center), m_background(a_background)
    {
    }

    void compute(const Cell<double> current_cell) const
    {
        const Coordinates<double> coords(current_cell, m_dx, m_center);
        // only check within the middle part of the grid
        // as it can't be more than 2GM (assumes M < 1)

        double horizon_distance = m_background.excise(current_cell);
        if (horizon_distance < 1.0)
        {
            current_cell.store_vars(0.0, c_gauss);
            current_cell.store_vars(0.0, c_Z);
        } // else do nothing
    }
};

#endif /* EXCISIONPROCATEST_HPP_ */
