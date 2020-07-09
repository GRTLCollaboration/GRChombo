/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISION_HPP_
#define EXCISION_HPP_

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
template <class matter_t, class background_t> class Excision
{
    // Use matter_t class
    using Vars = typename matter_t::template Vars<double>;

  protected:
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The BH center
    const FourthOrderDerivatives m_deriv;
    const background_t m_background;

  public:
    Excision(const double a_dx, const std::array<double, CH_SPACEDIM> a_center,
             background_t a_background)
        : m_dx(a_dx), m_deriv(m_dx), m_center(a_center),
          m_background(a_background)
    {
    }

    void compute(const Cell<double> current_cell) const
    {
        double horizon_distance = m_background.excise(current_cell);
        if (horizon_distance < 1.0)
        {
            // the matter rhs vars within the excision zone
            // recalculate them - for now set to decay to zero
            Vars rhs;
            rhs.enum_mapping(
                [this, &current_cell](const int &ivar, double &var) {
                    var = -0.01 * var;
                });

            // assign values of rhs in output box
            current_cell.store_vars(rhs);
        } // else do nothing
    }
};

#endif /* EXCISION_HPP_ */
