/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALCONDITIONS_HPP_
#define INITIALCONDITIONS_HPP_

#include "ADMFixedBGVars.hpp"
#include "BoxLoops.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "KerrSchildFixedBG.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

template <class background_t> class InitialConditions
{
  public:
    InitialConditions(const background_t a_background,
                      const double a_self_interaction, const double a_L,
                      const double a_dx,
                      const std::array<double, CH_SPACEDIM> a_center)
        : m_dx(a_dx), m_center(a_center), m_L(a_L),
          m_self_interaction(a_self_interaction), m_background(a_background)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;

        // assign random data to Avec_i
        Tensor<1, data_t> Avec;
        Avec[0] =
            0.1053 + 0.10494 * cos(x) + 0.82769 * sin(y) - 0.21199 * sin(z);
        Avec[1] =
            0.15054 - 0.50088 * cos(x) - 0.15428 * sin(y) + 0.16779 * sin(z);
        Avec[2] =
            -0.02174 - 0.36243 * cos(x) + 0.81531 * cos(y) + 0.34918 * sin(z);

        // Store the initial values of the variables
        current_cell.store_vars(Avec[0], c_Avec1);
        current_cell.store_vars(Avec[1], c_Avec2);
        current_cell.store_vars(Avec[2], c_Avec3);
    }

  protected:
    const background_t m_background;
    const double m_self_interaction;
    const double m_dx;
    const double m_L;
    const std::array<double, CH_SPACEDIM> m_center;
};

#endif /* INITIALCONDITIONS_HPP_ */
