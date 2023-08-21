/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POLYNOMIAL_HPP_
#define POLYNOMIAL_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "UserVariables.hpp"

// This compute class sets the grid to a certain polynomial
class Polynomial
{
  public:
    Polynomial(const std::array<double, CH_SPACEDIM> &a_center, double a_dx)
        : m_dx(a_dx), m_center(a_center)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        data_t x = coords.x;
        data_t y = coords.y;
        data_t z = coords.z;

        // A is even in x and z, but in y it's even on the upper boundary
        data_t poliA = 42. + x * x + y * y * z * z;
        // B is odd in x
        data_t poliB = pow(x, 3);

        current_cell.store_vars(poliA, c_A);
        current_cell.store_vars(poliB, c_B);
    }

  protected:
    const double m_dx;
    const std::array<double, CH_SPACEDIM> &m_center;
};

#endif /* POLYNOMIAL_HPP_ */
