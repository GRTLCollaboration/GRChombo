/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPHERICALHORIZON_HPP_
#define SPHERICALHORIZON_HPP_

#include "ADMConformalVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "simd.hpp"

//! This compute class calculates the outward expansion of null geodesics
//! given a center guess, and assuming spherical symmetry. The AH is where
//! Omega=0
class SphericalHorizon
{
    template <class data_t> using Vars = ADMConformalVars::VarsNoGauge<data_t>;

  protected:
    const FourthOrderDerivatives m_deriv;
    const double m_dx;
    const std::array<double, 3> m_center; // only works in 3+1D

  public:
    SphericalHorizon(double dx, std::array<double, 3> a_center)
        : m_center(a_center), m_dx(dx), m_deriv(dx)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    template <class data_t>
    data_t calculate_expansion(const Vars<data_t> &vars,
                               const Vars<Tensor<1, data_t>> &d1,
                               const Coordinates<data_t> coords) const;
};

#include "SphericalHorizon.impl.hpp"

#endif /* SPHERICALHORIZON_HPP_ */
