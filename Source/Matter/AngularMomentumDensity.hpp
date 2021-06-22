/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ANGULARMOMENTUMDENSITY_HPP_
#define ANGULARMOMENTUMDENSITY_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "MatterOnly.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the Angular Momentum Density with type matter_t and writes it to
//! the grid
template <class matter_t> class AngularMomentumDensity
{
    // Use the variable definition in MatterOnly (only need key vars)
    template <class data_t>
    using Vars = typename MatterOnly<matter_t>::template Vars<data_t>;

  protected:
    const FourthOrderDerivatives m_deriv; //!< An object for calculating derivs
    const matter_t my_matter;             //!< The matter object
    const double m_L;
    const double m_dx;

  public:
    AngularMomentumDensity(matter_t a_matter, double a_dx, double a_L)
        : my_matter(a_matter), m_deriv(a_dx), m_dx(a_dx), m_L(a_L)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and calc 1st
        // derivs
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);

        using namespace TensorAlgebra;
        const auto h_UU = compute_inverse_sym(vars.h);
        const auto chris = compute_christoffel(d1.h, h_UU);
        const emtensor_t<data_t> emtensor =
            my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

        // work out where we are on the grid
        std::array<double, 3> center;
        center.fill(0.5 * m_L);
        Coordinates<data_t> coords(current_cell, m_dx, center);

        // work out where we are on the grid
        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;

        Tensor<1, data_t> d_dphi;
        d_dphi[0] = -y;
        d_dphi[1] = x;
        d_dphi[2] = 0;

        // work out the component in phi direction
        data_t rhoJ = 0;
        FOR1(i) { rhoJ += emtensor.Si[i] * d_dphi[i]; }

        // store the density
        current_cell.store_vars(rhoJ, c_rhoJ);
    }
};

#endif /* ANGULARMOMENTUMDENSITY_HPP_ */
