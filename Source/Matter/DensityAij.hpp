/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DENSITYAIJ_HPP_
#define DENSITYAIJ_HPP_

#include "CCZ4.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates a quantity related to the local, non averaged Isaacson
//! 'energy density' and writes it to the grid
class DensityAij
{
    // Use the variable definition without gauge
    template <class data_t> using Vars = CCZ4::Vars<data_t>;

  protected:
    const double m_G_Newton; //!< The value of Newtons constant

  public:
    DensityAij(double a_G_Newton = 1.0) : m_G_Newton(a_G_Newton) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables
        const auto vars = current_cell.template load_vars<Vars>();

        using namespace TensorAlgebra;
        const auto h_UU = compute_inverse_sym(vars.h);
        const auto A_UU = raise_all(vars.A, h_UU);
        const data_t tr_A2 = compute_trace(vars.A, A_UU);

        // the value of rhoAij = A^ij A_ij /8 pi G
        data_t rhoAij = tr_A2 / 8.0 / M_PI / m_G_Newton;

        // assign values of density in output box
        current_cell.store_vars(rhoAij, c_rhoAij);
    }
};

#endif /* DENSITYAIJ_HPP_ */
