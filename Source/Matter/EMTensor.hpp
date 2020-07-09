/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EMTENSOR_HPP_
#define EMTENSOR_HPP_

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

//! Calculates the EMTensor with type matter_t and writes it to the grid
template <class matter_t> class EMTensor
{
    // Use the variable definition in MatterOnly - only require the key vars
    template <class data_t>
    using Vars = typename MatterOnly<matter_t>::template Vars<data_t>;

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const matter_t my_matter; //!< The matter object

  public:
    EMTensor(matter_t a_matter, double a_dx)
        : my_matter(a_matter), m_deriv(a_dx)
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

        // assign values of EMTensor in output box
        current_cell.store_vars(emtensor.rho, c_rho);
        current_cell.store_vars(emtensor.S, c_S);
        current_cell.store_vars(emtensor.Si[0], c_S1);
        current_cell.store_vars(emtensor.Si[1], c_S2);
        current_cell.store_vars(emtensor.Si[2], c_S3);
    }
};

#endif /* EMTENSOR_HPP_ */
