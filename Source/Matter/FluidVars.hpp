/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FLUIDVARS_HPP_
#define FLUIDVARS_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "PerfectFluid.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the density rho with type matter_t and writes it to the grid
template <class matter_t> class FluidVars
{
    // Use the variable definition in MatterOnly - only require the key vars
    template <class data_t>
    using Vars = typename PerfectFluid<eos_t>::template Vars<data_t>;
    

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const matter_t my_matter; //!< The matter object

  public:
    Density(matter_t a_matter, double a_dx) : my_matter(a_matter), m_deriv(a_dx)
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

        // assign values of density in output box
        current_cell.store_vars(emtensor.rho, c_rho);
    }
};

#endif /* FLUIDVARS_HPP_ */
