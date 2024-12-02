/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which sets the initial scalar field matter config
class InitialScalarData
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double amplitude; //!< Amplitude of initial constant scalar field
        double bh_mass;   //!< Mass on the initial BH
        std::array<double, CH_SPACEDIM>
            center; //!< Centre of grid, for working out coords if neeeded
    };

    //! The constructor
    InitialScalarData(const params_t a_params, const Potential a_potential,
                      const double a_dx, const double a_G_Newton = 1.0)
        : m_dx(a_dx), m_params(a_params), m_potential(a_potential),
          m_G_Newton(a_G_Newton)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am I?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        data_t r = coords.get_radius();

        // set the field values, as constant in space
        ScalarField<Potential>::Vars<data_t> scalar_vars;
        scalar_vars.phi = m_params.amplitude;
        scalar_vars.Pi = 0.0;

        // Set the background for a Schwazschild BH in isotropic coords
        data_t psi = 1.0 + 0.5 * m_params.bh_mass / r;
        data_t chi = pow(psi, -4.0);

        // calculate the appropriate value of K to solve the constraints
        data_t V_of_phi, dVdphi;
        m_potential.compute_potential(V_of_phi, dVdphi, scalar_vars);
        data_t K_squared = 24.0 * M_PI * m_G_Newton * V_of_phi;
        data_t K = sqrt(K_squared);

        // store the vars
        current_cell.store_vars(scalar_vars);
        current_cell.store_vars(K, c_K);
        current_cell.store_vars(chi, c_chi);
        current_cell.store_vars(1.0, c_h11);
        current_cell.store_vars(1.0, c_h22);
        current_cell.store_vars(1.0, c_h33);
        current_cell.store_vars(1.0, c_lapse);
    }

  protected:
    const double m_dx;
    const double m_G_Newton;
    const Potential m_potential;
    const params_t m_params; //!< The matter initial condition params
};

#endif /* INITIALSCALARDATA_HPP_ */
