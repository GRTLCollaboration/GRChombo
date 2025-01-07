/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "EMTensor.hpp" //< To compute matter quantit
#include "MatterCCZ4RHS.hpp"
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
        double amplitude; //!< Amplitude of bump in initial SF bubble
        std::array<double, CH_SPACEDIM>
            center; //!< Centre of perturbation in initial SF bubble
    };

    //! The constructor
    InitialScalarData(params_t a_params, double a_dx, double a_scalar_mass)
        : m_dx(a_dx), m_params(a_params), m_scalar_mass(a_scalar_mass)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

        // From rho = 1/2 (dphi/dx)^2 + V(phi) with V(phi) = 1/2 m^2 phi^2
        // ,we choose phi = A sin(2 n pi x/L) and m = 2 n pi/L such that initial
        // rho = constant
        // Calculate the field value
        data_t phi = m_params.amplitude * sin(m_scalar_mass * coords.x);

        // store the vars
        current_cell.store_vars(phi, c_phi);
        current_cell.store_vars(0.0, c_Pi);
        // Additional initial vars for flat cosmological background
        current_cell.store_vars(1.0, c_lapse);
        current_cell.store_vars(1.0, c_h11);
        current_cell.store_vars(1.0, c_h22);
        current_cell.store_vars(1.0, c_h33);
        current_cell.store_vars(1.0, c_chi);
    }

    // Function to analytically compute initial rho mean
    double compute_initial_rho_mean()
    {
        // Compute initial rho_mean here:
        // from rho = 1/2 m^2 phi^2
        // phi0 = A sin(2 pi n x /L)
        // mean of phi0^2 = A^2 sin^2 = 0.5*A^2
        double rho_mean = 0.5 * m_scalar_mass * m_scalar_mass *
                          (0.5 * m_params.amplitude * m_params.amplitude);
        return rho_mean;
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
    double m_scalar_mass;    //!< Mass of Scalar Field
};

#endif /* INITIALSCALARDATA_HPP_ */
