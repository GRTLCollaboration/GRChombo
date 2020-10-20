/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates a bubble of a scalar field given params for initial
//! matter config
class InitialScalarData
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double amplitude; //!< Amplitude of bump in initial SF bubble
        double omega;     // oscillation frequency of scalar
        std::array<double, CH_SPACEDIM>
            center; //!< Centre of perturbation in initial SF bubble
    };

  protected:
    const double m_dx;
    const params_t m_params; //!< The matter initial condition params

  public:
    //! The constructor
    InitialScalarData(params_t a_params, double a_dx)
        : m_dx(a_dx), m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

        // where am I?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        const data_t R = coords.get_radius();
        const data_t rho2 =
            simd_max(coords.x * coords.x + coords.y * coords.y, 1e-8);
        const data_t rho = sqrt(rho2);
        const data_t cos_phi = coords.x / rho;
        const data_t sin_phi = coords.y / rho;
        const data_t sin_theta = rho / R;
        const data_t cos_theta = coords.z / R;
        const double M = 0.5;
        const double a = 0.0;
        const double r_plus = M + sqrt(M * M - a * a);

        // the Boyer Lindquist coord (eqn (46) 1401.1548)
        const data_t r_BL = R + 0.5 * r_plus + 0.0625 * r_plus * r_plus / R;
        const double Mmu = m_params.omega;

        // set the field vars - so real and im parts out of phase
        data_t r0 = M / Mmu / Mmu;
        data_t R_of_r = 1.512 / cosh(r_BL / 22.5) / cosh(r_BL / 22.5) -
                        2.85 * exp(-r_BL / 2.43);
        data_t Theta_of_theta =
            0.86 * exp(0.9 * 0.0969013 * cos_theta) * sin_theta;
        data_t phi = m_params.amplitude * R_of_r * Theta_of_theta * cos_phi;
        data_t Pi = -m_params.omega * m_params.amplitude * R_of_r *
                    Theta_of_theta * sin_phi;

        // Store the initial values of the variables
        current_cell.store_vars(phi, c_phi);
        current_cell.store_vars(Pi, c_Pi);
    }
};

#endif /* INITIALSCALARDATA_HPP_ */
