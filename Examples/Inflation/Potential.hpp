/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"

class Potential
{
  public:
    struct params_t
    {
        double scalar_mass;
    };

  private:
    params_t m_params;
    double m_mode;
    double m_L;

  public:
    //! The constructor
    Potential(params_t a_params, double a_L, double a_mode)
        : m_params(a_params), m_L(a_L), m_mode(a_mode)
    {
    }

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        // From rho = 1/2 (dphi/dx)^2 + V(phi) with V(phi) = 1/2 m^2 phi^2
        // ,we choose phi = A sin(2 n pi x/L) and m = 2 n pi/L such that initial
        // rho = constant
        double mass = 2 * m_mode * M_PI / m_L;
        V_of_phi = 0.5 * pow(mass * vars.phi, 2.0);

        // The potential gradient at phi
        // m^2 phi
        dVdphi = pow(m_params.scalar_mass, 2.0) * vars.phi;
    }
};

#endif /* POTENTIAL_HPP_ */
