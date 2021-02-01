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
        double f_axion;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        // 1/2 m^2 phi^2
        //V_of_phi = 0.5 * pow(m_params.scalar_mass * vars.phi, 2.0);

        // Axion potential
        V_of_phi = pow(m_params.scalar_mass * m_params.f_axion, 2.0) *
                   (1.0 - cos(vars.phi / m_params.f_axion));

        // The potential gradient at phi
        // m^2 phi
        //dVdphi = pow(m_params.scalar_mass, 2.0) * vars.phi;

        // Axion potential
        dVdphi = m_params.f_axion * pow(m_params.scalar_mass, 2.0) *
                 sin(vars.phi / m_params.f_axion);
    }
};

#endif /* POTENTIAL_HPP_ */
