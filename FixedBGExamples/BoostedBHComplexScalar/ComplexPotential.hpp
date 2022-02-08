/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXPOTENTIAL_HPP_
#define COMPLEXPOTENTIAL_HPP_

#include "simd.hpp"

class ComplexPotential
{
  protected:
    const double m_mu;

  public:
    //! The constructor
    ComplexPotential(const double a_mu) : m_mu(a_mu) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi_re,
                           data_t &dVdphi_im, const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        // 1/2 m^2 phi^2
        V_of_phi = 0.5 * m_mu * m_mu * vars.phi_Re * vars.phi_Re +
                   0.5 * m_mu * m_mu * vars.phi_Im * vars.phi_Im;

        // The potential gradient at phi wrt the real and im field
        // m^2 phi
        dVdphi_re = m_mu * m_mu * vars.phi_Re;
        dVdphi_im = m_mu * m_mu * vars.phi_Im;
    }
};

#endif /* COMPLEXPOTENTIAL_HPP_ */
