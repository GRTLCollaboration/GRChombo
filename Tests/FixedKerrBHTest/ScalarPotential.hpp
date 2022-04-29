/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SCALARPOTENTIAL_HPP_
#define SCALARPOTENTIAL_HPP_

#include "simd.hpp"

class ScalarPotential
{
  protected:
    const double m_mu;

  public:
    //! The constructor
    ScalarPotential(const double a_mu) : m_mu(a_mu) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        // 1/2 m^2 phi^2
        V_of_phi = 0.5 * m_mu * m_mu * vars.phi * vars.phi;

        // The potential gradient at phi wrt the field
        // m^2 phi
        dVdphi = m_mu * m_mu * vars.phi;
    }
};

#endif /* SCALARXPOTENTIAL_HPP_ */
