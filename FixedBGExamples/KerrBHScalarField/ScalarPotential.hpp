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
    //    const double m_mu;
    const InitialScalarData::params_t m_initial_params;

  public:
    //! The constructor
    ScalarPotential(const InitialScalarData::params_t a_initial_params)
        : m_initial_params(a_initial_params)
    {
    }

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        // 1/2 m^2 phi^2
        const double mu = m_initial_params.mass;
        V_of_phi = 0.5 * mu * mu * vars.phi * vars.phi;

        // The potential gradient at phi wrt the field
        // m^2 phi
        dVdphi = mu * mu * vars.phi;
    }
};

#endif /* SCALARXPOTENTIAL_HPP_ */
