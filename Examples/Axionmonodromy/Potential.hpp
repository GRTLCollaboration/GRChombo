/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <iostream>
#include "parstream.H"

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"

class Potential
{
  public:
    struct params_t
    {
        double overall_normalization;
	double decay_constant;
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
        // 1/(8 pi) M^2[1-1/sqrt(1+8 pi phi^2/M^2)]
        double pi = M_PI;
        V_of_phi =(1/(8.0*M_PI))*pow(m_params.overall_normalization,2.0)*(1.0-1.0/sqrt(1.0+8.0*M_PI*pow(vars.phi,2.0)/pow(m_params.overall_normalization,2.0)));
        // The potential gradient at phi
        // phi/(1+8 pi phi^2/M^2)^(3/2)
        dVdphi = vars.phi/pow(1.0+8.0*M_PI*pow(vars.phi,2.0)/(pow(m_params.overall_normalization,2.0)),3.0/2.0);
    }
};

#endif /* POTENTIAL_HPP_ */
