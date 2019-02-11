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

    double MATH_PI = 3.14159265359; // CJ
    double Mp= 1.0/sqrt(8.0*MATH_PI);  // CJ
    // double Mass =  0.25 * 1e-11;  // CJ


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
        V_of_phi = 0.5 * pow(m_params.scalar_mass * vars.phi, 2.0);

        // The potential gradient at phi
        // m^2 phi
        dVdphi = pow(m_params.scalar_mass, 2.0) * vars.phi;

        // Starobinsky Higgs  // CJ
        V_of_phi =  (0.25 * 1e-11) * pow(1.0*Mp, 4) *
                   pow(1.0 - exp(-sqrt(2.0/3.0) * fabs(vars.phi)/Mp), 2); // CJ
        dVdphi =  2 * m_params.scalar_mass * pow(1.0*Mp,4) * sqrt(2.0/3.0) *
                  exp(- sqrt(2.0/3.0)* fabs(vars.phi) / Mp ) *
                  (1 - exp(- sqrt(2.0/3.0)* fabs(vars.phi) /Mp)) / Mp; // CJ

    }
};

#endif /* POTENTIAL_HPP_ */
