/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EOS_HPP_
#define EOS_HPP_

#include "simd.hpp"

class EquationOfState
{
  public:
    struct params_t
    {
      double omega;
      double mass;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    EquationOfState(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_eos(data_t &pressure, data_t &enthalpy,
                           const vars_t<data_t> &vars) const
   {

	       pressure = m_params.omega * vars.density;
	       enthalpy = m_params.mass + vars.energy + pressure / vars.density;

   }
};

#endif /* EquationOfState_HPP_ */
