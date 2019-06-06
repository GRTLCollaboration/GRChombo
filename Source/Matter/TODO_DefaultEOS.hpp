/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EOS_HPP_
#define EOSL_HPP_

#include "Tensor.hpp"
#include "simd.hpp"

class DefaultEquationOfState
{
  public:
    //! The constructor
    DefaultEquationOfStatel() {}

    //! Set the potential function for the scalar field here to zero
    template <class data_t, template <typename> class vars_t>
    void compute_eos(data_t &preassure, data_t &enthalpy,
                           const vars_t<data_t> &vars) const
    {
        // The default dust EOS
        pressure = 0.0;
        enthalpy = 0.0;  // density * velocity**2
    }
};

#endif /* DEFAULTPOTENTIAL_HPP_ */
