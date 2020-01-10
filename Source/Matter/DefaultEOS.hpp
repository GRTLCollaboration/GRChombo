/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DEFAULTEOS_HPP_
#define DEFAULTEOS_HPP_

#include "Tensor.hpp"
#include "simd.hpp"

class DefaultEOS
{
  public:
    //! The constructor
    DefaultEOS() {}

    //! Set the potential function for the scalar field here to zero
    template <class data_t, template <typename> class vars_t>
    void compute_eos(data_t &pressure, data_t &enthalpy,
                           const vars_t<data_t> &vars) const
    {
        // The default dust EOS
        pressure = 0.0;
        enthalpy = 1.0;  // mass + energy + pressure/density
    }
};

#endif /* DEFAULTPOTENTIAL_HPP_ */
