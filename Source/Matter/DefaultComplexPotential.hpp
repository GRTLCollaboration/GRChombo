/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DEFAULTCOMPLEXPOTENTIAL_HPP_
#define DEFAULTCOMPLEXPOTENTIAL_HPP_

#include "simd.hpp"

class DefaultComplexPotential
{
  public:
    //! The constructor
    DefaultComplexPotential() {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi_re,
                           data_t &dVdphi_im, const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        V_of_phi = 0.0;

        // The potential gradient at phi wrt real and im fields
        dVdphi_re = 0.0;
        dVdphi_im = 0.0;
    }
};

#endif /* DEFAULTCOMPLEXPOTENTIAL_HPP_ */
