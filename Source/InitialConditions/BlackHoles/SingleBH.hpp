/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SINGLEBH_HPP_
#define SINGLEBH_HPP_

#include "BoostedBH.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"
#include <array>

enum Lapse
{
    ONE,
    PRE_COLLAPSED,
    CHI
};

// Just a wrapper of the BoostedBH class, to be used by itself (instead of a
// BinaryBH)
class SingleBH
{
  protected:
    double m_dx;
    BoostedBH bh;
    int m_initial_lapse;

  public:
    SingleBH(BoostedBH::params_t a_bh_params, double a_dx,
             int a_initial_lapse = Lapse::PRE_COLLAPSED)
        : m_dx(a_dx), bh(a_bh_params), m_initial_lapse(a_initial_lapse)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    template <class data_t>
    data_t compute_chi(Coordinates<data_t> coords) const;

    template <class data_t>
    Tensor<2, data_t> compute_A(data_t chi, Coordinates<data_t> coords) const;
};

#include "SingleBH.impl.hpp"

#endif /* SINGLEBH_HPP_ */
