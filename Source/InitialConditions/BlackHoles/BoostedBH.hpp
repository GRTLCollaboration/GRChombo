/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOOSTEDBH_HPP_
#define BOOSTEDBH_HPP_
/**
 * BOOSTED SCHWARZSCHILD BLACK HOLE
 * Baumgarte & Shapiro, pp. 73-74
 * NB: \bar{A} as defined in the book is psi^{-6} * \bar{A}_{BSSN}
 */

#include "Coordinates.hpp"
#include "Tensor.hpp"
#include <array>

class BoostedBH
{

  public:
    struct params_t
    {
        double mass;
        std::array<double, CH_SPACEDIM> center;
        std::array<double, CH_SPACEDIM> momentum;
    };

    const params_t m_params;

    BoostedBH(params_t a_params);

    // conformal factor
    template <class data_t>
    data_t psi_minus_one(Coordinates<data_t> coords) const;

    // extrinsic curvature
    template <class data_t>
    Tensor<2, data_t> Aij(Coordinates<data_t> coords) const;

  private:
    template <class data_t>
    data_t center_dist(Coordinates<data_t> coords) const;

    template <class data_t> data_t psi0(data_t r) const;

    template <class data_t> data_t psi2(data_t r, data_t cos_theta) const;

    template <class data_t> data_t psi2_0(data_t r) const;

    template <class data_t> data_t psi2_2(data_t r) const;
};

#include "BoostedBH.impl.hpp"

#endif /*BOOSTEDBH_HPP_*/
