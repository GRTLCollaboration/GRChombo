/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(HARMONICTEST_HPP_)
#error "This file should only be included through HarmonicTest.hpp"
#endif

#ifndef HARMONICTEST_IMPL_HPP_
#define HARMONICTEST_IMPL_HPP_

#include "DebuggingTools.hpp"
#include "HarmonicTest.hpp"
#include "SphericalHarmonics.hpp"
#include "simd.hpp"

template <class data_t>
void HarmonicTest::compute(Cell<data_t> current_cell) const
{

    Coordinates<data_t> coords(current_cell, m_dx, m_center_vector);

    data_t phi = compute_harmonic(coords);

    data_t radius1 = coords.get_radius();
    data_t radius2 =
        Coordinates<data_t>::get_radius(current_cell, m_dx, m_center_vector);
    phi = phi / radius1 / radius2;

    current_cell.store_vars(phi, c_phi);
}

template <class data_t>
data_t HarmonicTest::compute_harmonic(Coordinates<data_t> coords) const
{

    // Add in el, em spherical harmonics here, spin weight es
    using namespace SphericalHarmonics;
    int es = -1;
    int el = 2;
    int em = -1;
    auto Y_lm = spin_Y_lm(coords.x, coords.y, coords.z, es, el, em);
    data_t out = Y_lm.Real;

    return out;
}

#endif /* HARMONICTEST_IMPL_HPP_ */
