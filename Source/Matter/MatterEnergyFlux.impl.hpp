/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MATTERENERGYFLUX_HPP)
#error "This file should only be included through MatterEnergyFlux.hpp"
#endif

#ifndef MATTERENERGYFLUX_IMPL_HPP
#define MATTERENERGYFLUX_IMPL_HPP

template <class matter_t>
MatterEnergyFlux<matter_t>::MatterEnergyFlux(
    matter_t a_matter, const std::array<double, CH_SPACEDIM> &a_center,
    const double a_dx, const int a_var_enum)
    : m_matter(a_matter), m_center(a_center), m_dx(a_dx), m_deriv(a_dx),
      m_var_enum(a_var_enum)
{
}

template <class matter_t>
template <class data_t>
void MatterEnergyFlux<matter_t>::compute(Cell<data_t> current_cell) const
{
    const auto matter_vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);

    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(matter_vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // Calculate EM tensor
    const auto emtensor =
        m_matter.compute_emtensor(matter_vars, d1, h_UU, chris.ULL);

    const Coordinates<data_t> coords(current_cell, m_dx, m_center);

    const auto energy_flux =
        compute_energy_flux(matter_vars, emtensor, h_UU, coords);

    current_cell.store_vars(energy_flux, m_var_enum);
}

template <class matter_t>
template <class data_t>
data_t MatterEnergyFlux<matter_t>::compute_energy_flux(
    const Vars<data_t> &vars, const emtensor_t<data_t> &emtensor,
    const Tensor<2, data_t> &h_UU, const Coordinates<data_t> &coords) const
{
    Tensor<1, data_t> flux_vector_cart_U = 0.0;

    // probably not necessary as this will be done at the end of every timestep
    // but no harm in it I guess...
    data_t chi_regularised = simd_max(vars.chi, 1.0e-6);

    FOR1(i)
    {
        FOR1(j)
        {
            flux_vector_cart_U[i] +=
                vars.lapse * chi_regularised * h_UU[i][j] * emtensor.Si[j];
        }
        flux_vector_cart_U[i] += -vars.shift[i] * emtensor.rho;
    }

    using namespace CoordinateTransformations;

    // only need the r component as this is contracted with the normal to sphere
    data_t flux_vector_spher_U_r = cartesian_to_spherical_U(
        flux_vector_cart_U, coords.x, coords.y, coords.z)[0];

    Tensor<2, data_t> h_phys;
    FOR2(i, j) { h_phys[i][j] = vars.h[i][j] / chi_regularised; }
    Tensor<2, data_t> h_phys_spher =
        cartesian_to_spherical_LL(h_phys, coords.x, coords.y, coords.z);
    data_t area_element = area_element_sphere(h_phys_spher);

    // divide by r^2 sin(theta) as this will be included by the
    // SphericalExtraction integration
    data_t rho2 = simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
    data_t r2sintheta = coords.get_radius() * sqrt(rho2);
    area_element /= r2sintheta;

    // the r component of the unit normal (lowered index) to the spheres.
    // since we are considering surfaces of constant r, the other components
    // vanish
    Tensor<2, data_t> h_phys_spher_UU =
        TensorAlgebra::compute_inverse_sym(h_phys_spher);
    data_t spherical_normal_r = 1.0 / sqrt(h_phys_spher_UU[0][0]);

    return area_element * spherical_normal_r * flux_vector_spher_U_r;
}

#endif /* MATTERENERGYFLUX__IMPL_HPP */
