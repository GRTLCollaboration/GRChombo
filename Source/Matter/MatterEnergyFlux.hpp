/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MATTERENERGYFLUX_HPP
#define MATTERENERGYFLUX_HPP

#include "Cell.hpp"
#include "CoordinateTransformations.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4.hpp"
#include "simd.hpp"

//! This compute class computes the integrand for the outgoing matter energy
//! flux on spheres
template <class matter_t> class MatterEnergyFlux
{
  public:
    template <class data_t>
    using Vars = typename MatterCCZ4<matter_t>::template Vars<data_t>;

    //! Constructor
    MatterEnergyFlux(matter_t a_matter,
                     const std::array<double, CH_SPACEDIM> &a_center,
                     const double a_dx, const int a_var_enum);

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    template <class data_t>
    data_t compute_energy_flux(const Vars<data_t> &vars,
                               const emtensor_t<data_t> &emtensor,
                               const Tensor<2, data_t> &h_UU,
                               const Coordinates<data_t> &coords) const;

    matter_t m_matter;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const int m_var_enum;
};

#include "MatterEnergyFlux.impl.hpp"

#endif /* MATTERENERGYFLUX_HPP */
