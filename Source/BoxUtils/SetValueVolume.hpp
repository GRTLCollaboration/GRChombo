/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETVALUEVOLUME_HPP_
#define SETVALUEVOLUME_HPP_

#include "Cell.hpp"
#include "Constraints.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Interval.H"
#include "simd.hpp"

/// This volume corresponds to the inside of a ball of specified radius and
/// center
class SphericalVolume
{
  private:
    const double m_radius;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_dx;

  public:
    //! Constructor
    SphericalVolume(const double a_radius,
                    const std::array<double, CH_SPACEDIM> &a_center,
                    const double a_dx)
        : m_radius(a_radius), m_center(a_center), m_dx(a_dx)
    {
    }

    template <class data_t>
    auto in_volume(const Cell<data_t> &current_cell) const
    {
        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        data_t r = coords.get_radius();
        return simd_compare_lt(r, m_radius);
    }
};

/// This volume corresponds to the complement of a ball (i.e. the whole grid
/// excluding a ball) with specified radius and center
class SphericalVolumeComplement
{
  private:
    const double m_radius;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_dx;

  public:
    //! Constructor
    SphericalVolumeComplement(const double a_radius,
                              const std::array<double, CH_SPACEDIM> &a_center,
                              const double a_dx)
        : m_radius(a_radius), m_center(a_center), m_dx(a_dx)
    {
    }

    template <class data_t>
    auto in_volume(const Cell<data_t> &current_cell) const
    {
        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        data_t r = coords.get_radius();
        return simd_compare_gt(r, m_radius);
    }
};

/// This compute class can set the value of an interval of variables to a
/// specified value within a volume
template <class volume_t> class SetValueVolume
{
  private:
    const double m_value;
    const Interval m_interval;
    const volume_t m_volume;

  public:
    //! constructor
    SetValueVolume(const double a_value, const Interval a_interval,
                   const volume_t a_volume)
        : m_value(a_value), m_interval(a_interval), m_volume(a_volume)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // ask the m_volume object if we're in the volume
        auto in_volume = m_volume.in_volume(current_cell);

        // set var to m_value if in_volume is true, otherwise just the existing
        // value
        for (int ivar = m_interval.begin(); ivar <= m_interval.end(); ++ivar)
        {
            auto in_val = current_cell.load_vars(ivar);
            data_t out_val = simd_conditional(in_volume, m_value, in_val);
            current_cell.store_vars(out_val, ivar);
        }
    }
};

#endif /* SETVALUEVOLUME_HPP_ */
