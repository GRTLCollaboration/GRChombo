/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHSPHERICALGEOMETRY_HPP_
#define _AHSPHERICALGEOMETRY_HPP_

#if CH_SPACEDIM != 2
#error "This file should only be included for 2+1D simulations"
#endif

// Chombo includes
#include "CH_Timer.H"
#include "MayDay.H"

// Other includes
#include "AHGeometryData.hpp"
#include "AlwaysInline.hpp"
#include <array>

// Chombo namespace
#include "UsingNamespace.H"

//! Class for coordinate system to use
//! String assumed to be along in 'x' direction
//! Periodic in x coordinate, y direction may be reflective
class AHStringGeometry
{
  private:
    double m_string_length;
    std::array<double, CH_SPACEDIM> fake_origin = {0., 0.};

  public:
    AHStringGeometry(double a_string_length) : m_string_length(a_string_length)
    {
    }

    ALWAYS_INLINE void
    set_origin(const std::array<double, CH_SPACEDIM> &a_origin)
    {
    }
    ALWAYS_INLINE const std::array<double, CH_SPACEDIM> &get_origin() const
    {
        return fake_origin;
    }

    ALWAYS_INLINE double get_domain_u_min() const { return 0.; }
    ALWAYS_INLINE double get_domain_u_max() const { return m_string_length; }

    ALWAYS_INLINE bool is_u_periodic() const { return true; }

    ALWAYS_INLINE std::string param_name() const { return "y"; }

    ALWAYS_INLINE std::string u_name() const { return "x"; }

    ALWAYS_INLINE double get_grid_coord(int a_dir, double a_y, double a_x) const
    {
        switch (a_dir)
        {
        case (0):
            return a_x;
        case (1):
            return a_y;
        default:
            MayDay::Error("AHStringGeometry: Direction not supported");
        }
    }

    AHGeometryData get_geometry_data(double y, double x) const
    {
        CH_TIME("AHStringGeometry::get_geometry_data");

        AHGeometryData out;

        out.du[0] = 1.;
        out.du[1] = 0.;

        double dfdx = 0.;
        double dfdy = 1.;

        out.ddu[0][0] = 0.;
        out.ddu[0][1] = 0.;
        out.ddu[1][1] = 0.;
        out.ddu[1][0] = out.ddu[0][1];

        double ddfdxdx = 0.;
        double ddfdxdy = 0.;
        double ddfdydy = 0.;

        out.df[0] = dfdx;
        out.df[1] = dfdy;

        out.ddf[0][0] = ddfdxdx;
        out.ddf[0][1] = ddfdxdy;
        out.ddf[1][1] = ddfdydy;
        out.ddf[1][0] = out.ddf[0][1];

        out.dxdu[0] = 1.;
        out.dxdu[1] = 0.;

        out.dxdf[0] = 0.;
        out.dxdf[1] = 1.;

        return out;
    }
};

#endif /* _AHSPHERICALGEOMETRY_HPP_ */
