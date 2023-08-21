/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef AHSPHERICALGEOMETRY_HPP_
#define AHSPHERICALGEOMETRY_HPP_

#include "DimensionDefinitions.hpp" // make sure GR_SPACEDIM exists

#if GR_SPACEDIM != 3
#error "This file should only be included for GR_SPACEDIM == 3."
#endif

// Chombo includes
#include "CH_Timer.H"

// Other includes
#include "AHGeometryData.hpp"
#include "AlwaysInline.hpp"
#include "SphericalGeometry.hpp"
#include <array>
#include <cmath>

// Chombo namespace
#include "UsingNamespace.H"

//! This SurfaceGeometry template class provides spherical shell geometry
//! implementation for the SurfaceExtraction class
//! u = theta, v = phi
class AHSphericalGeometry : public SphericalGeometry
{
  public:
    AHSphericalGeometry(const std::array<double, CH_SPACEDIM> &a_origin)
        : SphericalGeometry(a_origin)
    {
    }

    ALWAYS_INLINE void
    set_origin(const std::array<double, CH_SPACEDIM> &a_origin)
    {
        m_center = a_origin;
    }
    ALWAYS_INLINE const std::array<double, CH_SPACEDIM> &get_origin() const
    {
        return m_center;
    }

    AHGeometryData get_geometry_data(double r, double theta,
                                     double phi = 0.) const
    {
        // for 2D Cartoon spherical geometry (3D->2D), this should always be
        // called with no 3rd argument (so phi=0), which means the formulas
        // below are effectively polar coordinates

        CH_TIME("AHSphericalGeometry::get_geometry_data");

        double costheta = cos(theta);
        double sintheta = sin(theta);
        double tantheta = tan(theta);

        double cosphi = cos(phi);
        double sinphi = sin(phi);

        double cos2theta = cos(2. * theta);
        double sin2theta = sin(2. * theta);
        double cos2phi = cos(2. * phi);
        double sin2phi = sin(2. * phi);

        // in 3D, (x,y,z) are normal
        // in 2D Cartoon method, the spherical 'z' axis is the grid 'x' axis, so
        // the 2D coordinates (x,y) are (z,x) spherical coordinates
        // (that means commenting out everything related to 'y' and to 'v'
        // [=phi])
        int iz = m_up_dir;
        int ix = (m_up_dir + 1) % 3;
#if CH_SPACEDIM == 3
        int iy = (m_up_dir + 2) % 3;
#endif

        AHGeometryData out;

        out.du[ix] = (costheta * cosphi) / r;
#if CH_SPACEDIM == 3
        out.du[iy] = (costheta * sinphi) / r;
#endif
        out.du[iz] = -sintheta / r;

#if CH_SPACEDIM == 3
        out.dv[ix] = (-sinphi) / (r * sintheta);
        out.dv[iy] = (cosphi) / (r * sintheta);
        out.dv[iz] = 0;
#endif

        double dfdx = sintheta * cosphi;
#if CH_SPACEDIM == 3
        double dfdy = sintheta * sinphi;
#endif
        double dfdz = costheta;

        out.ddu[ix][ix] =
            (cos2theta * cosphi * cosphi - cos2phi) / (r * r * tantheta);
        out.ddu[ix][iz] = -(cos2theta * cosphi) / (r * r);
        out.ddu[iz][ix] = out.ddu[ix][iz];
        out.ddu[iz][iz] = (sin2theta) / (r * r);
#if CH_SPACEDIM == 3
        out.ddu[ix][iy] =
            ((cos2theta - 2.) * sin2phi) / (2. * r * r * tantheta);
        out.ddu[iy][iy] =
            (cos2theta * sinphi * sinphi + cos2phi) / (r * r * tantheta);
        out.ddu[iy][iz] = -(cos2theta * sinphi) / (r * r);
        out.ddu[iy][ix] = out.ddu[ix][iy];
        out.ddu[iz][iy] = out.ddu[iy][iz];
#endif

#if CH_SPACEDIM == 3
        out.ddv[ix][ix] = sin2phi / (r * r * sintheta * sintheta);
        out.ddv[ix][iz] = 0.;
        out.ddv[iz][ix] = out.ddv[ix][iz];
        out.ddv[iz][iz] = 0.;
        out.ddv[ix][iy] = -cos2phi / (r * r * sintheta * sintheta);
        out.ddv[iy][iy] = -sin2phi / (r * r * sintheta * sintheta);
        out.ddv[iy][iz] = 0.;
        out.ddv[iy][ix] = out.ddv[ix][iy];
        out.ddv[iz][iy] = out.ddv[iy][iz];
#endif

        double ddfdxdx =
            (costheta * costheta + sintheta * sintheta * sinphi * sinphi) / r;
        double ddfdxdz = -(sintheta * costheta * cosphi) / r;
        double ddfdzdz = (sintheta * sintheta) / r;
#if CH_SPACEDIM == 3
        double ddfdxdy = -(sintheta * sintheta * sinphi * cosphi) / r;
        double ddfdydy =
            (costheta * costheta + sintheta * sintheta * cosphi * cosphi) / r;
        double ddfdydz = -(sintheta * costheta * sinphi) / r;
#endif

        out.df[ix] = dfdx;
#if CH_SPACEDIM == 3
        out.df[iy] = dfdy;
#endif
        out.df[iz] = dfdz;

        out.ddf[ix][ix] = ddfdxdx;
        out.ddf[ix][iz] = ddfdxdz;
        out.ddf[iz][ix] = out.ddf[ix][iz];
        out.ddf[iz][iz] = ddfdzdz;
#if CH_SPACEDIM == 3
        out.ddf[ix][iy] = ddfdxdy;
        out.ddf[iy][iy] = ddfdydy;
        out.ddf[iy][iz] = ddfdydz;
        out.ddf[iy][ix] = out.ddf[ix][iy];
        out.ddf[iz][iy] = out.ddf[iy][iz];
#endif

        out.dxdu[ix] = r * costheta * cosphi;
#if CH_SPACEDIM == 3
        out.dxdu[iy] = r * costheta * sinphi;
#endif
        out.dxdu[iz] = -r * sintheta;

#if CH_SPACEDIM == 3
        out.dxdv[ix] = -r * sintheta * sinphi;
        out.dxdv[iy] = r * sintheta * cosphi;
        out.dxdv[iz] = 0.;
#endif

        out.dxdf[ix] = sintheta * cosphi;
#if CH_SPACEDIM == 3
        out.dxdf[iy] = sintheta * sinphi;
#endif
        out.dxdf[iz] = costheta;

        return out;
    }
};

#endif /* AHSPHERICALGEOMETRY_HPP_ */
