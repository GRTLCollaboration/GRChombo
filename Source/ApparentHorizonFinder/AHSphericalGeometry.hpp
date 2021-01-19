/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef AHSPHERICALGEOMETRY_HPP_
#define AHSPHERICALGEOMETRY_HPP_

#if CH_SPACEDIM != 3
#error "This file should only be included for 3+1D simulations"
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

    AHGeometryData get_geometry_data(double r, double theta, double phi) const
    {
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

        AHGeometryData out;

        out.du[0] = (costheta * cosphi) / r;
        out.du[1] = (costheta * sinphi) / r;
        out.du[2] = -sintheta / r;

        out.dv[0] = (-sinphi) / (r * sintheta);
        out.dv[1] = (cosphi) / (r * sintheta);
        out.dv[2] = 0;

        double dfdx = sintheta * cosphi;
        double dfdy = sintheta * sinphi;
        double dfdz = costheta;

        out.ddu[0][0] =
            (cos2theta * cosphi * cosphi - cos2phi) / (r * r * tantheta);
        out.ddu[0][1] = ((cos2theta - 2.) * sin2phi) / (2. * r * r * tantheta);
        out.ddu[0][2] = -(cos2theta * cosphi) / (r * r);
        out.ddu[1][1] =
            (cos2theta * sinphi * sinphi + cos2phi) / (r * r * tantheta);
        out.ddu[1][2] = -(cos2theta * sinphi) / (r * r);
        out.ddu[2][2] = (sin2theta) / (r * r);
        out.ddu[1][0] = out.ddu[0][1];
        out.ddu[2][0] = out.ddu[0][2];
        out.ddu[2][1] = out.ddu[1][2];

        out.ddv[0][0] = sin2phi / (r * r * sintheta * sintheta);
        out.ddv[0][1] = -cos2phi / (r * r * sintheta * sintheta);
        out.ddv[0][2] = 0.;
        out.ddv[1][1] = -sin2phi / (r * r * sintheta * sintheta);
        out.ddv[1][2] = 0.;
        out.ddv[2][2] = 0.;
        out.ddv[1][0] = out.ddv[0][1];
        out.ddv[2][0] = out.ddv[0][2];
        out.ddv[2][1] = out.ddv[1][2];

        double ddfdxdx =
            (costheta * costheta + sintheta * sintheta * sinphi * sinphi) / r;
        double ddfdxdy = -(sintheta * sintheta * sinphi * cosphi) / r;
        double ddfdxdz = -(sintheta * costheta * cosphi) / r;
        double ddfdydy =
            (costheta * costheta + sintheta * sintheta * cosphi * cosphi) / r;
        double ddfdydz = -(sintheta * costheta * sinphi) / r;
        double ddfdzdz = (sintheta * sintheta) / r;

        out.df[0] = dfdx;
        out.df[1] = dfdy;
        out.df[2] = dfdz;

        out.ddf[0][0] = ddfdxdx;
        out.ddf[0][1] = ddfdxdy;
        out.ddf[0][2] = ddfdxdz;
        out.ddf[1][1] = ddfdydy;
        out.ddf[1][2] = ddfdydz;
        out.ddf[2][2] = ddfdzdz;
        out.ddf[1][0] = out.ddf[0][1];
        out.ddf[2][0] = out.ddf[0][2];
        out.ddf[2][1] = out.ddf[1][2];

        out.dxdu[0] = r * costheta * cosphi;
        out.dxdu[1] = r * costheta * sinphi;
        out.dxdu[2] = -r * sintheta;

        out.dxdv[0] = -r * sintheta * sinphi;
        out.dxdv[1] = r * sintheta * cosphi;
        out.dxdv[2] = 0.;

        out.dxdf[0] = sintheta * cosphi;
        out.dxdf[1] = sintheta * sinphi;
        out.dxdf[2] = costheta;

        return out;
    }
};

#endif /* AHSPHERICALGEOMETRY_HPP_ */
