/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef AHSPHERICALGEOMETRYUNIFORM_HPP_
#define AHSPHERICALGEOMETRYUNIFORM_HPP_

#include "DimensionDefinitions.hpp" // make sure GR_SPACEDIM exists

#if GR_SPACEDIM != 3
#error "This file should only be included for GR_SPACEDIM == 3."
#endif

// Chombo includes
#include "CH_Timer.H"

// Other includes
#include "AHGeometryData.hpp"
#include "AlwaysInline.hpp"
#include "SphericalGeometryUniform.hpp"
#include <array>
#include <cmath>

// Chombo namespace
#include "UsingNamespace.H"

//! This SurfaceGeometry template class provides spherical shell geometry
//! implementation for the SurfaceExtraction class
//! u = theta, v = phi
class AHSphericalGeometryUniform : public SphericalGeometryUniform
{
  public:
    AHSphericalGeometryUniform(const std::array<double, CH_SPACEDIM> &a_origin)
        : SphericalGeometryUniform(a_origin)
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

    AHGeometryData get_geometry_data(double r, double u, double phi = 0.) const
    {
        // for 2D Cartoon spherical geometry (3D->2D), this should always be
        // called with no 3rd argument (so phi=0), which means the formulas
        // below are effectively polar coordinates

        CH_TIME("AHSphericalGeometryUniform::get_geometry_data");

        double sqrt1mu2 = sqrt(1. - u * u);

        double cosphi = cos(phi);
        double sinphi = sin(phi);

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

        out.du[ix] = -u * sqrt1mu2 * cosphi / r;
#if CH_SPACEDIM == 3
        out.du[iy] = -u * sqrt1mu2 * sinphi / r;
#endif
        out.du[iz] = (1. - u * u) / r;

#if CH_SPACEDIM == 3
        out.dv[ix] = (-sinphi) / (r * sqrt1mu2);
        out.dv[iy] = (cosphi) / (r * sqrt1mu2);
        out.dv[iz] = 0;
#endif

        double dfdx = sqrt1mu2 * cosphi;
#if CH_SPACEDIM == 3
        double dfdy = sqrt1mu2 * sinphi;
#endif
        double dfdz = u;

        out.ddu[ix][ix] = -u *
                          (-1. + 3. * u * u + 3. * (-1. + u * u) * cos2phi) /
                          (r * r * 2.);
        out.ddu[ix][iz] = sqrt1mu2 * (-1. + 3. * u * u) * cosphi / (r * r);
        out.ddu[iz][ix] = out.ddu[ix][iz];
        out.ddu[iz][iz] = 3. * u * (-1. + u * u) / (r * r);
#if CH_SPACEDIM == 3
        out.ddu[ix][iy] = -3. * u * (-1. + u * u) * cosphi * sinphi / (r * r);
        out.ddu[iy][iy] =
            u * (1. - 3. * u * u + 3. * (-1. + u * u) * cos2phi) / (2. * r * r);
        out.ddu[iy][iz] = sqrt1mu2 * (-1. + 3. * u * u) * sinphi / (r * r);
        out.ddu[iy][ix] = out.ddu[ix][iy];
        out.ddu[iz][iy] = out.ddu[iy][iz];
#endif

#if CH_SPACEDIM == 3
        out.ddv[ix][ix] = -sin2phi / (r * r * (-1. + u * u));
        out.ddv[ix][iz] = 0.;
        out.ddv[iz][ix] = out.ddv[ix][iz];
        out.ddv[iz][iz] = 0.;
        out.ddv[ix][iy] = cos2phi / (r * r * (-1. + u * u));
        out.ddv[iy][iy] = sin2phi / (r * r * (-1. + u * u));
        out.ddv[iy][iz] = 0.;
        out.ddv[iy][ix] = out.ddv[ix][iy];
        out.ddv[iz][iy] = out.ddv[iy][iz];
#endif

        double ddfdxdx = (u * u - (-1. + u * u) * sinphi * sinphi) / r;
        double ddfdxdz = -(u * sqrt1mu2 * cosphi) / r;
        double ddfdzdz = (1. - u * u) / r;
#if CH_SPACEDIM == 3
        double ddfdxdy = (-1. + u * u) * cosphi * sinphi / r;
        double ddfdydy = (u * u - (-1. + u * u) * cosphi * cosphi) / r;
        double ddfdydz = -(u * sqrt1mu2 * sinphi) / r;
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

        out.dxdu[ix] = -r * u / sqrt1mu2 * cosphi;
#if CH_SPACEDIM == 3
        out.dxdu[iy] = -r * u / sqrt1mu2 * sinphi;
#endif
        out.dxdu[iz] = r;

#if CH_SPACEDIM == 3
        out.dxdv[ix] = -r * sqrt1mu2 * sinphi;
        out.dxdv[iy] = r * sqrt1mu2 * cosphi;
        out.dxdv[iz] = 0.;
#endif

        out.dxdf[ix] = sqrt1mu2 * cosphi;
#if CH_SPACEDIM == 3
        out.dxdf[iy] = sqrt1mu2 * sinphi;
#endif
        out.dxdf[iz] = u;

        return out;
    }
};

#endif /* AHSPHERICALGEOMETRYUNIFORM_HPP_ */
