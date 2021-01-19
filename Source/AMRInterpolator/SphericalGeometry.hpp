/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPHERICALGEOMETRY_HPP_
#define SPHERICALGEOMETRY_HPP_

#if CH_SPACEDIM != 3
#error "This file should only be included for 3+1D simulations"
#endif

// Chombo includes
#include "AlwaysInline.hpp"
#include "MayDay.H"

// Other includes
#include "AlwaysInline.hpp"
#include "MayDay.H"
#include <array>
#include <cmath>
#include <string>

// Chombo namespace
#include "UsingNamespace.H"

//! This SurfaceGeometry template class provides spherical shell geometry
//! implementation for the SurfaceExtraction class
//! u = theta, v = phi
class SphericalGeometry
{
  protected:
    std::array<double, CH_SPACEDIM> m_center;

  public:
    SphericalGeometry(const std::array<double, CH_SPACEDIM> &a_center)
        : m_center(a_center)
    {
    }

    ALWAYS_INLINE double get_domain_u_min() const { return 0.; }
    ALWAYS_INLINE double get_domain_u_max() const { return M_PI; }
    ALWAYS_INLINE double get_domain_v_min() const { return 0.; }
    ALWAYS_INLINE double get_domain_v_max() const { return 2. * M_PI; }

    //! returns the grid spacing in theta
    ALWAYS_INLINE double du(int a_num_points_theta) const
    {
        return (get_domain_u_max() - get_domain_u_min()) /
               (double)(a_num_points_theta - 1);
    }

    //! returns the grid spacing in phi
    ALWAYS_INLINE double dv(int a_num_points_phi) const
    {
        return (get_domain_v_max() - get_domain_v_min()) /
               ((double)a_num_points_phi);
    }

    //! returns the theta coordinate associated to the theta/u index
    ALWAYS_INLINE double u(int a_itheta, int a_num_points_theta) const
    {
        return a_itheta * du(a_num_points_theta);
    }

    //! returns the phi coordinate associated to the phi/v index
    ALWAYS_INLINE double v(int a_iphi, int a_num_points_phi) const
    {
        return a_iphi * dv(a_num_points_phi);
    }

    ALWAYS_INLINE bool is_u_periodic() const { return false; }
    ALWAYS_INLINE bool is_v_periodic() const { return true; }

    ALWAYS_INLINE std::string param_name() const { return "r"; }

    ALWAYS_INLINE std::string u_name() const { return "theta"; }

    ALWAYS_INLINE std::string v_name() const { return "phi"; }

    //! returns the area element on a sphere with radius a_radius at the point
    //! (a_theta, a_phi)
    ALWAYS_INLINE double area_element(double a_radius, double a_theta,
                                      double a_phi) const
    {
        return a_radius * a_radius * sin(a_theta);
    }

    //! returns the Cartesian coordinate in direction a_dir with specified
    //! radius, theta and phi.
    ALWAYS_INLINE double get_grid_coord(int a_dir, double a_radius,
                                        double a_theta, double a_phi) const
    {
        switch (a_dir)
        {
        case (0):
            return m_center[0] + a_radius * sin(a_theta) * cos(a_phi);
        case (1):
            return m_center[1] + a_radius * sin(a_theta) * sin(a_phi);
        case (2):
            return m_center[2] + a_radius * cos(a_theta);
        default:
            MayDay::Error("SphericalGeometry: Direction not supported");
        }
    }
};

#endif /* SPHERICALGEOMETRY_HPP_ */
