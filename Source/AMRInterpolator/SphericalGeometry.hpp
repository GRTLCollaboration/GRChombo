/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPHERICALGEOMETRY_HPP_
#define SPHERICALGEOMETRY_HPP_

// Chombo includes
#include "MayDay.H"

// Other includes
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
  private:
    std::array<double, CH_SPACEDIM> m_center;

  public:
    SphericalGeometry(const std::array<double, CH_SPACEDIM> &a_center)
        : m_center(a_center)
    {
    }

    //! returns the grid spacing in theta
    inline double du(int a_num_points_theta) const
    {
        return M_PI / (double)(a_num_points_theta - 1);
    }

    //! returns the grid spacing in phi
    inline double dv(int a_num_points_phi) const
    {
        return 2.0 * M_PI / ((double)a_num_points_phi);
    }

    //! returns the theta coordinate associated to the theta/u index
    inline double u(int a_itheta, int a_num_points_theta) const
    {
        return a_itheta * du(a_num_points_theta);
    }

    //! returns the phi coordinate associated to the phi/v index
    inline double v(int a_iphi, int a_num_points_phi) const
    {
        return a_iphi * dv(a_num_points_phi);
    }

    inline bool is_u_periodic() const { return false; }
    inline bool is_v_periodic() const { return true; }

    //! returns the Cartesian coordinate in direction a_dir with specified
    //! radius, theta and phi.
    inline double get_grid_coord(int a_dir, double a_radius, double a_theta,
                                 double a_phi) const
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

    //! returns the area element on a sphere with radius a_radius at the point
    //! (a_theta, a_phi)
    inline double area_element(double a_radius, double a_theta,
                               double a_phi) const
    {
        return a_radius * a_radius * sin(a_theta);
    }

    inline std::string param_name() const { return "r"; }

    inline std::string u_name() const { return "theta"; }

    inline std::string v_name() const { return "phi"; }
};

#endif /* SPHERICALGEOMETRY_HPP_ */
