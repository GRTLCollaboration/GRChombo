/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPHERICALGEOMETRY_HPP_
#define SPHERICALGEOMETRY_HPP_

#include "SurfaceGeometry.hpp"
#include <array>
#include <cmath>

//! This class, derived from the SurfaceGeometry class, provides spherical shell
//! geometry implementation for the SurfaceExtraction class
//! u = theta, v = phi
class SphericalGeometry : public SurfaceGeometry
{
  private:
    std::array<double, CH_SPACEDIM> m_center;

  public:
    SphericalGeometry(std::array<double, CH_SPACEDIM> &a_center)
        : m_center(a_center)
    {
    }

    //! returns the grid spacing in theta
    virtual double du(int a_num_points_theta) const override
    {
        return M_PI / ((double)a_num_points_theta);
    }

    //! returns the grid spacing in phi
    virtual double dv(int a_num_points_phi) const override
    {
        return 2.0 * M_PI / ((double)a_num_points_phi);
    }

    //! returns the theta coordinate associated to the theta/u index
    virtual double u(int a_itheta, int a_num_points_theta) const override
    {
        return (a_itheta + 0.5) * du(a_num_points_theta);
    }

    //! returns the phi coordinate associated to the phi/v index
    virtual double v(int a_iphi, int a_num_points_phi) const override
    {
        return (a_iphi)*dv(a_num_points_phi);
    }

    //! returns the Cartesian coordinate in direction a_dir with specified
    //! radius, theta and phi.
    virtual double cartesian_coord(int a_dir, double a_radius, double a_theta,
                                   double a_phi) const override
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
    virtual double area_element(double a_radius, double a_theta,
                                double a_phi) const override
    {
        return a_radius * a_radius * sin(a_theta);
    }

    virtual std::string param_name() const override { return "r"; }

    virtual std::string u_name() const override { return "theta"; }

    virtual std::string v_name() const override { return "phi"; }
};

#endif /* SPHERICALGEOMETRY_HPP_ */
