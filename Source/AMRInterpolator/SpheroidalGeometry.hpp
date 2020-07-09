/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPHEROIDALGEOMETRY_HPP_
#define SPHEROIDALGEOMETRY_HPP_

#include <array>
#include <cmath>

//! This SurfaceGeometry template class provides spheroidal shell geometry
//! implementation for the SurfaceExtraction class, with the symmetry axis
//! fixed to be in the z direction
//! u = t, v = phi
//! note that t (called the eccentric anomaly in astronomy) is NOT the angle of
//! the point with the z-axis, but has a geometric meaning due to Philippe de La
//! Hire (see wikipedia). It is however, measured relative to the z axis.
class SpheroidalGeometry
{
  private:
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_zaxis_over_xaxis; // ratio of major/minor axes

  public:
    SpheroidalGeometry(const std::array<double, CH_SPACEDIM> &a_center,
                       const double a_zaxis_over_xaxis)
        : m_center(a_center), m_zaxis_over_xaxis(a_zaxis_over_xaxis)
    {
    }

    //! returns the grid spacing in t
    inline double du(int a_num_points_t) const
    {
        return M_PI / (double)(a_num_points_t - 1);
    }

    //! returns the grid spacing in phi
    inline double dv(int a_num_points_phi) const
    {
        return 2.0 * M_PI / ((double)a_num_points_phi);
    }

    //! returns the t coordinate associated to the t/u index
    inline double u(int a_it, int a_num_points_t) const
    {
        return a_it * du(a_num_points_t);
    }

    //! returns the phi coordinate associated to the phi/v index
    inline double v(int a_iphi, int a_num_points_phi) const
    {
        return a_iphi * dv(a_num_points_phi);
    }

    inline bool is_u_periodic() const { return false; }
    inline bool is_v_periodic() const { return true; }

    //! returns the Cartesian coordinate in direction a_dir with specified
    //! z_axis magnitude, t and phi.
    inline double get_grid_coord(int a_dir, double a_zaxis, double a_t,
                                 double a_phi) const
    {
        switch (a_dir)
        {
        case (0):
            return m_center[0] +
                   a_zaxis / m_zaxis_over_xaxis * sin(a_t) * cos(a_phi);
        case (1):
            return m_center[1] +
                   a_zaxis / m_zaxis_over_xaxis * sin(a_t) * sin(a_phi);
        case (2):
            return m_center[2] + a_zaxis * cos(a_t);
        default:
            MayDay::Error("SpheroidalGeometry: Direction not supported");
        }
    }

    //! returns the area element on a sphere with radius a_radius at the point
    //! (a_t, a_phi)
    inline double area_element(double a_zaxis, double a_t, double a_phi) const
    {
        double a_z = a_zaxis;
        double a_x = a_zaxis / m_zaxis_over_xaxis;
        double ds_dt = sqrt(a_x * a_x * cos(a_t) * cos(a_t) +
                            a_z * a_z * sin(a_t) * sin(a_t));
        double ds_dphi = a_x * sin(a_t);
        return ds_dt * ds_dphi;
    }

    inline std::string param_name() const { return "r"; }

    inline std::string u_name() const { return "t"; }

    inline std::string v_name() const { return "phi"; }
};

#endif /* SPHEROIDALGEOMETRY_HPP_ */
