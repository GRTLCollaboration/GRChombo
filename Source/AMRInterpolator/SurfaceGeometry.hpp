/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SURFACEGEOMETRY_HPP_
#define SURFACEGEOMETRY_HPP_

#include <string>

//! An abstract base class to provide required methods for geometries in the
//! SurfaceExtraction class
class SurfaceGeometry
{
  public:
    //! returns the grid spacing in u
    virtual double du(int a_num_points_u) const = 0;

    //! returns the grid spacing in v
    virtual double dv(int a_num_points_v) const = 0;

    //! returns the u coordinate associated to the u index
    virtual double u(int a_iu, int a_num_points_u) const = 0;

    //! returns the v coordinate associated to the v index
    virtual double v(int a_iv, int a_num_points_v) const = 0;

    //! returns the periodicity of the u coordinate
    virtual bool is_u_periodic() const = 0;

    //! returns the periodicity of the v coordinate
    virtual bool is_v_periodic() const = 0;

    //! returns the Cartesian coordinate in direction a_dir of the point on the
    //! surface with specified param value with given surface coordinates
    virtual double cartesian_coord(int a_dir, double a_surface_param_value,
                                   double a_u, double a_v) const = 0;

    //! returns the Cartesian coordinate in direction a_dir of the point on the
    //! surface with specified param value with given surface coordinate indices
    double cartesian_coord(int a_dir, double a_surface_param_value, int a_iu,
                           int a_num_points_u, int a_iv,
                           int a_num_points_v) const
    {
        return cartesian_coord(a_dir, a_surface_param_value,
                               u(a_iu, a_num_points_u),
                               v(a_iv, a_num_points_v));
    }

    //! returns the area element of the surface at a point
    virtual double area_element(double a_radius, double a_theta,
                                double a_phi) const = 0;

    //! returns the name associated to the parameter that labels the surface
    virtual std::string param_name() const { return "surface parameter"; }

    //! returns the name associated to the u coordinate
    virtual std::string u_name() const { return "u"; }

    //! returns the name associated to the v coordinate
    virtual std::string v_name() const { return "v"; }
};

#endif /* SURFACEGEOMETRY_HPP_ */
