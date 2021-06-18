/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FILEINTERPOLATOR1D_HPP_
#define FILEINTERPOLATOR1D_HPP_

#include "CubicSplineInterpolator.hpp"
#include "GRParmParse.hpp"

#include <cassert>
#include <string>

//! Calculate 'x's, 'y's elsewhere, import them to here
//! and this class interpolates with 4th order Lagrange
//! (left as a separate class to keep 'CubicSplineInterpolator' free from any
//! Chombo/GRChombo code)
class FileInterpolator1D
{
  public:
    FileInterpolator1D(const std::string &N_label, const std::string &x_label,
                       const std::string &y_label,
                       double y_default_if_missing_file = 0.)
    {
        GRParmParse pp;

        int n = 0;
        std::vector<double> xs, ys;
        if (pp.contains(N_label.c_str()))
        {
            pp.load(N_label.c_str(), n, 0);
            pp.load(x_label.c_str(), xs, n, {0.});
            pp.load(y_label.c_str(), ys, n, {0.});
        }

        assert(n >= 0);
        if (n >= 2)
            spline.set_points(xs, ys);
        else if (n == 1)
        {
            spline.allow_extrapolation(true);
            // *2 such that it avoid numerical erros if 'x' is very big (+1
            // would do little to x=1e20) +1 to still work for x=0
            spline.set_points({xs[0], xs[0] * 2. + 1.}, {ys[0], ys[0]});
        }
        else
        {
            pout()
                << "Parameter: " << N_label
                << " not found in parameter file. FileInterpolator1D will set "
                << y_label
                << " to its default value = " << y_default_if_missing_file
                << " for all " << x_label << std::endl;

            spline.allow_extrapolation(true);
            spline.set_points({-1.e5, 1.e5}, {y_default_if_missing_file,
                                              y_default_if_missing_file});
        }
    }

    inline void allow_extrapolation(bool allow)
    {
        if (spline.size() >= 2)
            spline.allow_extrapolation(allow);
    }

    // optional, default is 2nd derivative = 0
    inline void set_boundary_conditions(
        CubicSplineInterpolator::BoundaryType a_left, double a_left_value,
        CubicSplineInterpolator::BoundaryType a_right, double a_right_value)
    {
        spline.set_boundary_conditions(a_left, a_left_value, a_right,
                                       a_right_value);
    }
    inline double interpolate(double x, int derivative = 0)
    {
        return spline.interpolate(x, derivative);
    }

  private:
    CubicSplineInterpolator spline;
};

#endif /*  FILEINTERPOLATOR1D_HPP_ */
