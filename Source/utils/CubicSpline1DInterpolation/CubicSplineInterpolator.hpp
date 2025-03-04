/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CUBICSPLINEINTERPOLATOR_HPP_
#define CUBICSPLINEINTERPOLATOR_HPP_

#include <string>
#include <vector>

//! Calculate 'x's, 'y's and 'K's in Mathematica, import them to here
//! and this class interpolates with cubic spline
class CubicSplineInterpolator
{
  public:
    enum BoundaryType
    {
        first_deriv = 1,
        second_deriv
    };

    CubicSplineInterpolator(
        const std::vector<double> &a_x = std::vector<double>(),
        const std::vector<double> &a_y = std::vector<double>());

    void allow_extrapolation(bool allow);
    void set_points(const std::vector<double> &a_x,
                    const std::vector<double> &a_y, bool do_solve = true);
    // optional, default is 2nd derivative = 0
    void set_boundary_conditions(BoundaryType a_left, double a_left_value,
                                 BoundaryType a_right, double a_right_value,
                                 bool solve_if_points_set = true);

    double interpolate(double x, unsigned derivative = 0) const;

    void solve();

    inline unsigned size() const { return m_x.size(); }
    inline std::vector<double> &x() { return m_x; }
    inline std::vector<double> x() const { return m_x; }
    inline std::vector<double> &y() { return m_y; }
    inline std::vector<double> y() const { return m_y; }

    // return 2nd derivatives
    inline std::vector<double> &K() { return m_K; }
    inline std::vector<double> K() const { return m_K; }

  private:
    std::vector<double> m_x, m_y;
    bool m_allow_extrapolation;

    bool solved, points_set;
    BoundaryType m_bc_left, m_bc_right;
    double m_bc_left_value, m_bc_right_value;

    std::vector<double> m_K; // storing 2nd derivatives on each point
};

#endif /*  CUBICSPLINEINTERPOLATOR_HPP_ */
