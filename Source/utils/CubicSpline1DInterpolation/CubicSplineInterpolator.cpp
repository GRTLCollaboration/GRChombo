/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "parstream.H" //Gives us pout()

#include "CubicSplineInterpolator.hpp"
#include "TriDiagonalMatrix.hpp"

#include <cassert>

// Chombo namespace
#include "UsingNamespace.H"

CubicSplineInterpolator::CubicSplineInterpolator(const std::vector<double> &a_x,
                                                 const std::vector<double> &a_y)
    : solved(false), points_set(false), m_bc_left(second_deriv),
      m_bc_right(second_deriv), m_bc_left_value(0.), m_bc_right_value(0.)
{
    allow_extrapolation(false);
    set_points(a_x, a_y);
}

void CubicSplineInterpolator::allow_extrapolation(bool allow)
{
    m_allow_extrapolation = allow;
}

void CubicSplineInterpolator::set_points(const std::vector<double> &a_x,
                                         const std::vector<double> &a_y,
                                         bool do_solve)
{
    if (a_x.size() > 0)
    {
        assert(a_x.size() == a_y.size());
        m_x = a_x;
        m_y = a_y;
        points_set = true;
        if (do_solve)
            solve();
    }
}

void CubicSplineInterpolator::set_boundary_conditions(BoundaryType a_left,
                                                      double a_left_value,
                                                      BoundaryType a_right,
                                                      double a_right_value,
                                                      bool solve_if_points_set)
{
    m_bc_left = a_left;
    m_bc_right = a_right;
    m_bc_left_value = a_left_value;
    m_bc_right_value = a_right_value;

    if (points_set && solve_if_points_set)
        solve();
}

void CubicSplineInterpolator::solve()
{
    unsigned n = size();

    TriDiagonalMatrix mat(n);
    std::vector<double> rhs(n);

    for (unsigned i = 1; i < n - 1; i++)
    {
        mat.upper(i) = m_x[i] - m_x[i + 1];
        mat.diag(i) = 2. * (m_x[i - 1] - m_x[i + 1]);
        mat.lower(i - 1) = m_x[i - 1] - m_x[i];
        rhs[i] = 6. * ((m_y[i - 1] - m_y[i]) / (m_x[i - 1] - m_x[i]) -
                       (m_y[i] - m_y[i + 1]) / (m_x[i] - m_x[i + 1]));
    }

    // left boundary conditions
    if (m_bc_left == second_deriv)
    {
        // K[0] = f''
        mat.diag(0) = 1.;
        mat.upper(0) = 0.;
        rhs[0] = m_bc_left_value;
    }
    else if (m_bc_left == first_deriv)
    {
        // f', needs to be re-expressed in terms of K:
        // (2K[0]+K[1])(x[0]-x[1]) = 6 ( f' - (y[0]-y[1])/(x[0]-x[1]) )
        mat.diag(0) = 2. * (m_x[0] - m_x[1]);
        mat.upper(0) = (m_x[0] - m_x[1]);
        rhs[0] = 6. * (m_bc_left_value - (m_y[0] - m_y[1]) / (m_x[0] - m_x[1]));
    }
    else
    {
        assert(false);
    }

    // right boundary conditions
    if (m_bc_right == second_deriv)
    {
        // K[n-1] = f''
        mat.diag(n - 1) = 1.;
        mat.lower(n - 2) = 0.;
        rhs[n - 1] = m_bc_right_value;
    }
    else if (m_bc_right == first_deriv)
    {
        // f', needs to be re-expressed in terms of K:
        // (K[n-2]+2K[n-1])(x[n-2]-x[n-1]) =
        // = 6 ( (y[n-2]-y[n-1])/(x[n-2]-x[n-1]) - f' )
        mat.diag(n - 1) = 2. * (m_x[n - 2] - m_x[n - 1]);
        mat.lower(n - 2) = (m_x[n - 2] - m_x[n - 1]);
        rhs[0] = 6. * ((m_y[n - 2] - m_y[n - 1]) / (m_x[n - 2] - m_x[n - 1]) -
                       m_bc_right_value);
    }
    else
    {
        assert(false);
    }

    m_K = mat.lu_solve(rhs);
    solved = true;
}

double CubicSplineInterpolator::interpolate(double x, unsigned derivative) const
{
    assert(solved);

    unsigned n = size();

    unsigned i;
    for (i = 0; i < n; ++i)
        if (x < m_x[i])
            break;

    if (i == 0 || i == n)
    {
        std::string msg =
            "CubicSplineInterpolator::interpolate - extrapolating at " +
            std::to_string(x) + " outside of valid range [" +
            std::to_string(m_x[0]) + "," + std::to_string(m_x[n - 1]) + "]";
        pout() << msg << std::endl;
        assert(m_allow_extrapolation);

        if (i == 0) // use segment 0
            ++i;
        else // use segment n-1
            --i;
    }

    double d_left = (x - m_x[i - 1]);
    double d_right = (x - m_x[i]);
    double delta = (m_x[i - 1] - m_x[i]);

    double answer;
    switch (derivative)
    {
    case 0:
        answer =
            m_K[i - 1] / 6. *
                (d_right * d_right * d_right / delta - d_right * delta) -
            m_K[i] / 6. * (d_left * d_left * d_left / delta - d_left * delta) +
            (m_y[i - 1] * d_right - m_y[i] * d_left) / delta;
        break;
    case 1:
        answer = m_K[i - 1] / 6. * (3. * d_right * d_right / delta - delta) -
                 m_K[i] / 6. * (3 * d_left * d_left / delta - delta) +
                 (m_y[i - 1] - m_y[i]) / delta;
        break;
    case 2:
        answer = (m_K[i - 1] * d_right - m_K[i] * d_left) / delta;
        break;
    case 3:
        answer = (m_K[i - 1] - m_K[i]) / delta;
        break;
    default: // all others
        answer = 0.;
    }

    return answer;
}
