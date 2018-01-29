/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef QUINTICCONVOLUTION_IMPL_HPP_
#define QUINTICCONVOLUTION_IMPL_HPP_

const string QuinticConvolution::TAG = "\x1b[36;1m[QuinticConvolution]\x1b[0m ";

QuinticConvolution::QuinticConvolution(const InterpSource &source,
                                       bool verbosity)
    : m_source(source), m_verbosity(verbosity)
{
    CH_assert(CH_SPACEDIM <= 3);
}

void QuinticConvolution::setup(const std::array<int, CH_SPACEDIM> &deriv,
                               const std::array<double, CH_SPACEDIM> &dx,
                               const std::array<double, CH_SPACEDIM> &evalCoord,
                               const IntVect &nearest)
{
    m_interp_points.clear();
    m_interp_weights.clear();

    double weights_1d[CH_SPACEDIM][6];

    for (int dim = 0; dim < CH_SPACEDIM; ++dim)
    {
        double s = evalCoord[dim] - floor(evalCoord[dim]);

        if (deriv[dim] == 0)
        {
            weights_1d[dim][0] =
                s *
                (0.046875 +
                 s * (-0.1875 + s * (0.28125 + (-0.1875 + 0.046875 * s) * s)));
            weights_1d[dim][1] =
                s *
                (-0.59375 +
                 s * (1.25 + s * (-0.5625 + (-0.296875 + 0.203125 * s) * s)));
            weights_1d[dim][2] =
                1. + (s * s) * (-2.125 + (1.96875 - 0.84375 * s) * (s * s));
            weights_1d[dim][3] =
                s * (0.59375 +
                     s * (1.25 + s * (0.5625 + (-2.25 + 0.84375 * s) * s)));
            weights_1d[dim][4] =
                s *
                (-0.046875 +
                 s * (-0.1875 + s * (-0.28125 + (0.71875 - 0.203125 * s) * s)));
            weights_1d[dim][5] =
                (0.046875 - 0.046875 * s) * ((s * s) * (s * s));
        }
        else if (deriv[dim] == 1)
        {
            weights_1d[dim][0] =
                (0.046875 +
                 s * (-0.375 + s * (0.84375 + (-0.75 + 0.234375 * s) * s))) /
                dx[dim];
            weights_1d[dim][1] =
                (-0.59375 +
                 s * (2.5 + s * (-1.6875 + s * (-1.1875 + 1.015625 * s)))) /
                dx[dim];
            weights_1d[dim][2] =
                (s * (-4.25 + (7.875 - 4.21875 * s) * (s * s))) / dx[dim];
            weights_1d[dim][3] =
                (0.59375 + s * (2.5 + s * (1.6875 + s * (-9. + 4.21875 * s)))) /
                dx[dim];
            weights_1d[dim][4] =
                (-0.046875 +
                 s * (-0.375 + s * (-0.84375 + (2.875 - 1.015625 * s) * s))) /
                dx[dim];
            weights_1d[dim][5] =
                ((0.1875 - 0.234375 * s) * s * (s * s)) / dx[dim];
        }
        else if (deriv[dim] == 2)
        {
            weights_1d[dim][0] =
                (-0.375 + s * (1.6875 + (-2.25 + 0.9375 * s) * s)) /
                (dx[dim] * dx[dim]);
            weights_1d[dim][1] =
                (2.5 + s * (-3.375 + s * (-3.5625 + 4.0625 * s))) /
                (dx[dim] * dx[dim]);
            weights_1d[dim][2] =
                (-4.25 + (23.625 - 16.875 * s) * (s * s)) / (dx[dim] * dx[dim]);
            weights_1d[dim][3] = (2.5 + s * (3.375 + s * (-27. + 16.875 * s))) /
                                 (dx[dim] * dx[dim]);
            weights_1d[dim][4] =
                (-0.375 + s * (-1.6875 + (8.625 - 4.0625 * s) * s)) /
                (dx[dim] * dx[dim]);
            weights_1d[dim][5] =
                ((0.5625 - 0.9375 * s) * (s * s)) / (dx[dim] * dx[dim]);
        }
        else
        {
            MayDay::Error("Quintic convolution algorithm only supports up to "
                          "second derivative");
        }
    }

    std::array<double, CH_SPACEDIM> interp_coord;

#if CH_SPACEDIM >= 3
    for (int z = 0; z < 6; ++z)
    {
        interp_coord[2] = floor(evalCoord[2]) + z - 2;
#endif
        for (int y = 0; y < 6; ++y)
        {
            interp_coord[1] = floor(evalCoord[1]) + y - 2;

            for (int x = 0; x < 6; ++x)
            {
                interp_coord[0] = floor(evalCoord[0]) + x - 2;
                CH_assert(m_source.contains(interp_coord));

                m_interp_points.push_back(IntVect(D_DECL6(
                    interp_coord[0], interp_coord[1], interp_coord[2],
                    interp_coord[3], interp_coord[4], interp_coord[5])));
                m_interp_weights.push_back(D_TERM6(weights_1d[0][x],
                                                   *weights_1d[1][y],
                                                   *weights_1d[2][z], , , ));
            }
        }
#if CH_SPACEDIM >= 3
    }
#endif
}

double QuinticConvolution::interpData(const FArrayBox &fab, int comp)
{
    long double accum = 0.0;

    for (int i = 0; i < m_interp_points.size(); ++i)
    {
        double data = m_interp_weights[i] * fab.get(m_interp_points[i], comp);
        accum += data;
    }

    return accum;
}

#endif /* QUINTICCONVOLUTION_IMPL_HPP_ */

/*
HIGHER ORDER WEIGHTS

s*(-0.009114583333333334 + s*(0.036458333333333336 + s*(-0.0546875 +
(0.036458333333333336 - 0.009114583333333334*s)*s))) s*(0.11979166666666667 +
s*(-0.2604166666666667 + s*(0.13541666666666666 + (0.040364583333333336 -
0.03515625*s)*s))) s*(-0.7122395833333334 + s*(1.2135416666666667 +
s*(-0.10677083333333333 + (-0.6979166666666666 + 0.3033854166666667*s)*s))) 1. +
(s*s)*(-1.9791666666666667 + (1.6497395833333333 - 0.6705729166666666*s)*(s*s))
s*(0.7122395833333334 + s*(1.2135416666666667 + s*(0.10677083333333333 +
(-1.703125 + 0.6705729166666666*s)*s))) s*(-0.11979166666666667 +
s*(-0.2604166666666667 + s*(-0.13541666666666666 + (0.8190104166666666 -
0.3033854166666667*s)*s))) s*(0.009114583333333334 + s*(0.036458333333333336 +
s*(0.0546875 + (-0.13541666666666666 + 0.03515625*s)*s)))
(-0.009114583333333334 + 0.009114583333333334*s)*((s*s)*(s*s))

-0.009114583333333334 + s*(0.07291666666666667 + s*(-0.1640625 +
(0.14583333333333334 - 0.045572916666666664*s)*s)) 0.11979166666666667 +
s*(-0.5208333333333334 + s*(0.40625 + (0.16145833333333334 - 0.17578125*s)*s))
-0.7122395833333334 + s*(2.4270833333333335 + s*(-0.3203125 +
s*(-2.7916666666666665 + 1.5169270833333333*s))) s*(-3.9583333333333335 +
(6.598958333333333 - 3.3528645833333335*s)*(s*s)) -0.7122395833333334 +
s*(-2.4270833333333335 + s*(-0.3203125 + (6.8125 - 3.3528645833333335*s)*s))
0.11979166666666667 + s*(0.5208333333333334 + s*(0.40625 +
s*(-3.2760416666666665 + 1.5169270833333333*s))) -0.009114583333333334 +
s*(-0.07291666666666667 + s*(-0.1640625 + (0.5416666666666666 -
0.17578125*s)*s)) (0.036458333333333336 - 0.045572916666666664*s)*s*(s*s)

0.07291666666666667 + s*(-0.328125 + (0.4375 - 0.18229166666666666*s)*s)
-0.5208333333333334 + s*(0.8125 + (0.484375 - 0.703125*s)*s)
2.4270833333333335 + s*(-0.640625 + s*(-8.375 + 6.067708333333333*s))
-3.9583333333333335 + (19.796875 - 13.411458333333334*s)*(s*s)
2.4270833333333335 + s*(0.640625 + s*(-20.4375 + 13.411458333333334*s))
-0.5208333333333334 + s*(-0.8125 + (9.828125 - 6.067708333333333*s)*s)
0.07291666666666667 + s*(0.328125 + (-1.625 + 0.703125*s)*s)
(-0.109375 + 0.18229166666666666*s)*(s*s)
*/
