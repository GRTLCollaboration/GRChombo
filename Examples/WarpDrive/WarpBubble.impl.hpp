/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(WARPBUBBLE_HPP_)
#error "This file should only be included through WarpBubble.hpp"
#endif

#ifndef WARPBUBBLE_IMPL_HPP_
#define WARPBUBBLE_IMPL_HPP_

// Compute the value of the initial vars on the grid
template <class data_t>
void WarpBubble::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<WarpField>::Vars<data_t> vars;
    VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx, m_params.bubble_center);

    double t = 0.0;
    double vs = m_params.warp_speed;
    double xs = vs * t;
    double sigma = m_params.sigma_wall;
    double R0 = m_params.bubble_size;

    // coords
    data_t x = coords.x;
    double y = coords.y;
    double z = coords.z;
    double rho2 = y * y + z * z;
    data_t r_s = sqrt((x - xs) * (x - xs) + rho2);

    // useful functions
    data_t Xp = (r_s + R0) * sigma;
    data_t Xm = (r_s - R0) * sigma;
    data_t drsdx = (x - xs) / r_s;
    data_t drsdy = y / r_s;
    data_t drsdz = z / r_s;
    data_t fp = tanh(Xp);
    data_t dfp = 1 - fp * fp;
    data_t fm = tanh(Xm);
    data_t dfm = 1 - fm * fm;
    data_t A = -0.5 * vs / tanh(R0 * sigma);

    // set the shift and metric vars
    vars.shift[0] = A * (fp - fm);
    vars.lapse = 1;
    // spatial metric is flat, K_ij is not
    vars.chi = 1;
    FOR2(i, j) { vars.h[i][j] = TensorAlgebra::delta(i, j); }
    data_t dbetadx = A * (dfp - dfm) * drsdx * sigma;
    data_t dbetady = A * (dfp - dfm) * drsdy * sigma;
    data_t dbetadz = A * (dfp - dfm) * drsdz * sigma;
    vars.K = dbetadx;

    // set Aij
    FOR2(i, j) vars.A[i][j] = 0.0;
    vars.A[0][0] = vars.K;
    vars.A[0][1] = 0.5 * dbetady;
    vars.A[0][2] = 0.5 * dbetadz;
    FOR1(i) vars.A[i][i] += -1.0 / 3.0 * vars.K;
    vars.A[1][0] = vars.A[0][1];
    vars.A[2][0] = vars.A[0][2];

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* WARPBUBBLE_IMPL_HPP_ */
