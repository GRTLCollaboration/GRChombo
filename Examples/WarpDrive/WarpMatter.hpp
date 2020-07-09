/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WARPMATTER_HPP_
#define WARPMATTER_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "WarpBubble.hpp"
#include "WarpField.hpp"
#include "simd.hpp"

//! Class which creates a bubble of a Warp field given params for initial
//! matter config
class WarpMatter
{
  public:
    //! The constructor
    WarpMatter(WarpBubble::params_t a_params, double a_dx)
        : m_dx(a_dx), m_params(a_params)
    {
    }

  protected:
    const double m_dx;
    const WarpBubble::params_t
        m_params; //!< The matter initial condition params

  public:
    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // OLD CODE: values set during setup from constraints
        //        auto Ham = current_cell.load_vars(c_Ham);
        //        auto Mom1 = current_cell.load_vars(c_Mom1);
        //        auto Mom2 = current_cell.load_vars(c_Mom2);
        //        auto Mom3 = current_cell.load_vars(c_Mom3);
        Coordinates<data_t> coords(current_cell, m_dx, m_params.bubble_center);

        double t = 0.0;
        double vs = m_params.warp_speed;
        double xs = vs * t;
        double sigma = m_params.sigma_wall;
        double sigma2 = sigma * sigma;
        double R0 = m_params.bubble_size;

        // coords
        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;
        data_t rho2 = y * y + z * z;
        data_t r_s = sqrt((x - xs) * (x - xs) + rho2);
        data_t r_s2 = r_s * r_s;
        data_t r_s3 = r_s2 * r_s;

        // derivatives
        data_t drsdx = (x - xs) / r_s;
        data_t drsdy = y / r_s;
        data_t drsdz = z / r_s;
        data_t drsdt = -vs * (x - xs) / r_s;
        data_t d2rsdx2 = 1.0 / r_s - (x - xs) * (x - xs) / r_s3;
        data_t d2rsdy2 = 1.0 / r_s - y * y / r_s3;
        data_t d2rsdz2 = 1.0 / r_s - z * z / r_s3;
        data_t d2rsdt2 = vs * vs / r_s - vs * vs * (x - xs) * (x - xs) / r_s3;
        data_t d2rsdxdy = -(x - xs) * y / r_s3;
        data_t d2rsdxdz = -(x - xs) * z / r_s3;
        data_t d2rsdydz = -y * z / r_s3;
        data_t d2rsdtdx = -vs / r_s + vs * (x - xs) * (x - xs) / r_s3;
        data_t d2rsdtdy = vs * (x - xs) * y / r_s3;
        data_t d2rsdtdz = vs * (x - xs) * z / r_s3;

        // useful functions
        data_t Xp = (r_s + R0) * sigma;
        data_t Xm = (r_s - R0) * sigma;
        data_t fp = tanh(Xp);
        data_t dfp = 1 - fp * fp;
        data_t fm = tanh(Xm);
        data_t dfm = 1 - fm * fm;
        data_t A = -0.5 * vs / tanh(R0 * sigma);
        data_t betax = A * (fp - fm);
        data_t dbetadx = A * (dfp - dfm) * drsdx * sigma;
        data_t dbetady = A * (dfp - dfm) * drsdy * sigma;
        data_t dbetadz = A * (dfp - dfm) * drsdz * sigma;
        data_t dbetadt = A * (dfp - dfm) * drsdt * sigma;
        data_t d2betadx2 =
            2.0 * A * (fm * dfm - fp * dfp) * drsdx * drsdx * sigma2 +
            A * (dfp - dfm) * d2rsdx2 * sigma;
        data_t d2betady2 =
            2.0 * A * (fm * dfm - fp * dfp) * drsdy * drsdy * sigma2 +
            A * (dfp - dfm) * d2rsdy2 * sigma;
        data_t d2betadz2 =
            2.0 * A * (fm * dfm - fp * dfp) * drsdz * drsdz * sigma2 +
            A * (dfp - dfm) * d2rsdz2 * sigma;
        data_t d2betadt2 =
            2.0 * A * (fm * dfm - fp * dfp) * drsdt * drsdt * sigma2 +
            A * (dfp - dfm) * d2rsdt2 * sigma;
        data_t d2betadxdy =
            2.0 * A * (fm * dfm - fp * dfp) * drsdx * drsdy * sigma2 +
            A * (dfp - dfm) * d2rsdxdy * sigma;
        data_t d2betadxdz =
            2.0 * A * (fm * dfm - fp * dfp) * drsdx * drsdz * sigma2 +
            A * (dfp - dfm) * d2rsdxdz * sigma;
        data_t d2betadtdx =
            2.0 * A * (fm * dfm - fp * dfp) * drsdt * drsdx * sigma2 +
            A * (dfp - dfm) * d2rsdtdx * sigma;
        data_t d2betadtdy =
            2.0 * A * (fm * dfm - fp * dfp) * drsdt * drsdy * sigma2 +
            A * (dfp - dfm) * d2rsdtdy * sigma;
        data_t d2betadtdz =
            2.0 * A * (fm * dfm - fp * dfp) * drsdt * drsdz * sigma2 +
            A * (dfp - dfm) * d2rsdtdz * sigma;

        // analytic solutions for the T_ab comps
        // Note that the index 3 it the time index, 0-2 is spatial
        std::array<std::array<data_t, 4>, 4> G;
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                G[i][j] = 0;
            }
        }
        // time comps first
        data_t dbetadrho2 = dbetady * dbetady + dbetadz * dbetadz;
        data_t d2betadrho2 = d2betady2 + d2betadz2;
        G[3][0] = 0.25 * (-2.0 * d2betadrho2 - 3.0 * betax * dbetadrho2);
        G[3][1] = 0.5 * (-(betax * d2betadtdy) + d2betadxdy +
                         betax * betax * d2betadxdy +
                         2.0 * betax * dbetadx * dbetady);
        G[3][2] = 0.5 * (-(betax * d2betadtdz) + d2betadxdz +
                         betax * betax * d2betadxdz +
                         2.0 * betax * dbetadx * dbetadz);
        G[3][3] = 0.25 * (d2betadrho2 * (-4.0 * betax) -
                          dbetadrho2 * (1.0 + 3.0 * betax * betax));
        // Now the x ones
        G[0][0] = -0.75 * dbetadrho2;
        G[0][1] =
            0.5 * (-d2betadtdy + betax * d2betadxdy + 2.0 * dbetadx * dbetady);
        G[0][2] =
            0.5 * (-d2betadtdz + betax * d2betadxdz + 2.0 * dbetadx * dbetadz);
        G[0][3] = G[3][0];

        // Now the y ones
        G[1][0] = G[0][1];
        G[1][1] = 0.25 * (4.0 * d2betadtdx - 4.0 * betax * d2betadx2 -
                          4.0 * dbetadx * dbetadx + dbetady * dbetady -
                          dbetadz * dbetadz);
        G[1][2] = 0.5 * dbetady * dbetadz;
        G[1][3] = G[3][1];

        // Now the z ones
        G[2][0] = G[0][2];
        G[2][1] = G[1][2];
        G[2][2] = 0.25 * (4.0 * d2betadtdx - 4.0 * betax * d2betadx2 -
                          4.0 * dbetadx * dbetadx - dbetady * dbetady +
                          dbetadz * dbetadz);
        G[2][3] = G[3][2];

        // Now assign the stress energy tensor
        data_t rhoM =
            (G[3][3] - 2.0 * betax * G[3][0] + betax * betax * G[0][0]) / 8.0 /
            M_PI;
        Tensor<1, data_t> Si;
        FOR1(i) { Si[i] = (-G[3][i] + betax * G[0][i]) / 8.0 / M_PI; }
        /*
                // DEBUG - NB need to disable simd in WarpFieldLevel call
                if((coords.y == 1.5) && (coords.z == 2.5)
                    && (coords.x == -3.5))
                {
                    pout() << "x, y, z " << x << " " << y << " " << z << endl;
                    FOR2(i,j)
                    {
                        pout() << "G0i " << i << " " << G[3][i] << endl;
                        pout() << "G0j " << j << " " << G[3][j] << endl;
                        pout() << "Gij " << G[i][j] << endl;
                    }
                    pout() << "betax " << betax << endl;
                    pout() << "dbetadx " << dbetadx << endl;
                    pout() << "dbetady " << dbetady << endl;
                    pout() << "dbetadz " << dbetadz << endl;
                    pout() << "dbetadt " << dbetadt << endl;
                    pout() << "d2betadx2 " << d2betadx2 << endl;
                    pout() << "d2betady2 " << d2betady2 << endl;
                    pout() << "d2betadz2 " << d2betadz2 << endl;
                    pout() << "d2betadt2 " << d2betadt2 << endl;
                    pout() << "d2betadxdz " << d2betadxdz << endl;
                    pout() << "d2betadxdy " << d2betadxdy << endl;
                    pout() << "d2betadtdx " << d2betadtdx << endl;
                    pout() << "d2betadtdy " << d2betadtdy << endl;
                    pout() << "d2betadtdz " << d2betadtdz << endl;
                }
        */
        // Now store them
        current_cell.store_vars(rhoM, c_rho);
        current_cell.store_vars(Si[0], c_S1);
        current_cell.store_vars(Si[1], c_S2);
        current_cell.store_vars(Si[2], c_S3);
        current_cell.store_vars(G[0][0] / 8.0 / M_PI, c_S11);
        current_cell.store_vars(G[0][1] / 8.0 / M_PI, c_S12);
        current_cell.store_vars(G[0][2] / 8.0 / M_PI, c_S13);
        current_cell.store_vars(G[1][1] / 8.0 / M_PI, c_S22);
        current_cell.store_vars(G[1][2] / 8.0 / M_PI, c_S23);
        current_cell.store_vars(G[2][2] / 8.0 / M_PI, c_S33);
    }
};

#endif /* WARPMATTER_HPP_ */
