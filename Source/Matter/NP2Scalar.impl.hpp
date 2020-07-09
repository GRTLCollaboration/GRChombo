/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(NP2SCALAR_HPP_)
#error "This file should only be included through NP2Scalar.hpp"
#endif

#ifndef NP2SCALAR_IMPL_HPP_
#define NP2SCALAR_IMPL_HPP_

template <class data_t> void NP2Scalar::compute(Cell<data_t> current_cell) const
{
    // copy data from chombo gridpoint into local variables
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Diff1Vars>(current_cell);

    // Get the coordinates, about centre of grid
    std::array<double, 3> center;
    center.fill(0.5 * m_L);
    Coordinates<data_t> coords(current_cell, m_dx, center);

    // work out the Newman Penrose scalar
    NPScalar_t<data_t> out = compute_NP_scalar_2(vars, d1, coords);

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(out.Real, c_ReNP2);
    current_cell.store_vars(out.Im, c_ImNP2);
}

// Calculation of the Newman Penrose scalar for a vector field Avec
// Phi_2 = F_{\mu \nu} * l^\mu \bar(m)^\nu, see e.g. arxiv 1212.0551 eqn 32
// NB for now this assumes the field is purely real
template <class data_t>
NPScalar_t<data_t>
NP2Scalar::compute_NP_scalar_2(const Vars<data_t> &vars,
                               const Diff1Vars<Tensor<1, data_t>> &d1,
                               const Coordinates<data_t> &coords) const
{
    NPScalar_t<data_t> out;

    // Real and Imaginary parts of the scalar
    out.Real = 0.0;
    out.Im = 0.0;

    // set up the tetrad
    Tetrad_t<data_t> tetrad = compute_null_tetrad(vars, coords);

    // calculate the projections
    FOR1(i)
    {
        // Assumes no Imaginary components to the fields
        out.Real += -0.5 * (vars.Evec[i] * tetrad.v[i]);
        out.Im += 0.5 * (vars.Evec[i] * tetrad.w[i]);

        FOR1(j)
        {
            out.Real += -0.5 * tetrad.u[i] * tetrad.v[j] *
                        (d1.Avec[j][i] - d1.Avec[i][j]);
            out.Im += 0.5 * tetrad.u[i] * tetrad.w[j] *
                      (d1.Avec[j][i] - d1.Avec[i][j]);
        }
    }

    out.magnitude = sqrt(out.Real * out.Real + out.Im * out.Im);

    return out;
}

// Calculation of the null tetrad
template <class data_t>
Tetrad_t<data_t>
NP2Scalar::compute_null_tetrad(const Vars<data_t> &vars,
                               const Coordinates<data_t> &coords) const
{
    Tetrad_t<data_t> out;

    // compute position on grid relative to center
    const data_t x = coords.x;
    const double y = coords.y;
    const double z = coords.z;
    const data_t r = coords.get_radius();

    // metric in raised form
    const auto h_UU = TensorAlgebra::compute_inverse(vars.h);

    // the levi civita symbol
    const Tensor<3, double> epsilon = TensorAlgebra::epsilon();

    // calculate tetrad
    // u is the radial vector
    out.u[0] = x;
    out.u[1] = y;
    out.u[2] = z;

    // v is the azimuthal vector
    out.v[0] = -y;
    out.v[1] = x;
    out.v[2] = 0.0;

    // floor on chi
    const data_t chi = simd_max(vars.chi, 1e-6);

    // w is the polar vector
    FOR1(i)
    {
        out.w[i] = 0.0;
        FOR3(j, k, m)
        {
            out.w[i] += pow(chi, -0.5) * h_UU[i][j] * epsilon[j][k][m] *
                        out.v[k] * out.u[m];
        }
    }

    // Gram Schmitt orthonormalisation
    // Usual to start with azimuthal direction (because of frame dragging)
    // then radial then polar
    data_t omega_11 = 0.0;
    FOR2(i, j) { omega_11 += out.v[i] * out.v[j] * vars.h[i][j] / chi; }
    FOR1(i) { out.v[i] = out.v[i] / sqrt(omega_11); }

    data_t omega_12 = 0.0;
    FOR2(i, j) { omega_12 += out.v[i] * out.u[j] * vars.h[i][j] / chi; }
    FOR1(i) { out.u[i] += -omega_12 * out.v[i]; }

    data_t omega_22 = 0.0;
    FOR2(i, j) { omega_22 += out.u[i] * out.u[j] * vars.h[i][j] / chi; }
    FOR1(i) { out.u[i] = out.u[i] / sqrt(omega_22); }

    data_t omega_13 = 0.0;
    data_t omega_23 = 0.0;
    FOR2(i, j)
    {
        omega_13 += out.v[i] * out.w[j] * vars.h[i][j] / chi;
        omega_23 += out.u[i] * out.w[j] * vars.h[i][j] / chi;
    }
    FOR1(i) { out.w[i] += -(omega_13 * out.v[i] + omega_23 * out.u[i]); }

    data_t omega_33 = 0.0;
    FOR2(i, j) { omega_33 += out.w[i] * out.w[j] * vars.h[i][j] / chi; }

    FOR1(i) { out.w[i] = out.w[i] / sqrt(omega_33); }

    return out;
}

#endif /* NP2SCALAR_HPP_ */
