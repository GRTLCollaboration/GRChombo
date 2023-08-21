/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TENSORALGEBRA_HPP_
#define TENSORALGEBRA_HPP_

#include "AlwaysInline.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include <array>

template <class data_t> struct chris_t
{
    Tensor<3, data_t> ULL;        //!< standard christoffel symbols
    Tensor<3, data_t> LLL;        //!< 3 lower indices
    Tensor<1, data_t> contracted; //!< contracted christoffel
};

namespace TensorAlgebra
{
/// Computes the determinant of a general 1x1 matrix
template <class data_t>
ALWAYS_INLINE data_t compute_determinant(const Tensor<2, data_t, 1> &matrix)
{
    return matrix[0][0];
}

/// Computes the determinant of a general 2x2 matrix
template <class data_t>
ALWAYS_INLINE data_t compute_determinant(const Tensor<2, data_t, 2> &matrix)
{
    data_t det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    return det;
}

/// Computes the determinant of a general 3x3 matrix.
/// Note: for a symmetric matrix use the simplified function
template <class data_t>
ALWAYS_INLINE data_t compute_determinant(const Tensor<2, data_t, 3> &matrix)
{
    data_t det =
        matrix[0][0] *
            (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
        matrix[0][1] *
            (matrix[2][2] * matrix[1][0] - matrix[1][2] * matrix[2][0]) +
        matrix[0][2] *
            (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
    return det;
}

/// Computes determinant of a symmetric 3x3 matrix
template <class data_t>
ALWAYS_INLINE data_t compute_determinant_sym(const Tensor<2, data_t, 3> &matrix)
{
    data_t det = matrix[0][0] * matrix[1][1] * matrix[2][2] +
                 2 * matrix[0][1] * matrix[0][2] * matrix[1][2] -
                 matrix[0][0] * matrix[1][2] * matrix[1][2] -
                 matrix[1][1] * matrix[0][2] * matrix[0][2] -
                 matrix[2][2] * matrix[0][1] * matrix[0][1];

    return det;
}

template <class data_t>
ALWAYS_INLINE data_t
compute_determinant_sym(const Tensor<2, data_t, 2>
                            &matrix) // This function only works for 2D matrix
{
    data_t det = matrix[1][1] * matrix[0][0] - matrix[0][1] * matrix[0][1];
    return det;
}

/// Computes the inverse of a symmetric 3x3 matrix directly.
template <class data_t>
Tensor<2, data_t, 3> compute_inverse_sym(const Tensor<2, data_t, 3> &matrix)
{
    data_t deth = compute_determinant_sym(matrix);
    data_t deth_inverse = 1. / deth;
    Tensor<2, data_t, 3> h_UU;
    h_UU[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[1][2]) *
                 deth_inverse;
    h_UU[0][1] = (matrix[0][2] * matrix[1][2] - matrix[0][1] * matrix[2][2]) *
                 deth_inverse;
    h_UU[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) *
                 deth_inverse;
    h_UU[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[0][2]) *
                 deth_inverse;
    h_UU[1][2] = (matrix[0][1] * matrix[0][2] - matrix[0][0] * matrix[1][2]) *
                 deth_inverse;
    h_UU[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[0][1]) *
                 deth_inverse;
    h_UU[1][0] = h_UU[0][1];
    h_UU[2][0] = h_UU[0][2];
    h_UU[2][1] = h_UU[1][2];

    return h_UU;
}

template <class data_t>
Tensor<2, data_t, 2>
compute_inverse_sym(const Tensor<2, data_t, 2>
                        &matrix) // This function only works for 2D matrix
{
    data_t deth = compute_determinant_sym(matrix);
    data_t deth_inverse = 1. / deth;
    Tensor<2, data_t, 2> h_UU;
    h_UU[0][0] = matrix[1][1] * deth_inverse;
    h_UU[0][1] = -matrix[0][1] * deth_inverse;
    h_UU[1][1] = matrix[0][0] * deth_inverse;
    h_UU[1][0] = h_UU[0][1];

    return h_UU;
}

/// Computes the inverse of a general 3x3 matrix.
/// Note: for a symmetric matrix use the simplified function
template <class data_t>
Tensor<2, data_t, 3> compute_inverse(const Tensor<2, data_t, 3> &matrix)
{
    data_t deth = compute_determinant(matrix);
    data_t deth_inverse = 1. / deth;
    Tensor<2, data_t, 3> h_UU;
    h_UU[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) *
                 deth_inverse;
    h_UU[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) *
                 deth_inverse;
    h_UU[2][2] = (matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]) *
                 deth_inverse;
    h_UU[1][0] = (matrix[2][0] * matrix[1][2] - matrix[1][0] * matrix[2][2]) *
                 deth_inverse;
    h_UU[0][1] = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) *
                 deth_inverse;
    h_UU[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) *
                 deth_inverse;
    h_UU[0][2] = (matrix[0][1] * matrix[1][2] - matrix[1][1] * matrix[0][2]) *
                 deth_inverse;
    h_UU[2][1] = (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) *
                 deth_inverse;
    h_UU[1][2] = (matrix[1][0] * matrix[0][2] - matrix[0][0] * matrix[1][2]) *
                 deth_inverse;

    return h_UU;
}

template <class data_t>
Tensor<2, data_t, 2>
compute_inverse(const Tensor<2, data_t, 2>
                    &matrix) // This function only works for 2D matrix
{
    data_t deth = compute_determinant(matrix);
    data_t deth_inverse = 1. / deth;
    Tensor<2, data_t, 2> h_UU;
    h_UU[0][0] = matrix[1][1] * deth_inverse;
    h_UU[1][1] = matrix[0][0] * deth_inverse;
    h_UU[0][1] = -matrix[0][1] * deth_inverse;
    h_UU[1][0] = -matrix[1][0] * deth_inverse;

    return h_UU;
}

/// Computes the trace of a 2-Tensor with lower indices given an inverse metric.
template <class data_t>
ALWAYS_INLINE data_t compute_trace(const Tensor<2, data_t> &tensor_LL,
                                   const Tensor<2, data_t> &inverse_metric)
{
    data_t trace = 0.;
    FOR(i, j) { trace += inverse_metric[i][j] * tensor_LL[i][j]; }
    return trace;
}

/// Computes the trace of a 1,1 Tensor (a matrix) - no metric required.
template <class data_t>
ALWAYS_INLINE data_t compute_trace(const Tensor<2, data_t> &tensor_UL)
{
    data_t trace = 0.;
    FOR(i) trace += tensor_UL[i][i];
    return trace;
}

template <class data_t>
ALWAYS_INLINE data_t
compute_trace(const Tensor<1, Tensor<1, data_t>> &tensor_UL)
{
    data_t trace = 0.;
    FOR(i) trace += tensor_UL[i][i];
    return trace;
}

/// Computes dot product of a vector and a covector (no metric required)
template <class data_t>
ALWAYS_INLINE data_t compute_dot_product(const Tensor<1, data_t> &vector_U,
                                         const Tensor<1, data_t> &covector_L)
{
    data_t dot_product = 0.;
    FOR(i) dot_product += vector_U[i] * covector_L[i];
    return dot_product;
}

/// Computes dot product of two covectors given an inverse metric or
/// the dot product of two vectors given a metric.
template <class data_t>
ALWAYS_INLINE data_t compute_dot_product(
    const Tensor<1, data_t> &covector1_L, const Tensor<1, data_t> &covector2_L,
    const Tensor<2, data_t> &inverse_metric)
{
    data_t dot_product = 0.;
    FOR(m, n)
    {
        dot_product += inverse_metric[m][n] * covector1_L[m] * covector2_L[n];
    }
    return dot_product;
}

/// Removes the trace of a 2-Tensor with lower indices given a metric and an
/// inverse metric.  Or a Tensor with upper indices given an inverse metric and
/// a metric.
template <class data_t>
ALWAYS_INLINE void make_trace_free(Tensor<2, data_t> &tensor_LL,
                                   const Tensor<2, data_t> &metric,
                                   const Tensor<2, data_t> &inverse_metric)
{
    auto trace = compute_trace(tensor_LL, inverse_metric);
    double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    FOR(i, j)
    {
        tensor_LL[i][j] -= one_over_gr_spacedim * metric[i][j] * trace;
    }
}

/// Makes a 2-Tensor symmetric
template <class data_t, int size>
ALWAYS_INLINE void make_symmetric(Tensor<2, data_t, size> &tensor_LL)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            tensor_LL[i][j] = 0.5 * (tensor_LL[i][j] + tensor_LL[j][i]);
            tensor_LL[j][i] = tensor_LL[i][j];
        }
    }
}

/// Raises the index of a covector
template <class data_t>
ALWAYS_INLINE Tensor<1, data_t>
raise_all(const Tensor<1, data_t> &tensor_L,
          const Tensor<2, data_t> &inverse_metric)
{
    Tensor<1, data_t> tensor_U = 0.;
    FOR(i, j) { tensor_U[i] += inverse_metric[i][j] * tensor_L[j]; }
    return tensor_U;
}

/// Raises the indices of a 2-Tensor
template <class data_t>
ALWAYS_INLINE Tensor<2, data_t>
raise_all(const Tensor<2, data_t> &tensor_LL,
          const Tensor<2, data_t> &inverse_metric)
{
    Tensor<2, data_t> tensor_UU = 0.;
    FOR(i, j, k, l)
    {
        tensor_UU[i][j] +=
            inverse_metric[i][k] * inverse_metric[j][l] * tensor_LL[k][l];
    }
    return tensor_UU;
}

/// Lowers the indices of a vector
/// Note: same functionality as raise; included to improve readibility
template <class data_t>
ALWAYS_INLINE Tensor<1, data_t> lower_all(const Tensor<1, data_t> &tensor_U,
                                          const Tensor<2, data_t> &metric)
{ // The code for lowering is exactly the same as for raising
    return raise_all(tensor_U, metric);
}

/// Lowers the indices of a 2-Tensor
/// Note: same functionality as raise; included to improve readibility
template <class data_t>
ALWAYS_INLINE Tensor<2, data_t> lower_all(const Tensor<2, data_t> &tensor_UU,
                                          const Tensor<2, data_t> &metric)
{ // The code for lowering is exactly the same as for raising
    return raise_all(tensor_UU, metric);
}

/// Computes the (i,j) component of the Kronecker delta
constexpr int delta(int i, int j) { return (i == j); }

/// Computes the levi-civita symbol (3D, NB, symbol, not the Tensor)
inline Tensor<3, double, 3> epsilon()
{
    Tensor<3, double, 3> epsilon = {0.};
    epsilon[0][1][2] = 1.0;
    epsilon[1][2][0] = 1.0;
    epsilon[2][0][1] = 1.0;
    epsilon[0][2][1] = -1.0;
    epsilon[2][1][0] = -1.0;
    epsilon[1][0][2] = -1.0;

    return epsilon;
}

/// Computes the levi-civita symbol (4D, NB, symbol, not the Tensor)
inline Tensor<4, double, 4> epsilon4D()
{
    Tensor<4, double, 4> epsilon4D = {0.0};
    epsilon4D[0][1][2][3] = 1.0;
    epsilon4D[0][1][3][2] = -1.0;
    epsilon4D[0][3][1][2] = 1.0;
    epsilon4D[0][3][2][1] = -1.0;
    epsilon4D[0][2][1][3] = -1.0;
    epsilon4D[0][2][3][1] = 1.0;

    epsilon4D[1][0][2][3] = -1.0;
    epsilon4D[1][2][0][3] = 1.0;
    epsilon4D[1][2][3][0] = -1.0;
    epsilon4D[1][3][2][0] = 1.0;
    epsilon4D[1][3][0][2] = -1.0;
    epsilon4D[1][0][3][2] = 1.0;

    epsilon4D[2][0][1][3] = 1.0;
    epsilon4D[2][0][3][1] = -1.0;
    epsilon4D[2][3][0][1] = 1.0;
    epsilon4D[2][3][1][0] = -1.0;
    epsilon4D[2][1][3][0] = 1.0;
    epsilon4D[2][1][0][3] = -1.0;

    epsilon4D[3][0][1][2] = -1.0;
    epsilon4D[3][1][0][2] = 1.0;
    epsilon4D[3][1][2][0] = -1.0;
    epsilon4D[3][2][1][0] = 1.0;
    epsilon4D[3][2][0][1] = -1.0;
    epsilon4D[3][0][2][1] = 1.0;

    return epsilon4D;
}

/// Computes the conformal christoffel symbol
template <class data_t>
chris_t<data_t>
compute_christoffel(const Tensor<2, Tensor<1, data_t>> &d1_metric,
                    const Tensor<2, data_t> &h_UU)
{
    chris_t<data_t> out;

    FOR(i, j, k)
    {
        out.LLL[i][j][k] = 0.5 * (d1_metric[j][i][k] + d1_metric[k][i][j] -
                                  d1_metric[j][k][i]);
    }
    FOR(i, j, k)
    {
        out.ULL[i][j][k] = 0;
        FOR(l) { out.ULL[i][j][k] += h_UU[i][l] * out.LLL[l][j][k]; }
    }
    FOR(i)
    {
        out.contracted[i] = 0;
        FOR(j, k) { out.contracted[i] += h_UU[j][k] * out.ULL[i][j][k]; }
    }

    return out;
}

template <class data_t>
Tensor<3, data_t> compute_phys_chris(const Tensor<1, data_t> &d1_chi,
                                     const data_t &vars_chi,
                                     const Tensor<2, data_t> &vars_h,
                                     const Tensor<2, data_t> &h_UU,
                                     const Tensor<3, data_t> &chris_ULL)
{
    Tensor<3, data_t> chris_phys;
    FOR(i, j, k)
    {
        chris_phys[i][j][k] =
            chris_ULL[i][j][k] -
            0.5 / vars_chi *
                (delta(i, k) * d1_chi[j] + delta(i, j) * d1_chi[k]);
        FOR(m)
        {
            chris_phys[i][j][k] +=
                0.5 / vars_chi * vars_h[j][k] * h_UU[i][m] * d1_chi[m];
        }
    }
    return chris_phys;
}
} // namespace TensorAlgebra

#endif /* TENSORALGEBRA_HPP_ */
