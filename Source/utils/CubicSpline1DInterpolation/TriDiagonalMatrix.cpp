/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "TriDiagonalMatrix.hpp"
#include <cassert>

void TriDiagonalMatrix::resize(unsigned a_size)
{
    assert(a_size > 0);
    m_upper.resize(a_size - 1);
    m_diag.resize(a_size);
    m_lower.resize(a_size - 1);
}

double &TriDiagonalMatrix::upper(unsigned i)
{
    assert(i < size() - 1);
    return m_upper[i];
}
double TriDiagonalMatrix::upper(unsigned i) const
{
    assert(i < size() - 1);
    return m_upper[i];
}

double &TriDiagonalMatrix::diag(unsigned i)
{
    assert(i < size());
    return m_diag[i];
}
double TriDiagonalMatrix::diag(unsigned i) const
{
    assert(i < size());
    return m_diag[i];
}

double &TriDiagonalMatrix::lower(unsigned i)
{
    assert(i < size() - 1);
    return m_lower[i];
}
double TriDiagonalMatrix::lower(unsigned i) const
{
    assert(i < size() - 1);
    return m_lower[i];
}

std::vector<double> TriDiagonalMatrix::lu_solve(const std::vector<double> &b)
{
    lu_decomposition();
    return back_sub(forward_sub(b));
}

void TriDiagonalMatrix::lu_decomposition()
{
    for (unsigned i = 1; i < size(); i++)
    {
        lower(i - 1) /= diag(i - 1);
        diag(i) -= lower(i - 1) * upper(i - 1);
    }
}

std::vector<double> TriDiagonalMatrix::back_sub(const std::vector<double> &b)
{
    assert(size() == b.size());
    std::vector<double> out(b);
    out[size() - 1] /= diag(size() - 1);
    for (int i = size() - 2; i >= 0; i--)
        out[i] = (out[i] - upper(i) * out[i + 1]) / diag(i);
    return out;
}
std::vector<double> TriDiagonalMatrix::forward_sub(const std::vector<double> &d)
{
    assert(size() == d.size());
    std::vector<double> out(d);
    for (unsigned i = 1; i < size(); i++)
        out[i] -= lower(i - 1) * out[i - 1];
    return out;
}
