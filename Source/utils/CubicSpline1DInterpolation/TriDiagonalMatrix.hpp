/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TRIDIAGONALMATRIX_HPP_
#define TRIDIAGONALMATRIX_HPP_

#include <vector>

class TriDiagonalMatrix
{
  public:
    TriDiagonalMatrix() {} // default constructor
    TriDiagonalMatrix(unsigned a_size) { resize(a_size); }

    inline unsigned size() const { return m_diag.size(); }

    void resize(unsigned a_size);

    double &upper(unsigned i);
    double upper(unsigned i) const;

    double &diag(unsigned i);
    double diag(unsigned i) const;

    double &lower(unsigned i);
    double lower(unsigned i) const;

    std::vector<double> lu_solve(const std::vector<double> &b);

  private:
    void lu_decomposition();
    std::vector<double> back_sub(const std::vector<double> &);
    std::vector<double> forward_sub(const std::vector<double> &);

  private:
    std::vector<double> m_upper; // upper band
    std::vector<double> m_diag;  // diagonal
    std::vector<double> m_lower; // lower band
};

#endif /* TRIDIAGONALMATRIX_HPP_ */
