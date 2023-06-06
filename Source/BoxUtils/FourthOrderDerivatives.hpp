/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FOURTHORDERDERIVATIVES_HPP_
#define FOURTHORDERDERIVATIVES_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp"
#include <array>

class FourthOrderDerivatives
{
  public:
    const double m_dx;

  private:
    const double m_one_over_dx;
    const double m_one_over_dx2;

  public:
    FourthOrderDerivatives(double dx)
        : m_dx(dx), m_one_over_dx(1 / dx), m_one_over_dx2(1 / (dx * dx))
    {
    }

    template <class data_t>
    ALWAYS_INLINE data_t diff1(const double *in_ptr, const int idx,
                               const int stride) const
    {
        auto in = SIMDIFY<data_t>(in_ptr);

        data_t weight_far = 8.33333333333333333333e-2;
        data_t weight_near = 6.66666666666666666667e-1;

        // NOTE: if you have been sent here by the debugger because of
        // EXC_BAD_ACCESS  or something similar you might be trying to take
        // derivatives without ghost points.
        return (weight_far * in[idx - 2 * stride] -
                weight_near * in[idx - stride] +
                weight_near * in[idx + stride] -
                weight_far * in[idx + 2 * stride]) *
               m_one_over_dx;
    }

    // Writes directly into the vars object - use this wherever possible
    template <class data_t, template <typename> class vars_t>
    void diff1(vars_t<Tensor<1, data_t>> &d1, const Cell<data_t> &current_cell,
               int direction) const
    {
        const int stride =
            current_cell.get_box_pointers().m_in_stride[direction];
        const int in_index = current_cell.get_in_index();
        d1.enum_mapping(
            [&](const int &ivar, Tensor<1, data_t> &var)
            {
                var[direction] = diff1<data_t>(
                    current_cell.get_box_pointers().m_in_ptr[ivar], in_index,
                    stride);
            });
    }

    /// Calculates all first derivatives and returns as variable type specified
    /// by the template parameter
    template <template <typename> class vars_t, class data_t>
    auto diff1(const Cell<data_t> &current_cell) const
    {
        const auto in_index = current_cell.get_in_index();
        const auto strides = current_cell.get_box_pointers().m_in_stride;
        vars_t<Tensor<1, data_t>> d1;
        d1.enum_mapping(
            [&](const int &ivar, Tensor<1, data_t> &var)
            {
                FOR(idir)
                {
                    var[idir] = diff1<data_t>(
                        current_cell.get_box_pointers().m_in_ptr[ivar],
                        in_index, strides[idir]);
                }
            });
        return d1;
    }

    template <class data_t>
    void diff1(Tensor<1, data_t> &diff_value, const Cell<data_t> &current_cell,
               int direction, int ivar) const
    {
        const int stride =
            current_cell.get_box_pointers().m_in_stride[direction];
        const int in_index = current_cell.get_in_index();
        diff_value[direction] = diff1<data_t>(
            current_cell.get_box_pointers().m_in_ptr[ivar], in_index, stride);
    }

    template <class data_t, int num_vars>
    void diff1(Tensor<1, data_t> (&diff_array)[num_vars],
               const Cell<data_t> &current_cell, int direction,
               int start_var = 0) const
    {
        const int stride =
            current_cell.get_box_pointers().m_in_stride[direction];
        const int in_index = current_cell.get_in_index();
        for (int i = start_var; i < start_var + num_vars; ++i)
        {
            diff_array[i][direction] = diff1<data_t>(
                current_cell.get_box_pointers().m_in_ptr[i], in_index, stride);
        }
    }

    template <class data_t>
    ALWAYS_INLINE data_t diff2(const double *in_ptr, const int idx,
                               const int stride) const
    {
        auto in = SIMDIFY<data_t>(in_ptr);

        data_t weight_far = 8.33333333333333333333e-2;
        data_t weight_near = 1.33333333333333333333e+0;
        data_t weight_local = 2.50000000000000000000e+0;

        return (-weight_far * in[idx - 2 * stride] +
                weight_near * in[idx - stride] - weight_local * in[idx] +
                weight_near * in[idx + stride] -
                weight_far * in[idx + 2 * stride]) *
               m_one_over_dx2;
    }

    // Writes 2nd deriv directly into the vars object - use this wherever
    // possible
    template <class data_t, template <typename> class vars_t>
    void diff2(vars_t<Tensor<2, data_t>> &d2, const Cell<data_t> &current_cell,
               int direction) const
    {
        const int stride =
            current_cell.get_box_pointers().m_in_stride[direction];
        const int in_index = current_cell.get_in_index();
        d2.enum_mapping(
            [&](const int &ivar, Tensor<2, data_t> &var)
            {
                var[direction][direction] = diff2<data_t>(
                    current_cell.get_box_pointers().m_in_ptr[ivar], in_index,
                    stride);
            });
    }

    template <class data_t, int num_vars>
    void diff2(Tensor<2, data_t> (&diffArray)[num_vars],
               const Cell<data_t> &current_cell, int direction) const
    {
        const int stride =
            current_cell.get_box_pointers().m_in_stride[direction];
        const int in_index = current_cell.get_in_index();
        for (int ivar = 0; ivar < num_vars; ++ivar)
        {
            diffArray[ivar][direction][direction] =
                diff2<data_t>(current_cell.get_box_pointers().m_in_ptr[ivar],
                              in_index, stride);
        }
    }

    template <class data_t>
    ALWAYS_INLINE data_t mixed_diff2(const double *in_ptr, const int idx,
                                     const int stride1, const int stride2) const
    {
        auto in = SIMDIFY<data_t>(in_ptr);

        data_t weight_far_far = 6.94444444444444444444e-3;
        data_t weight_near_far = 5.55555555555555555556e-2;
        data_t weight_near_near = 4.44444444444444444444e-1;

        return (weight_far_far * in[idx - 2 * stride1 - 2 * stride2] -
                weight_near_far * in[idx - 2 * stride1 - stride2] +
                weight_near_far * in[idx - 2 * stride1 + stride2] -
                weight_far_far * in[idx - 2 * stride1 + 2 * stride2]

                - weight_near_far * in[idx - stride1 - 2 * stride2] +
                weight_near_near * in[idx - stride1 - stride2] -
                weight_near_near * in[idx - stride1 + stride2] +
                weight_near_far * in[idx - stride1 + 2 * stride2]

                + weight_near_far * in[idx + stride1 - 2 * stride2] -
                weight_near_near * in[idx + stride1 - stride2] +
                weight_near_near * in[idx + stride1 + stride2] -
                weight_near_far * in[idx + stride1 + 2 * stride2]

                - weight_far_far * in[idx + 2 * stride1 - 2 * stride2] +
                weight_near_far * in[idx + 2 * stride1 - stride2] -
                weight_near_far * in[idx + 2 * stride1 + stride2] +
                weight_far_far * in[idx + 2 * stride1 + 2 * stride2]) *
               m_one_over_dx2;
    }

    template <class data_t, template <typename> class vars_t>
    void mixed_diff2(vars_t<Tensor<2, data_t>> &d2,
                     const Cell<data_t> &current_cell, int direction1,
                     int direction2) const
    {
        const int stride1 =
            current_cell.get_box_pointers().m_in_stride[direction1];
        const int stride2 =
            current_cell.get_box_pointers().m_in_stride[direction2];
        const int in_index = current_cell.get_in_index();
        d2.enum_mapping(
            [&](const int &ivar, Tensor<2, data_t> &var)
            {
                auto tmp = mixed_diff2<data_t>(
                    current_cell.get_box_pointers().m_in_ptr[ivar], in_index,
                    stride1, stride2);
                var[direction1][direction2] = tmp;
                var[direction2][direction1] = tmp;
            });
    }

    template <class data_t, int num_vars>
    void mixed_diff2(Tensor<2, data_t> (&diffArray)[num_vars],
                     const Cell<data_t> &current_cell, int direction1,
                     int direction2) const
    {
        const int stride1 =
            current_cell.get_box_pointers().m_in_stride[direction1];
        const int stride2 =
            current_cell.get_box_pointers().m_in_stride[direction2];
        const int in_index = current_cell.get_in_index();
        for (int ivar = 0; ivar < num_vars; ++ivar)
        {
            data_t diff2_value = mixed_diff2<data_t>(
                current_cell.get_box_pointers().m_in_ptr[ivar], in_index,
                stride1, stride2);
            diffArray[ivar][direction1][direction2] = diff2_value;
            diffArray[ivar][direction2][direction1] = diff2_value;
        }
    }

    template <typename data_t> struct Detector;

    /// Calculates all second derivatives and returns as variable type specified
    /// by the template parameter
    template <template <typename> class vars_t, class data_t>
    auto diff2(const Cell<data_t> &current_cell) const
    {
        vars_t<Tensor<2, data_t>> d2;
        const auto in_index = current_cell.get_in_index();
        const auto strides = current_cell.get_box_pointers().m_in_stride;
        d2.enum_mapping(
            [&](const int &ivar, Tensor<2, data_t> &var)
            {
                FOR(dir1) // First calculate the repeated derivatives
                {
                    var[dir1][dir1] = diff2<data_t>(
                        current_cell.get_box_pointers().m_in_ptr[ivar],
                        in_index, strides[dir1]);
                    for (int dir2 = 0; dir2 < dir1; ++dir2)
                    {
                        auto tmp = mixed_diff2<data_t>(
                            current_cell.get_box_pointers().m_in_ptr[ivar],
                            in_index, strides[dir1], strides[dir2]);
                        var[dir1][dir2] = tmp;
                        var[dir2][dir1] = tmp;
                    }
                }
            });
        return d2;
    }

  protected: // Let's keep this protected ... we may want to change the
             // advection calculation
    template <class data_t, class mask_t>
    ALWAYS_INLINE data_t advection_term(const double *in_ptr, const int idx,
                                        const data_t &vec_comp,
                                        const int stride,
                                        const mask_t shift_positive) const
    {
        const auto in = SIMDIFY<data_t>(in_ptr);
        const data_t in_left = in[idx - stride];
        const data_t in_centre = in[idx];
        const data_t in_right = in[idx + stride];

        data_t weight_0 = -2.50000000000000000000e-1;
        data_t weight_1 = -8.33333333333333333333e-1;
        data_t weight_2 = +1.50000000000000000000e+0;
        data_t weight_3 = -5.00000000000000000000e-1;
        data_t weight_4 = +8.33333333333333333333e-2;

        data_t upwind;
        upwind = vec_comp *
                 (weight_0 * in_left + weight_1 * in_centre +
                  weight_2 * in_right + weight_3 * in[idx + 2 * stride] +
                  weight_4 * in[idx + 3 * stride]) *
                 m_one_over_dx;

        data_t downwind;
        downwind = vec_comp *
                   (-weight_4 * in[idx - 3 * stride] -
                    weight_3 * in[idx - 2 * stride] - weight_2 * in_left -
                    weight_1 * in_centre - weight_0 * in_right) *
                   m_one_over_dx;

        return simd_conditional(shift_positive, upwind, downwind);
    }

  public:
    template <class data_t, template <typename> class vars_t>
    void add_advection(vars_t<data_t> &vars, const Cell<data_t> &current_cell,
                       const data_t &vec_comp, const int dir) const
    {
        const int stride = current_cell.get_box_pointers().m_in_stride[dir];
        auto shift_positive = simd_compare_gt(vec_comp, 0.0);
        const int in_index = current_cell.get_in_index();
        vars.enum_mapping(
            [&](const int &ivar, data_t &var)
            {
                var += advection_term(
                    current_cell.get_box_pointers().m_in_ptr[ivar], in_index,
                    vec_comp, stride, shift_positive);
            });
    }

    template <class data_t, int num_vars>
    void add_advection(data_t (&out)[num_vars],
                       const Cell<data_t> &current_cell, const data_t &vec_comp,
                       const int dir) const
    {
        const int stride = current_cell.get_box_pointers().m_in_stride[dir];
        auto shift_positive = simd_compare_gt(vec_comp, 0.0);
        const int in_index = current_cell.get_in_index();
        for (int ivar = 0; ivar < num_vars; ++ivar)
        {
            out[ivar] +=
                advection_term(current_cell.get_box_pointers().m_in_ptr[ivar],
                               in_index, vec_comp, stride, shift_positive);
        }
    }

    /// Calculates all second derivatives and returns as variable type specified
    /// by the template parameter
    template <template <typename> class vars_t, class data_t>
    auto advection(const Cell<data_t> &current_cell,
                   const Tensor<1, data_t> &vector) const
    {
        const auto in_index = current_cell.get_in_index();
        const auto strides = current_cell.get_box_pointers().m_in_stride;
        vars_t<data_t> advec;
        advec.enum_mapping(
            [&](const int &ivar, data_t &var)
            {
                var = 0.;
                FOR(dir)
                {
                    const auto shift_positive =
                        simd_compare_gt(vector[dir], 0.0);
                    var += advection_term(
                        current_cell.get_box_pointers().m_in_ptr[ivar],
                        in_index, vector[dir], strides[dir], shift_positive);
                }
            });
        return advec;
    }

    template <class data_t>
    ALWAYS_INLINE data_t dissipation_term(const double *in_ptr, const int idx,
                                          const int stride) const
    {
        const auto in = SIMDIFY<data_t>(in_ptr);
        data_t weight_vfar = 1.56250e-2;
        data_t weight_far = 9.37500e-2;
        data_t weight_near = 2.34375e-1;
        data_t weight_local = 3.12500e-1;

        return (weight_vfar * in[idx - 3 * stride] -
                weight_far * in[idx - 2 * stride] +
                weight_near * in[idx - stride] - weight_local * in[idx] +
                weight_near * in[idx + stride] -
                weight_far * in[idx + 2 * stride] +
                weight_vfar * in[idx + 3 * stride]) *
               m_one_over_dx;
    }

    template <class data_t, template <typename> class vars_t>
    void add_dissipation(vars_t<data_t> &vars, const Cell<data_t> &current_cell,
                         const double factor, const int direction) const
    {
        const int stride =
            current_cell.get_box_pointers().m_in_stride[direction];
        const int in_index = current_cell.get_in_index();
        vars.enum_mapping(
            [&](const int &ivar, data_t &var)
            {
                var +=
                    factor * dissipation_term<data_t>(
                                 current_cell.get_box_pointers().m_in_ptr[ivar],
                                 in_index, stride);
            });
    }

    template <class data_t, template <typename> class vars_t>
    void add_dissipation(vars_t<data_t> &vars, const Cell<data_t> &current_cell,
                         const double factor) const
    {
        const auto in_index = current_cell.get_in_index();
        vars.enum_mapping(
            [&](const int &ivar, data_t &var)
            {
                FOR(dir)
                {
                    const auto stride =
                        current_cell.get_box_pointers().m_in_stride[dir];
                    var += factor *
                           dissipation_term<data_t>(
                               current_cell.get_box_pointers().m_in_ptr[ivar],
                               in_index, stride);
                }
            });
    }

    template <class data_t, int num_vars>
    void add_dissipation(data_t (&out)[num_vars],
                         const Cell<data_t> &current_cell, const double factor,
                         const int direction) const
    {
        const int stride =
            current_cell.get_box_pointers().m_in_stride[direction];
        const int in_index = current_cell.get_in_index();
        for (int ivar = 0; ivar < num_vars; ++ivar)
        {
            out[ivar] +=
                factor * dissipation_term<data_t>(
                             current_cell.get_box_pointers().m_in_ptr[ivar],
                             in_index, stride);
        }
    }
};

#endif /* FOURTHORDERDERIVATIVES_HPP_ */
