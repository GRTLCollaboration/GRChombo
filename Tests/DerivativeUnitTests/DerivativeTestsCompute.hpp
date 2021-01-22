/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DERIVATIVETESTSCOMPUTE_HPP_
#define DERIVATIVETESTSCOMPUTE_HPP_

#include "VarsTools.hpp"

template <class deriv_t> class DerivativeTestsCompute
{
  protected:
    const deriv_t m_deriv;

  public:
    DerivativeTestsCompute(double dx) : m_deriv(dx) {}

    template <class data_t> struct Vars
    {
        data_t d1;
        data_t d2;
        data_t d2_mixed;
        data_t diss;
        data_t advec_up;
        data_t advec_down;

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_d1, d1);
            define_enum_mapping(mapping_function, c_d2, d2);
            define_enum_mapping(mapping_function, c_d2_mixed, d2_mixed);
            define_enum_mapping(mapping_function, c_diss, diss);
            define_enum_mapping(mapping_function, c_advec_up, advec_up);
            define_enum_mapping(mapping_function, c_advec_down, advec_down);
        }
    };

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto out_d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto out_d2 = m_deriv.template diff2<Vars>(current_cell);

        Vars<data_t> out_diss;
        VarsTools::assign(out_diss, 0.);
        m_deriv.add_dissipation(out_diss, current_cell, 1.0);

        Tensor<1, data_t> shift_down = {-2., 0., -3.};
        const auto out_advec_down =
            m_deriv.template advection<Vars>(current_cell, shift_down);

        Tensor<1, data_t> shift_up = {2., 0., 3.};
        const auto out_advec_up =
            m_deriv.template advection<Vars>(current_cell, shift_up);

        current_cell.store_vars(out_d1.d1[2], c_d1);
        current_cell.store_vars(out_d2.d2[2][2], c_d2);
        current_cell.store_vars(out_d2.d2[0][2], c_d2_mixed);
        current_cell.store_vars(out_diss.diss, c_diss);
        current_cell.store_vars(out_advec_down.advec_down, c_advec_down);
        current_cell.store_vars(out_advec_up.advec_up, c_advec_up);
    }
};

#endif /* DERIVATIVETESTSCOMPUTE_HPP_ */
