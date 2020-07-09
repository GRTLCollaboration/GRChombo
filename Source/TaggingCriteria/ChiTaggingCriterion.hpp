/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHITAGGINGCRITERION_HPP_
#define CHITAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

class ChiTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;

    /// Vars object for chi
    template <class data_t> struct Vars
    {
        data_t chi; //!< Conformal factor

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
        }
    };

  public:
    ChiTaggingCriterion(double dx) : m_dx(dx), m_deriv(dx){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        //        auto chi = current_cell.load_vars(c_chi);
        //        Tensor<1, data_t> d1_chi;
        //        FOR1(idir) m_deriv.diff1(d1_chi, current_cell, idir, c_chi);

        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Vars>(current_cell);

        data_t mod_d2_chi = 0;
        FOR2(idir, jdir)
        {
            mod_d2_chi += d2.chi[idir][jdir] * d2.chi[idir][jdir] /
                          (1e-2 + abs(d1.chi[idir] * d1.chi[jdir]));
        }

        data_t criterion = m_dx * sqrt(mod_d2_chi);

        // data_t mod_d1_chi = 0;
        // FOR1(idir) mod_d1_chi += d1_chi[idir] * d1_chi[idir];
        // data_t criterion = m_dx * sqrt(mod_d1_chi) / pow(chi, 2);

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* CHITAGGINGCRITERION_HPP_ */
