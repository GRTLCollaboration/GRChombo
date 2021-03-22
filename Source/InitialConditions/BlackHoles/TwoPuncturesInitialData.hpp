/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef USE_TWOPUNCTURES

#ifndef TWOPUNCTURESINITIALDATA_HPP_
#define TWOPUNCTURESINITIALDATA_HPP_

#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "TwoPunctures.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"
#include <array>

//! This compute class sets the initial data computed by TwoPunctures on the
//! grid
class TwoPuncturesInitialData
{
  protected:
    double m_dx;
    std::array<double, CH_SPACEDIM> m_center;
    const TP::TwoPunctures &m_two_punctures;

  public:
    template <class data_t> using Vars = CCZ4Vars::VarsWithGauge<data_t>;

    TwoPuncturesInitialData(const double a_dx,
                            const std::array<double, CH_SPACEDIM> a_center,
                            const TP::TwoPunctures &a_two_punctures)
        : m_dx(a_dx), m_center(a_center), m_two_punctures(a_two_punctures)
    {
    }

    void compute(Cell<double> current_cell) const;

  protected:
    void interpolate_tp_vars(const Coordinates<double> &coords,
                             Tensor<2, double> &out_h_phys,
                             Tensor<2, double> &out_extrinsic_K,
                             double &out_lapse, Tensor<1, double> &out_shift,
                             double &out_Theta,
                             Tensor<1, double> &out_Z3) const;
};

#include "TwoPuncturesInitialData.impl.hpp"

#endif /* TWOPUNCTURESINITIALDATA_HPP_ */
#endif /* USE_TWOPUNCTURES */
