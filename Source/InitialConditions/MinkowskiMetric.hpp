/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MINKOWSKIMETRIC_HPP_
#define MINKOWSKIMETRIC_HPP_

#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs c_NUM
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which sets CCZ4 vars to Minkowski
class MinkowskiMetric
{
  public:
    MinkowskiMetric() {}

  public:
    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // Set Minkowski vars
        CCZ4Vars::VarsWithGauge<data_t> vars;
        VarsTools::assign(vars, 0.0);

        vars.lapse = 1.0;
        vars.chi = 1.0;
        FOR1(i) { vars.h[i][i] = 1.0; }

        // store vars
        current_cell.store_vars(vars);
    }
};

#endif /* MINKOWSKIMETRIC_HPP_ */
