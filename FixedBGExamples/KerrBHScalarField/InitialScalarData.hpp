/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

#include "FixedBGScalarField.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates a constant scalar field given params for initial
//! matter config
class InitialScalarData
{
  public:
    struct params_t
    {
        double mass;
        double amplitude;
    };

    //! The constructor for the class
    InitialScalarData(params_t a_params) : m_params(a_params) {}

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

        FixedBGScalarField<>::Vars<data_t> vars;
        VarsTools::assign(vars, 0.);

        // set the field vars
        vars.phi = m_params.amplitude;
        vars.Pi = 0;

        current_cell.store_vars(vars);
    }

  protected:
    const params_t m_params;
};

#endif /* INITIALSCALARDATA_HPP_ */
