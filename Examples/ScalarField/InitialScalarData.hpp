/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "Potential.hpp"

//! Class which sets the initial scalar field matter config
class InitialScalarData
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double amplitude; //!< Amplitude of bump in initial SF bubble
        std::array<double, CH_SPACEDIM>
            center;   //!< Centre of perturbation in initial SF bubble
        double width; //!< Width of bump in initial SF bubble
        double mass;
    };

    //! The constructor
    InitialScalarData(params_t a_params, double a_dx)
        : m_dx(a_dx), m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

        // calculate and store the scalar field value
        const data_t phi = m_params.amplitude;
        const data_t phidot = -0.00162846;

        current_cell.store_vars(phi, c_phi);
        current_cell.store_vars(phidot, c_Pi);

        //Planck's constant, set to change the units of the scalar field
        double m_pl = 1.0;//1.220890e19; // (GeV)
        data_t V;
        data_t dV;

        //calculate and store gauge variables
        data_t lapse = 1.0;
        data_t shift[3] = {0.0, 0.0, 0.0};

        current_cell.store_vars(lapse, c_lapse);
        current_cell.store_vars(shift[0], c_shift1);
        current_cell.store_vars(shift[1], c_shift2);
        current_cell.store_vars(shift[2], c_shift3);

        //calculate and store scalar metric variables
        data_t chi = 1.0;
        data_t K = -3.0*sqrt((8*M_PI/3/m_pl)*(0.5*phidot*phidot + 0.5*pow(m_params.mass * phi, 2.0))); //This needs grad energy with perturbations...

        current_cell.store_vars(chi, c_chi);
        current_cell.store_vars(K, c_K);

        //store tensor metric variables
        current_cell.store_vars(1.0, c_h11);
        current_cell.store_vars(0.0, c_h12);
        current_cell.store_vars(0.0, c_h13);
        current_cell.store_vars(1.0, c_h22);
        current_cell.store_vars(0.0, c_h23);
        current_cell.store_vars(1.0, c_h33);

        current_cell.store_vars(0.0, c_A11);
        current_cell.store_vars(0.0, c_A12);
        current_cell.store_vars(0.0, c_A13);
        current_cell.store_vars(0.0, c_A22);
        current_cell.store_vars(0.0, c_A23);
        current_cell.store_vars(0.0, c_A33);
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
};

#endif /* INITIALSCALARDATA_HPP_ */
