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

#include "IntVect.H"
#include "MayDay.H"
#include <fstream>

//! Class which sets the initial scalar field matter config
class InitialScalarData
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double amplitude; //!< Amplitude of bump in initial SF bubble
        double velocity;
        std::array<double, CH_SPACEDIM>
            center;   //!< Centre of perturbation in initial SF bubble
        double width; //!< Width of bump in initial SF bubble
        double mass;
        int N_init;
    };

    //! The constructor
    InitialScalarData(params_t a_params, double a_dx, std::vector<std::vector<double> > a_h, double a_hdot)
        : m_params(a_params), m_dx(a_dx), m_h(a_h), m_hdot(a_hdot)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center); // note: coords.x, etc. are in program units
        auto current_cell_index = current_cell.get_in_index(); // pulls the unitless coordinate, or index

        data_t rr = coords.get_radius();
        data_t rr2 = rr * rr;

        // Pull out the grid parametersß
        int N = m_params.N_init;
        double L = N*m_dx;

        // Coordinate of this cell in program units
        data_t x = coords.x + L/2;
        double y = coords.y + L/2;
        double z = coords.z + L/2;

        // Coordinates of this cell, unitless
        int i = static_cast<int>(x / m_dx);
        int j = static_cast<int>(y / m_dx);
        int k = static_cast<int>(z / m_dx);

        // This is to guard against ghost cells that can take you outside 
        // the domain of dependence of the box. Uses periodic BCs.
        if(i < 0)
        {
            i = N + i;
        }
        else if(i >= N)
        {
            i = i - N;
        }

        if(j < 0)
        {
            j = N + j;
        }
        else if(j >= N)
        {
            j = j - N;
        }

        if(k < 0)
        {
            k = N + k;
        }
        else if(k >= N)
        {
            k = k - N;
        }

        // The flattened position (leading with z?)
        int r = k + N*(j + N*i);

        if (current_cell_index < 0)
        {
            cout << current_cell_index << endl;
            MayDay::Error("Cell index value below zero.");
        }
        else if(current_cell_index > pow(m_params.N_init, 3.))
        {
            cout << current_cell_index << endl;
            MayDay::Error("Cell index greater than resolution^3 at coarsest level.");
        }

        // calculate and store the scalar field value
        const data_t phi = m_params.amplitude;
        const data_t phidot = m_params.velocity;

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

        //store tensor metric variables, g_ij = delta_ij + 1/2 h_ij
        current_cell.store_vars(1. + 0.5*m_h[r][0], c_h11);
        current_cell.store_vars(0.5*m_h[r][1], c_h12);
        current_cell.store_vars(0.5*m_h[r][2], c_h13);
        current_cell.store_vars(1. + 0.5*m_h[r][3], c_h22);
        current_cell.store_vars(0.5*m_h[r][4], c_h23);
        current_cell.store_vars(1. + 0.5*m_h[r][5], c_h33);

        current_cell.store_vars(-m_hdot, c_A11);
        current_cell.store_vars(-m_hdot, c_A12);
        current_cell.store_vars(-m_hdot, c_A13);
        current_cell.store_vars(-m_hdot, c_A22);
        current_cell.store_vars(-m_hdot, c_A23);
        current_cell.store_vars(-m_hdot, c_A33);
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
    std::vector< std::vector<double>> m_h;
    double m_hdot;
};

#endif /* INITIALSCALARDATA_HPP_ */
