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
    InitialScalarData(params_t a_params, double a_dx)
        : m_dx(a_dx), m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center); // note: coords.x, etc. are in program units
        auto current_cell_index = current_cell.get_in_index(); // pulls the unitless coordinate, or index

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

        /*array<data_t, 6> h;
        array<data_t, 6> hdot;

        ifstream gw_pos;
        ifstream gw_vel;
        gw_pos.open("./gw-re-position.dat", ios::in); //open the file with the waves in it
        gw_vel.open("./gw-re-velocity.dat", ios::in);

        if (gw_pos.fail() || gw_vel.fail())
        {
            MayDay::Error("GW position or velocity file failed to open.");
        }

        std::string junk;
        std::string hs;
        std::string hdots;
        for(int m=0; m<=current_cell_index; m++)
        {
            junk = "";

            if(m == current_cell_index)
            {
                for(int s=0; s<6; s++)
                {
                    std::getline(gw_pos, hs); //choose that row from the data file and save it
                    std::getline(gw_vel, hdots);
                }
            }

            std::getline(gw_pos, junk);

            if (gw_pos.eof() || gw_vel.eof()) { MayDay::Error("End of file reached, exiting."); }
        }

        gw_pos.close();
        gw_vel.close();

        if (current_cell_index == 100)
        {
            cout << current_cell_index << "," << h[0] << "," << hdot[0] << endl;
            MayDay::Error("100 index reached");
        }*/

        //store tensor metric variables, g_ij = delta_ij + 1/2 h_ij
        current_cell.store_vars(1. + 0.5*m_h[current_cell_index][0], c_h11);
        current_cell.store_vars(0.5*m_h[current_cell_index][1], c_h12);
        current_cell.store_vars(0.5*m_h[current_cell_index][2], c_h13);
        current_cell.store_vars(1. + 0.5*m_h[current_cell_index][3], c_h22);
        current_cell.store_vars(0.5*m_h[current_cell_index][4], c_h23);
        current_cell.store_vars(1. + 0.5*m_h[current_cell_index][5], c_h33);

        /*current_cell.store_vars(-hdot[0], c_A11);
        current_cell.store_vars(-hdot[1], c_A12);
        current_cell.store_vars(-hdot[2], c_A13);
        current_cell.store_vars(-hdot[3], c_A22);
        current_cell.store_vars(-hdot[4], c_A23);
        current_cell.store_vars(-hdot[5], c_A33);*/

        /*current_cell.store_vars(1., c_h11);
        current_cell.store_vars(0., c_h12);
        current_cell.store_vars(0., c_h13);
        current_cell.store_vars(1., c_h22);
        current_cell.store_vars(0., c_h23);
        current_cell.store_vars(1., c_h33);*/

        current_cell.store_vars(0., c_A11);
        current_cell.store_vars(0., c_A12);
        current_cell.store_vars(0., c_A13);
        current_cell.store_vars(0., c_A22);
        current_cell.store_vars(0., c_A23);
        current_cell.store_vars(0., c_A33);
    }

    void load_gws(std::vector<std::vector<double>> h) const
    {
        ifstream gw_pos;
        ifstream gw_vel;
        gw_pos.open("./gw-re-position.dat", ios::in); //open the file with the waves in it
        gw_vel.open("./gw-re-velocity.dat", ios::in);

        if (gw_pos.fail() || gw_vel.fail())
        {
            MayDay::Error("GW position or velocity file failed to open.");
        }

        std::string junk;
        std::string hs;
        std::string hdots;
        int count = 0;

        /*while (!gw_pos.eof())
        {
            for(int s=0; s<6; s++)
            {
                gw_pos >> h[count][s];
            }
            count++;

            if (count == 100)
            {
                cout << count << "," << h[count][0] << endl;
                MayDay::Error("100 index reached");
            }
        }*/
        
        /*for(int m=0; m<=current_cell_index; m++)
        {
            junk = "";

            if(m == current_cell_index)
            {
                for(int s=0; s<6; s++)
                {
                    std::getline(gw_pos, hs); //choose that row from the data file and save it
                    std::getline(gw_vel, hdots);
                }
            }

            std::getline(gw_pos, junk);

            if (gw_pos.eof() || gw_vel.eof()) { MayDay::Error("End of file reached, exiting."); }
        }*/

        gw_pos.close();
        gw_vel.close();
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
    std::vector< std::vector<double>> m_h;
};

#endif /* INITIALSCALARDATA_HPP_ */
