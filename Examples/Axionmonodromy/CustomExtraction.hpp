/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CUSTOMEXTRACTION_HPP_
#define CUSTOMEXTRACTION_HPP_

#include "AMRInterpolator.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "SimulationParametersBase.hpp"
#include "SphericalHarmonics.hpp"
#include "UserVariables.hpp"
#include "SmallDataIO.hpp"
#include <fstream>
#include <iostream>

//!  The class allows extraction of the values of components at
//!  specified point
class CustomExtraction
{
  private:
    //! Params for extraction
    const int m_comp1 = c_chi;
    const int m_comp2 = c_phi;
    const int m_comp3 = c_lapse;
    const double m_L;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_time;
    const double m_dt;
    const double m_restart_time;

  public:
    //! The constructor
    CustomExtraction(double a_L,
                     std::array<double, CH_SPACEDIM> a_center,
                     double a_time, double a_dt, double a_restart_time = 0.0) :
                    m_center(a_center), m_restart_time(a_restart_time),
                    m_L(a_L), m_time(a_time), m_dt(a_dt)
    {
    }

    //! Destructor
    ~CustomExtraction() {}

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator,
                       std::string a_filename) const
    {
        // extraction data holders
        int num_points = 1;
        int num_comps = 3;
        std::vector<double> interp_x(num_points);
        std::vector<double> interp_y(num_points);
        std::vector<double> interp_z(num_points);
        std::vector<double> interp_comp1(num_points);
        std::vector<double> interp_comp2(num_points);
        std::vector<double> interp_comp3(num_points);

        // extract at the centre
        interp_x[0] = m_center[0];
        interp_x[1] = m_center[1];
        interp_x[2] = m_center[2];

        // setup query
        InterpolationQuery query(num_points);
        query.setCoords(0, interp_x.data())
             .setCoords(1, interp_y.data())
             .setCoords(2, interp_z.data())
             .addComp(m_comp1, interp_comp1.data())
             .addComp(m_comp2, interp_comp2.data())
             .addComp(m_comp3, interp_comp3.data());

        // engage!
        a_interpolator->interp(query);

        // write out
        SmallDataIO extraction_file(a_filename, m_dt, m_time, m_restart_time,
                                           SmallDataIO::APPEND);
        extraction_file.remove_duplicate_time_data();
        if (m_time == m_dt)
        {
            // make header strings if at first timestep
            std::vector<std::string> header1_strings(num_comps*num_points);
            header1_strings[0] = "chi";
            header1_strings[1] = "phi";
            header1_strings[2] = "lapse";
            extraction_file.write_header_line(header1_strings);
        }
        std::vector<double> data_for_writing(num_comps*num_points);
        data_for_writing[0] = interp_comp1[0];
        data_for_writing[1] = interp_comp2[0];
        data_for_writing[2] = interp_comp3[0];
        extraction_file.write_time_data_line(data_for_writing);
    }
};

#endif /* CUSTOMEXTRACTION_HPP_ */
