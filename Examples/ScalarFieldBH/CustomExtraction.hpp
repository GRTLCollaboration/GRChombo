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
#include "SmallDataIO.hpp"
#include "SphericalHarmonics.hpp"
#include "UserVariables.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

// Function to convert double to string with precision
std::string to_string_with_precision(double value, int precision)
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

//!  The class allows extraction of the values of components at
//!  specified points

//!  In this example, we do a lineout data extraction along x-axis

class CustomExtraction
{
  private:
    //! Params for extraction
    const int m_comp1;
    const int m_comp2;
    const int m_num_points;
    const double m_L;
    const std::array<double, CH_SPACEDIM>
        m_origin; // Origin of an extraction line
    const double m_dt;
    const double m_time;

  public:
    //! The constructor
    CustomExtraction(int a_comp1, int a_comp2, int a_num_points, double a_L,
                     std::array<double, CH_SPACEDIM> a_origin, double a_dt,
                     double a_time)
        : m_comp1(a_comp1), m_comp2(a_comp2), m_num_points(a_num_points),
          m_origin(a_origin), m_L(a_L), m_dt(a_dt), m_time(a_time)
    {
    }

    //! Destructor
    ~CustomExtraction() {}

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<2>> *a_interpolator,
                       std::string a_file_prefix) const
    {
        CH_TIME("CustomExtraction::execute_query");
        if (a_interpolator == nullptr)
        {
            MayDay::Error("Interpolator has not been initialised.");
        }
        std::vector<double> interp_ham_data(m_num_points);
        std::vector<double> interp_mom_data(m_num_points);
        std::vector<double> interp_x(m_num_points);
        std::vector<double> interp_y(m_num_points);
        std::vector<double> interp_z(m_num_points);

        // Work out the coordinates
        // go out along x-axis from 0 to L
        for (int idx = 0; idx < m_num_points; ++idx)
        {
            interp_x[idx] = (double(idx) / double(m_num_points) * m_L);
            interp_y[idx] = m_origin[1];
            interp_z[idx] = m_origin[2];
        }

        // set up the query
        InterpolationQuery query_ham(m_num_points);
        query_ham.setCoords(0, interp_x.data())
            .setCoords(1, interp_y.data())
            .setCoords(2, interp_z.data())
            .addComp(m_comp1, interp_ham_data.data(), Derivative::LOCAL,
                     VariableType::diagnostic); // evolution or diagnostic
        InterpolationQuery query_mom(m_num_points);
        query_mom.setCoords(0, interp_x.data())
            .setCoords(1, interp_y.data())
            .setCoords(2, interp_z.data())
            .addComp(m_comp2, interp_mom_data.data(), Derivative::LOCAL,
                     VariableType::diagnostic); // evolution or diagnostic

        // submit the query
        a_interpolator->interp(query_ham);
        a_interpolator->interp(query_mom);

        // now write out
        bool first_step = (m_time == 0.0);
        double restart_time = 0.0;
        SmallDataIO output_file(a_file_prefix, m_dt, m_time, restart_time,
                                SmallDataIO::APPEND, first_step);

        std::vector<std::string> header_line(2);
        if (first_step)
        {
            header_line[0] = "Ham";
            header_line[1] = "Mom";
            output_file.write_header_line(header_line, "x");
        }

        for (int idx = 0; idx < m_num_points; ++idx)
        {
            std::vector<double> data(2);
            data[0] = interp_ham_data[idx];
            data[1] = interp_mom_data[idx];
            output_file.write_data_line(data, interp_x[idx]);
        }
    }
};

#endif /* CUSTOMEXTRACTION_HPP_ */
