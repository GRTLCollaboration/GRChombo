/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "SPMD.H" // for Chombo_MPI

// Other includes
#include "SmallDataIO.hpp"
#include <cmath>
#include <random>
// (MR): if it were up to me, I'd be using the C++17 filesystems library
// instead of cstdio but I'm sure someone would tell me off for not maintaining
// backwards compatability.
#include <cstdio> // for std::rename and std::remove
#include <iomanip>
#include <iostream>
#include <sstream>

// Chombo namespace
#include "UsingNamespace.H"

// ------------ Constructors -----------------

// This has to be initialised outside the class declaration in C++14
const std::string SmallDataIO::s_default_file_extension = ".dat";

SmallDataIO::SmallDataIO(std::string a_filename_prefix, double a_dt,
                         double a_time, double a_restart_time, Mode a_mode,
                         bool a_first_step, std::string a_file_extension,
                         int a_data_precision, int a_coords_precision,
                         int a_filename_steps_width)
    : m_filename(a_filename_prefix + a_file_extension), m_dt(a_dt),
      m_time(a_time), m_restart_time(a_restart_time), m_mode(a_mode),
      m_first_step(a_first_step), m_data_precision(a_data_precision),
      // data columns need extra space for scientific notation
      // compared to coords columns
      m_data_width(m_data_precision + 10),
      m_data_epsilon(std::pow(10.0, -a_data_precision)),
      m_coords_precision(a_coords_precision),
      m_coords_width(m_coords_precision + 5),
      m_coords_epsilon(std::pow(10.0, -a_coords_precision))
{
    if (procID() == 0)
    {
        if (m_coords_precision >= m_data_precision)
        {
            MayDay::Error("SmallDataIO: m_coords_precision should be smaller "
                          "than m_data_precision for SmallDataIOReader");
        }
        std::ios::openmode file_openmode;
        if (m_mode == APPEND)
        {
            if (m_first_step)
            {
                // overwrite any existing file if this is the first step
                file_openmode = std::ios::out;
            }
            else if (m_restart_time >= 0. &&
                     m_time < m_restart_time + m_dt + m_coords_epsilon)
            {
                // allow reading in the restart case so that duplicate time
                // data may be removed
                file_openmode = std::ios::app | std::ios::in;
            }
            else
            {
                // default mode is just appending to existing file
                file_openmode = std::ios::app;
            }
        }
        else if (m_mode == NEW)
        {
            file_openmode = std::ios::out;
            m_filename =
                get_new_filename(a_filename_prefix, m_dt, m_time,
                                 a_file_extension, a_filename_steps_width);
        }
        else if (m_mode == READ)
        {
            file_openmode = std::ios::in;
        }
        else
        {
            MayDay::Error("SmallDataIO: mode not supported");
        }
        m_file.open(m_filename, file_openmode);
        if (!m_file)
        {
            MayDay::Error((std::string("SmallDataIO::error opening file '") +
                           m_filename + "' for writing")
                              .c_str());
        }
    }
}

SmallDataIO::SmallDataIO(std::string a_filename_prefix, double a_dt,
                         double a_time, double a_restart_time, Mode a_mode,
                         std::string a_file_extension, int a_data_precision,
                         int a_coords_precision, int a_filename_steps_width)
    : SmallDataIO(a_filename_prefix, a_dt, a_time, a_restart_time, a_mode,
                  (a_time == a_dt), a_file_extension, a_data_precision,
                  a_coords_precision, a_filename_steps_width)
{
}

SmallDataIO::SmallDataIO(std::string a_filename_prefix,
                         std::string a_file_extension, int a_data_precision,
                         int a_coords_precision)
    : SmallDataIO(a_filename_prefix, 0.0, 0.0, 0.0, READ, false,
                  a_file_extension, a_data_precision, a_coords_precision, 0)
{
}

//! Destructor (closes file)
SmallDataIO::~SmallDataIO()
{
    if (procID() == 0)
    {
        m_file.close();
    }
}

// ------------ Writing Functions ------------

void SmallDataIO::write_header_line(
    const std::vector<std::string> &a_header_strings,
    const std::string &a_pre_header_string)
{
    std::vector<std::string> pre_header_strings;
    if (a_pre_header_string != "")
    {
        pre_header_strings.push_back(a_pre_header_string);
    }
    write_header_line(a_header_strings, pre_header_strings);
}

void SmallDataIO::write_header_line(
    const std::vector<std::string> &a_header_strings,
    const std::vector<std::string> &a_pre_header_strings)
{
    if (procID() == 0)
    {
        // all header lines start with a '#'.
        m_file << "#";
        for (int istr = 0; istr < a_pre_header_strings.size(); ++istr)
        {
            // first column header is shorter due to preceeding #
            if (istr == 0)
            {
                m_file << std::setw(m_coords_width - 1)
                       << a_pre_header_strings[istr];
            }
            else
            {
                m_file << std::setw(m_coords_width)
                       << a_pre_header_strings[istr];
            }
        }
        for (std::string header_item : a_header_strings)
        {
            m_file << std::setw(m_data_width) << header_item;
        }
        m_file << "\n";
    }
}

void SmallDataIO::write_data_line(const std::vector<double> &a_data,
                                  const double a_coord)
{
    const std::vector<double> coords(1, a_coord);
    write_data_line(a_data, coords);
}

void SmallDataIO::write_time_data_line(const std::vector<double> &a_data)
{
    write_data_line(a_data, m_time);
}

void SmallDataIO::write_data_line(const std::vector<double> &a_data,
                                  const std::vector<double> &a_coords)
{
    if (procID() == 0)
    {
        m_file << std::fixed << std::setprecision(m_coords_precision);
        for (double coord : a_coords)
        {
            m_file << std::setw(m_coords_width) << coord;
        }
        m_file << std::scientific << std::setprecision(m_data_precision);
        for (double data : a_data)
        {
            m_file << std::setw(m_data_width) << data;
        }
        m_file << "\n";
    }
}

void SmallDataIO::line_break()
{
    if (procID() == 0)
    {
        m_file << "\n\n";
    }
}

void SmallDataIO::remove_duplicate_time_data(const bool keep_m_time_data)
{
    if (procID() == 0 && m_restart_time >= 0. && m_mode == APPEND &&
        m_time < m_restart_time + m_dt + m_coords_epsilon)
    {
        // copy lines with time < m_time into a temporary file
        m_file.seekg(0);
        std::string line;
        // adding a random integer might make this a little more robust...
        const int random_int = std::default_random_engine()();
        std::string temp_filename =
            m_filename + ".temp" + std::to_string(random_int);
        std::ofstream temp_file(temp_filename);
        int sign = -1;
        if (keep_m_time_data)
        {
            sign = 1;
        }
        while (std::getline(m_file, line))
        {
            if (!(line.find("#") == std::string::npos))
            {
                temp_file << line << "\n";
            }
            else if (std::stod(line.substr(0, m_coords_width)) <
                     m_time + sign * m_coords_epsilon)
            {
                temp_file << line << "\n";
            }
        }

        m_file.close();
        temp_file.close();

        // now delete the original file and rename the temporary file with the
        // original filename
        std::remove(m_filename.data());
        std::rename(temp_filename.data(), m_filename.data());
        // reopen the file in append mode
        m_file.open(m_filename, std::ios::app);
    }
}

// ------------ Reading Functions ------------

void SmallDataIO::get_specific_data_line(std::vector<double> &a_out_data,
                                         const std::vector<double> a_coords)
{
    if (procID() == 0)
    {
        bool line_found = false;
        // first set the current position to the beginning of the file
        m_file.seekg(0);

        // get a string of the coords as they are in the file
        std::stringstream coords_ss;
        coords_ss << std::fixed << std::setprecision(m_coords_precision);
        for (double coord : a_coords)
        {
            coords_ss << std::setw(m_coords_width) << coord;
        }
        std::string coords_string = coords_ss.str();

        // now search for lines that start with coords_string and put the data
        // in a_out_data
        std::string line;
        while (std::getline(m_file, line))
        {
            if (!(line.find(coords_string) == std::string::npos))
            {
                for (int ichar = a_coords.size() * m_coords_width;
                     ichar < line.size(); ichar += m_data_width)
                {
                    double data_value =
                        std::stod(line.substr(ichar, m_data_width));
                    a_out_data.push_back(data_value);
                }
                line_found = true;
                // only want the first occurrence so break the while loop
                break;
            }
        }
        if (!line_found)
        {
            MayDay::Error(
                "SmallDataIO : Data to be read in at coord not found in file");
        }
    }
    // now broadcast the vector to all ranks using Chombo broadcast function
    // need to convert std::vector to Vector first
    Vector<double> data_Vect(a_out_data);
    int broadcast_rank = 0;
    broadcast(data_Vect, broadcast_rank);
    a_out_data = data_Vect.stdVector();
}

void SmallDataIO::get_specific_data_line(std::vector<double> &a_out_data,
                                         const double a_coord)
{
    std::vector<double> coords(1, a_coord);
    get_specific_data_line(a_out_data, coords);
}

// ------------ Other Functions --------------

std::string SmallDataIO::get_new_filename(const std::string &a_file_prefix,
                                          double a_dt, double a_time,
                                          const std::string &a_file_extension,
                                          int a_filename_steps_width)
{
    CH_assert(a_dt > 0);
    const int step = std::round(a_time / a_dt);

    // append step number to filename if in NEW mode
    return a_file_prefix + pad_number(step) + a_file_extension;
}

std::string SmallDataIO::pad_number(int step, int a_filename_steps_width)
{
    // append step number to filename (pad to make it
    // a_filename_steps_width digits).
    std::string step_string = std::to_string(step);
    if (a_filename_steps_width < step_string.length())
    {
        MayDay::Error("SmallDataIO: a_filename_steps_width too small "
                      "for step number");
    }
    std::string step_string_padded =
        std::string(a_filename_steps_width - step_string.length(), '0') +
        step_string;
    return step_string_padded;
}

// returns m_data_epsilon
double SmallDataIO::get_data_epsilon() const { return m_data_epsilon; }

// returns the default data_epsilon
double SmallDataIO::get_default_data_epsilon()
{
    return pow(10.0, -s_default_data_precision);
}

// returns m_coords_epsilon
double SmallDataIO::get_coords_epsilon() const { return m_coords_epsilon; }

// returns the default coords epsilon
double SmallDataIO::get_default_coords_epsilon()
{
    return pow(10.0, -s_default_coords_precision);
}
