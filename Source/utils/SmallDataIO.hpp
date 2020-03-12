/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SMALLDATAIO_HPP_
#define SMALLDATAIO_HPP_

// (MR): if it were up to me, I'd be using the C++17 filesystems library
// instead of cstdio but I'm sure someone would tell me off for not maintaining
// backwards compatability.
#include "SPMD.H" // for Chombo_MPI
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

//! A class for reading and writing small data to a file in ASCII format.
/*!
    A class for reading and writing small data, usually 0D, 1D or 2D, in ASCII
    format. For an example on how to use it, see the WeylExtraction class.
*/
class SmallDataIO
{
  public:
    //! Choose between appending data to the same file, writing to a new file
    //! at each timestep or reading a file.
    enum Mode
    {
        APPEND, // data is APPENDed to the same file at each timestep
        NEW,    // data is written to a NEW file at each timestep
        READ    // read data
    };

  protected:
    std::string m_filename; // this is the filename prefix in the NEW Mode
    const double m_dt;
    const double m_time;
    const double m_restart_time;
    int m_step;
    const Mode m_mode;
    const bool m_first_step; // this should be set to true if this is the first
                             // timestep
    const std::string m_file_extension;
    const int m_data_precision;
    const int m_data_width;
    const int m_coords_precision;
    const int m_coords_width;

    std::fstream m_file;
    int m_rank; // only rank 0 does the write out

  public:
    //! Constructor (opens file)
    SmallDataIO(std::string a_filename, double a_dt, double a_time,
                double a_restart_time, Mode a_mode, bool a_first_step,
                std::string a_file_extension = ".dat",
                int a_data_precision = 10, int a_coords_precision = 7)
        : m_filename(a_filename), m_dt(a_dt), m_time(a_time),
          m_restart_time(a_restart_time), m_mode(a_mode),
          m_first_step(a_first_step), m_file_extension(a_file_extension),
          m_data_precision(a_data_precision),
          // data columns need extra space for scientific notation
          // compared to coords columns
          m_data_width(m_data_precision + 10),
          m_coords_precision(a_coords_precision),
          m_coords_width(m_coords_precision + 5)
    {
        constexpr double epsilon = 1.0e-8;
#ifdef CH_MPI
        MPI_Comm_rank(Chombo_MPI::comm, &m_rank);
#else
        m_rank = 0;
#endif
        if (m_rank == 0)
        {
            std::ios::openmode file_openmode;
            if (m_mode == APPEND)
            {
                if (m_first_step)
                {
                    // overwrite any existing file if this is the first step
                    file_openmode = std::ios::out;
                }
                else if (m_restart_time > 0. &&
                         m_time < m_restart_time + m_dt + epsilon)
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
                m_step = std::round(m_time / m_dt);
                // append step number to filename if in NEW mode
                m_filename += std::to_string(m_step);
            }
            else if (m_mode == READ)
            {
                file_openmode = std::ios::in;
            }
            else
            {
                MayDay::Error("SmallDataIO: mode not supported");
            }
            m_filename += m_file_extension;
            m_file.open(m_filename, file_openmode);
            if (!m_file)
            {
                MayDay::Error("SmallDataIO::error opening file for writing");
            }
        }
    }

    //! Old constructor which assumes SmallDataIO is called in
    //! specificPostTimeStep
    SmallDataIO(std::string a_filename, double a_dt, double a_time,
                double a_restart_time, Mode a_mode,
                std::string a_file_extension = ".dat",
                int a_data_precision = 10, int a_coords_precision = 7)
        : SmallDataIO(a_filename, a_dt, a_time, a_restart_time, a_mode,
                      (a_time == a_dt), a_file_extension, a_data_precision,
                      a_coords_precision)
    {
    }

    //! Constructor for reading when m_time, m_dt, m_restart_time are irrelevant
    SmallDataIO(std::string a_filename, std::string a_file_extension = ".dat",
                int a_data_precision = 10, int a_coords_precision = 7)
        : SmallDataIO(a_filename, 0.0, 0.0, 0.0, READ, false, a_file_extension,
                      a_data_precision, a_coords_precision)
    {
    }

    //! Destructor (closes file)
    ~SmallDataIO()
    {
        if (m_rank == 0)
        {
            m_file.close();
        }
    }

    // disable default copy constructor and assignment operator
    SmallDataIO(const SmallDataIO &) = delete;
    SmallDataIO &operator=(const SmallDataIO &) = delete;

    // ------------ Writing Functions ------------

    //! Writes a header_line
    //! Use this for 0D or 1D data, where the first column is either the time
    //! or another coordinate whose name should be provided in
    //! a_pre_header_string.
    void write_header_line(const std::vector<std::string> &a_header_strings,
                           const std::string &a_pre_header_string = "time");

    //! Writes a header line
    //! Use this for 1D or 2D data when the first two or more columns are
    //! coordinates whose names should be provided in the vector of strings
    //! a_pre_header_strings
    void
    write_header_line(const std::vector<std::string> &a_header_strings,
                      const std::vector<std::string> &a_pre_header_strings);

    //! Writes a data line
    //! Use this for 0D or 1D data, where the first column is either the time or
    //! another coordinate.
    void write_data_line(const std::vector<double> &a_data,
                         const double a_coord);

    //! Writes a data line for a specific time.
    void write_time_data_line(const std::vector<double> &a_data);

    //! Writes a data line
    //! Use this for 1D or 2D data when the first two or more columns are
    //! coordinates.
    void write_data_line(const std::vector<double> &a_data,
                         const std::vector<double> &a_coords);

    //! This just adds a double line break to the file.
    void line_break();

    //! if restarting from an earlier checkpoint file, this function removes
    //! any time data that will be replaced.
    void remove_duplicate_time_data(const bool keep_m_time_data = false);

    // ------------ Reading Functions ------------

    //! Get the data associated to specific coordinates from the file
    //! Note only the first line with the given coordinates is obtained
    void get_specific_data_line(std::vector<double> &a_out_data,
                                const std::vector<double> a_coords);

    //! Get the data associated to a specific coordinate (e.g. time) from the
    //! file
    void get_specific_data_line(std::vector<double> &a_out_data,
                                const double a_coord);
};

#endif /* SMALLDATAIO_HPP_ */
