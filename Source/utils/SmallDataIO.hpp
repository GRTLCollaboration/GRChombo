/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SMALLDATAIO_HPP_
#define SMALLDATAIO_HPP_

#include <fstream>
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
    std::string m_filename;
    const double m_dt;
    const double m_time;
    const double m_restart_time;
    int m_step;
    const Mode m_mode;
    const bool m_first_step; // this should be set to true if this is the first
                             // timestep
    static const std::string s_default_file_extension;
    static constexpr int s_default_data_precision = 10;
    const int m_data_precision;
    const int m_data_width;
    const double m_data_epsilon; //!< the maximum data precision error
    static constexpr int s_default_coords_precision = 7;
    const int m_coords_precision;
    const int m_coords_width;
    const double m_coords_epsilon; //!< the maximum coords precision error
    static constexpr int s_default_filename_steps_width = 6;

    std::fstream m_file;
    int m_rank; // only rank 0 does the write out

  public:
    //! Constructor (opens file)
    SmallDataIO(std::string a_filename_prefix, double a_dt, double a_time,
                double a_restart_time, Mode a_mode, bool a_first_step,
                std::string a_file_extension = s_default_file_extension,
                int a_data_precision = s_default_data_precision,
                int a_coords_precision = s_default_coords_precision,
                int a_filename_steps_width = s_default_filename_steps_width);

    //! Old constructor which assumes SmallDataIO is called in
    //! specificPostTimeStep
    SmallDataIO(std::string a_filename_prefix, double a_dt, double a_time,
                double a_restart_time, Mode a_mode,
                std::string a_file_extension = s_default_file_extension,
                int a_data_precision = s_default_data_precision,
                int a_coords_precision = s_default_coords_precision,
                int a_filename_steps_width = s_default_filename_steps_width);

    //! Constructor for reading when m_time, m_dt, m_restart_time are irrelevant
    SmallDataIO(std::string a_filename_prefix,
                std::string a_file_extension = s_default_file_extension,
                int a_data_precision = s_default_data_precision,
                int a_coords_precision = s_default_coords_precision);

    //! Destructor (closes file)
    ~SmallDataIO();

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
                         const std::vector<double> &a_coords = {});

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

    // ------------ Other Functions --------------

    //! returns the full filename of a file created in NEW mode at time=a_time
    //! with dt=a_dt
    static std::string get_new_filename(
        const std::string &a_filename_prefix, double a_dt, double a_time,
        const std::string &a_file_extension = s_default_file_extension,
        int a_filename_steps_width = s_default_filename_steps_width);

    //! returns m_data_epsilon
    double get_data_epsilon() const;

    //! returns the default data_epsilon
    static double get_default_data_epsilon();

    //! returns m_coords_epsilon
    double get_coords_epsilon() const;

    //! returns the default coords epsilon
    static double get_default_coords_epsilon();
};

#endif /* SMALLDATAIO_HPP_ */
