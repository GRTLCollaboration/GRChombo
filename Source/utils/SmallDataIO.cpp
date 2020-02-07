/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SmallDataIO.hpp"

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
    if (m_rank == 0)
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
    if (m_rank == 0)
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
    if (m_rank == 0)
    {
        m_file << "\n\n";
    }
}

void SmallDataIO::remove_duplicate_time_data(const bool keep_m_time_data)
{
    constexpr double epsilon = 1.0e-8;
    if (m_rank == 0 && m_restart_time > 0. && m_mode == APPEND &&
        m_time < m_restart_time + m_dt + epsilon)
    {
        // copy lines with time < m_time into a temporary file
        m_file.seekg(0);
        std::string line;
        std::string temp_filename = m_filename + ".temp";
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
                     m_time + sign * epsilon)
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
    if (m_rank == 0)
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
