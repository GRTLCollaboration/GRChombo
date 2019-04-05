/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SMALLDATAIO_HPP_)
#error "This file should only be included through SmallDataIO.hpp"
#endif

#ifndef SMALLDATAIO_IMPL_HPP_
#define SMALLDATAIO_IMPL_HPP_

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
    if (m_rank == 0 && m_time == m_dt)
    {
        // all header lines start with a '#'.
        m_file << "#";
        // the coords are only given a width of 12 characters
        for (int istr = 0; istr < a_pre_header_strings.size(); ++istr)
        {
            if (istr == 0)
            {
                m_file << std::setw(11) << a_pre_header_strings[istr];
            }
            else
            {
                m_file << std::setw(12) << a_pre_header_strings[istr];
            }
        }
        // the data is given a width of 20 characters
        for (std::string header_item : a_header_strings)
        {
            m_file << std::setw(20) << header_item;
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
        m_file << std::fixed << std::setprecision(7);
        for (double coord : a_coords)
        {
            m_file << std::setw(12) << coord;
        }
        m_file << std::scientific << std::setprecision(10);
        for (double data : a_data)
        {
            m_file << std::setw(20) << data;
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

void SmallDataIO::remove_duplicate_time_data()
{
    constexpr double epsilon = 1.0e-8;
    if (m_rank == 0 && m_restart_time > 0. && m_mode == APPEND &&
        m_time < m_restart_time + m_dt + epsilon)
    {
        // copy lines with time <= restart_time into a temporary file
        std::string line;
        std::string temp_filename = m_filename + ".temp";
        std::ofstream temp_file(temp_filename);
        while (std::getline(m_file, line))
        {
            if (!(line.find("#") == std::string::npos))
            {
                temp_file << line << "\n";
            }
            else if (std::stod(line.substr(0, 12)) <= m_restart_time + epsilon)
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
        m_file.open(m_filename, std::ios::app);
    }
}

#endif /* SMALLDATAIO_IMPL_HPP_ */
