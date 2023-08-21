/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SmallDataIOReader.hpp"
#include "AlwaysAssert.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>

// Chombo headers
#include "SPMD.H"

// Chombo namespace
#include "UsingNamespace.H"

void SmallDataIOReader::file_structure_t::clear()
{
    block_starts.clear();
    num_data_rows.clear();
    num_header_rows.clear();
    num_coords_columns.clear();
    coords_width = 0;
    num_data_columns.clear();
    data_width = 0;
    num_columns.clear();
}

// Constructor
SmallDataIOReader::SmallDataIOReader() : m_structure_defined(false) {}

// Destructor
SmallDataIOReader::~SmallDataIOReader()
{
    if (m_file.is_open())
    {
        close();
    }
}

void SmallDataIOReader::abort_if_not_rank_zero() const
{
    always_assert(procID() == 0);
}

// Opens the file and sets m_filename. Note that this does not determine the
// file structure
void SmallDataIOReader::open(const std::string &a_filename)
{
    m_filename = a_filename;
    m_structure_defined = false;
    if (procID() == 0)
    {
        m_file.open(m_filename);

        // check file opening successful
        if (!m_file)
        {
            std::string error_message =
                "SmallDataIOReader::open: unable to open file ";
            error_message += m_filename;
            MayDay::Error(error_message.c_str(), 1);
        }
    }
}

// Closes the file
void SmallDataIOReader::close()
{
    if (procID() == 0 && m_file.is_open())
    {
        m_file.close();
    }
    m_filename.clear();
    m_structure_defined = false;
    m_file_structure.clear();
}

// Parses the file and determines its structure
void SmallDataIOReader::determine_file_structure()
{
    // can't broadcast bool
    int structure_defined_int;
    if (procID() == 0)
    {
        // first check file is open
        always_assert(m_file.is_open());

        // move file stream position to start of file
        m_file.clear();
        m_file.seekg(0, std::ios::beg);

        // go through each line and determine structure
        std::string line;
        int consecutive_empty_line_count = 2;
        int block_counter = 0; // assume we always have one block
        int current_position = m_file.tellg();
        int block_start_position = current_position;
        int header_row_counter = 0;
        int data_row_counter = 0;
        while (std::getline(m_file, line))
        {
            if (!line.empty())
            {
                if (consecutive_empty_line_count == 2)
                {
                    // start of new block
                    block_start_position = current_position;
                }
                consecutive_empty_line_count = 0;

                // header rows start with '#'
                if (line.find("#") != std::string::npos)
                    ++header_row_counter;
                else
                {
                    if (data_row_counter++ == 0)
                    {
                        // only count a new block if it contains a data row
                        m_file_structure.block_starts.push_back(
                            block_start_position);
                        ++block_counter;
                        // determine column structure from first data row in
                        // block get a vector of the widths of the columns
                        // including preceeding whitespace
                        std::vector<int> widths;
                        std::string::size_type start_whitespace = 0;
                        while (!(start_whitespace == std::string::npos))
                        {
                            std::string::size_type start_non_whitespace =
                                line.find_first_not_of(' ', start_whitespace);
                            std::string::size_type next_start_whitespace =
                                line.find_first_of(' ', start_non_whitespace);
                            int width;
                            if (next_start_whitespace == std::string::npos)
                            {
                                width = line.length() - start_whitespace;
                            }
                            else
                            {
                                width =
                                    next_start_whitespace - start_whitespace;
                            }
                            widths.push_back(width);
                            start_whitespace = next_start_whitespace;
                        }

                        if (block_counter == 1)
                        {
                            // first data row in file so get coord and data
                            // width from this. assume min width is coord
                            // width and max is data width
                            auto widths_minmax_it = std::minmax_element(
                                widths.begin(), widths.end());
                            m_file_structure.coords_width =
                                *(widths_minmax_it.first);
                            m_file_structure.data_width =
                                *(widths_minmax_it.second);
                        }
                        int num_coords_columns =
                            std::count(widths.begin(), widths.end(),
                                       m_file_structure.coords_width);
                        int num_data_columns =
                            std::count(widths.begin(), widths.end(),
                                       m_file_structure.data_width);
                        m_file_structure.num_coords_columns.push_back(
                            num_coords_columns);
                        m_file_structure.num_data_columns.push_back(
                            num_data_columns);
                        m_file_structure.num_columns.push_back(widths.size());
                    }
                }
            }
            else
            {
                if (consecutive_empty_line_count++ == 0 &&
                    (header_row_counter > 0 || data_row_counter > 0))
                {
                    // end of previous block
                    m_file_structure.num_header_rows.push_back(
                        header_row_counter);
                    m_file_structure.num_data_rows.push_back(data_row_counter);
                    header_row_counter = 0;
                    data_row_counter = 0;
                }
            }
            current_position = m_file.tellg();
        }
        // Just in case the file ends without a line break:
        if (data_row_counter > 0)
        {
            m_file_structure.num_header_rows.push_back(header_row_counter);
            m_file_structure.num_data_rows.push_back(data_row_counter);
            header_row_counter = 0;
            data_row_counter = 0;
        }

        m_file_structure.num_blocks = block_counter;

        always_assert(m_file_structure.num_data_rows.size() ==
                          m_file_structure.num_blocks &&
                      m_file_structure.num_header_rows.size() ==
                          m_file_structure.num_blocks &&
                      m_file_structure.num_coords_columns.size() ==
                          m_file_structure.num_blocks &&
                      m_file_structure.num_data_columns.size() ==
                          m_file_structure.num_blocks &&
                      m_file_structure.num_columns.size() ==
                          m_file_structure.num_blocks &&
                      m_file_structure.block_starts.size() ==
                          m_file_structure.num_blocks);
    }
    m_structure_defined = true;
}

// Set structure if known already (e.g. same as another file already
// determined)
void SmallDataIOReader::set_file_structure(
    const SmallDataIOReader::file_structure_t &a_file_structure)
{
    m_file_structure = a_file_structure;
    m_structure_defined = true;
}

// File structure getter
const SmallDataIOReader::file_structure_t &
SmallDataIOReader::get_file_structure(bool a_broadcast_to_all_ranks)
{
    if (a_broadcast_to_all_ranks)
    {
        // First convert std::vectors to Chombo Vectors to use Chombo broadcast
        // functions
        Vector<unsigned long long> block_starts_Vect =
            m_file_structure.block_starts;
        Vector<int> num_data_rows_Vect = m_file_structure.num_data_rows;
        Vector<int> num_header_rows_Vect = m_file_structure.num_header_rows;
        Vector<int> num_coords_columns_Vect =
            m_file_structure.num_coords_columns;
        Vector<int> num_data_columns_Vect = m_file_structure.num_data_columns;
        Vector<int> num_columns_Vect = m_file_structure.num_columns;

        broadcast(m_file_structure.num_blocks, 0);
        broadcast(block_starts_Vect, 0);
        broadcast(num_data_rows_Vect, 0);
        broadcast(num_header_rows_Vect, 0);
        broadcast(num_coords_columns_Vect, 0);
        broadcast(m_file_structure.coords_width, 0);
        broadcast(num_data_columns_Vect, 0);
        broadcast(m_file_structure.data_width, 0);
        broadcast(num_columns_Vect, 0);

        m_file_structure.block_starts = block_starts_Vect.stdVector();
        m_file_structure.num_data_rows = num_data_rows_Vect.stdVector();
        m_file_structure.num_header_rows = num_header_rows_Vect.stdVector();
        m_file_structure.num_coords_columns =
            num_coords_columns_Vect.stdVector();
        m_file_structure.num_data_columns = num_data_columns_Vect.stdVector();
        m_file_structure.num_columns = num_columns_Vect.stdVector();
    }
    else
    {
        always_assert(procID() == 0);
    }

    return m_file_structure;
}

bool SmallDataIOReader::contains_data()
{
    if (!m_structure_defined)
    {
        determine_file_structure();
    }
    bool contains_data = false;

    if (procID() == 0)
    {
        for (int iblock = 0; iblock < m_file_structure.num_blocks; ++iblock)
        {
            contains_data |= (m_file_structure.num_data_rows[iblock] > 0);
        }
    }
#ifdef CH_MPI
    int contains_data_int = static_cast<int>(contains_data);
    MPI_Bcast(&contains_data_int, 1, MPI_INT, 0, Chombo_MPI::comm);
    if (procID() != 0)
    {
        contains_data = static_cast<bool>(contains_data_int);
    }
#endif /* CH_MPI */
    return contains_data;
}

// Get an interval of columns (inclusive) from a block
std::vector<SmallDataIOReader::column_t>
SmallDataIOReader::get_columns(int a_min_column, int a_max_column, int a_block,
                               bool a_broadcast_to_all_ranks)
{
    const int num_columns = a_max_column - a_min_column + 1;
    std::vector<column_t> out(num_columns);
    if (!m_structure_defined)
    {
        determine_file_structure();
    }
    if (procID() == 0)
    {
        always_assert(m_file.is_open());
        always_assert((0 <= a_min_column) &&
                      (a_max_column < m_file_structure.num_columns[a_block]) &&
                      (a_min_column <= a_max_column));
        for (auto &column : out)
        {
            column.resize(m_file_structure.num_data_rows[a_block]);
        }

        // how many characters across is the start of each column
        std::vector<int> start_position(num_columns);
        std::vector<int> column_widths(num_columns);
        for (int ifcolumn = a_min_column; ifcolumn <= a_max_column; ++ifcolumn)
        {
            // ifcolumn is column index in file and icolumn is column index in
            // output
            int icolumn = ifcolumn - a_min_column;
            if (ifcolumn < m_file_structure.num_coords_columns[a_block])
            {
                start_position[icolumn] =
                    ifcolumn * m_file_structure.coords_width;
                column_widths[icolumn] = m_file_structure.coords_width;
            }
            else
            {
                start_position[icolumn] =
                    m_file_structure.num_coords_columns[a_block] *
                        m_file_structure.coords_width +
                    (ifcolumn - m_file_structure.num_coords_columns[a_block]) *
                        m_file_structure.data_width;
                column_widths[icolumn] = m_file_structure.data_width;
            }
        }

        // move stream position to start of block
        m_file.clear();
        m_file.seekg(m_file_structure.block_starts[a_block], std::ios::beg);
        std::string line;

        // assume header rows are all at the top of the block so skip these
        for (int irow = 0; irow < m_file_structure.num_header_rows[a_block];
             ++irow)
        {
            m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        for (int irow = 0; irow < m_file_structure.num_data_rows[a_block];
             ++irow)
        {
            std::getline(m_file, line);
            for (int icolumn = 0; icolumn < num_columns; ++icolumn)
            {
                out[icolumn][irow] = std::stod(line.substr(
                    start_position[icolumn], column_widths[icolumn]));
            }
        }
    }
    if (a_broadcast_to_all_ranks)
    {
        for (auto &column : out)
        {
            Vector<double> out_Vect = column;
            broadcast(out_Vect, 0);
            if (procID() != 0)
            {
                column = out_Vect.stdVector();
            }
        }
    }
    else
    {
        always_assert(procID() == 0);
    }

    return out;
}

// Get all data columns from a block
std::vector<SmallDataIOReader::column_t>
SmallDataIOReader::get_all_data_columns(int a_block,
                                        bool a_broadcast_to_all_ranks)
{
    if (!m_structure_defined)
    {
        determine_file_structure();
    }
    int min_data_column = 0;
    int max_data_column = 0;
    if (procID() == 0)
    {
        min_data_column = m_file_structure.num_coords_columns[a_block];
        max_data_column =
            min_data_column + m_file_structure.num_data_columns[a_block] - 1;
    }
    if (a_broadcast_to_all_ranks)
    {
        Vector<int> min_max_data_columns(2);
        min_max_data_columns[0] = min_data_column;
        min_max_data_columns[1] = max_data_column;
        broadcast(min_max_data_columns, 0);
        if (procID() != 0)
        {
            min_data_column = min_max_data_columns[0];
            max_data_column = min_max_data_columns[1];
        }
    }
    return get_columns(min_data_column, max_data_column, a_block,
                       a_broadcast_to_all_ranks);
}

// Get all columns from a block
std::vector<SmallDataIOReader::column_t>
SmallDataIOReader::get_all_columns(int a_block, bool a_broadcast_to_all_ranks)
{
    if (!m_structure_defined)
    {
        determine_file_structure();
    }
    int min_column = 0;
    int max_column = 0;
    if (procID() == 0)
    {
        max_column = m_file_structure.num_columns[a_block] - 1;
    }
    if (a_broadcast_to_all_ranks)
    {
        broadcast(max_column, 0);
    }
    return get_columns(min_column, max_column, a_block,
                       a_broadcast_to_all_ranks);
}

SmallDataIOReader::column_t
SmallDataIOReader::get_column(int a_column, int a_block,
                              bool a_broadcast_to_all_ranks)
{
    auto out_vect =
        get_columns(a_column, a_column, a_block, a_broadcast_to_all_ranks);
    return out_vect[0];
}

// Returns a vector of numeric values from a header row
std::vector<double>
SmallDataIOReader::get_data_from_header(int a_header_row_number, int a_block,
                                        bool a_broadcast_to_all_ranks)
{
    std::vector<double> out;
    if (!m_structure_defined)
    {
        determine_file_structure();
    }
    if (procID() == 0)
    {
        always_assert(m_file.is_open());
        always_assert(a_header_row_number <
                      m_file_structure.num_header_rows[a_block]);

        // move stream to start of block
        m_file.clear();
        m_file.seekg(m_file_structure.block_starts[a_block], std::ios::beg);
        std::string line;

        // assume header rows are at start of block
        for (int irow = 0; irow < a_header_row_number; ++irow)
        {
            // skip lines before desired row
            m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }

        // get desired header line
        std::getline(m_file, line);

        // find numbers in header using regex
        // I think this takes a long time to compile...
        std::regex number("[+-]?([0-9]*\\.)?[0-9]+");
        auto numbers_begin =
            std::sregex_iterator(line.begin(), line.end(), number);
        auto numbers_end = std::sregex_iterator();
        for (std::sregex_iterator rit = numbers_begin; rit != numbers_end;
             ++rit)
        {
            // put matches in vector
            out.push_back(std::stod((*rit).str()));
        }
    }

    if (a_broadcast_to_all_ranks)
    {
        Vector<double> out_Vect = out;
        broadcast(out_Vect, 0);
        if (procID() != 0)
        {
            out = out_Vect.stdVector();
        }
    }
    else
    {
        always_assert(procID() == 0);
    }

    return out;
}

std::vector<std::string>
SmallDataIOReader::get_header_strings(int a_header_row_number, int a_block,
                                      bool a_broadcast_to_all_ranks)
{
    std::vector<std::string> out(m_file_structure.num_data_columns[a_block]);
    if (!m_structure_defined)
    {
        determine_file_structure();
    }
    if (procID() == 0)
    {
        always_assert(m_file.is_open());
        always_assert(a_header_row_number <
                      m_file_structure.num_header_rows[a_block]);

        // move stream to start of block
        m_file.clear();
        m_file.seekg(m_file_structure.block_starts[a_block], std::ios::beg);

        // assume header rows are at start of block
        for (int irow = 0; irow < a_header_row_number; ++irow)
        {
            // skip lines before desired row
            m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }

        std::string line;
        std::getline(m_file, line);

        for (int icol = 0; icol < out.size(); ++icol)
        {
            int start_idx = m_file_structure.num_coords_columns[a_block] *
                                m_file_structure.coords_width +
                            icol * m_file_structure.data_width;
            if (start_idx < line.size())
            {
                std::string column_header =
                    line.substr(start_idx, m_file_structure.data_width);
                int first_non_whitespace_char = std::distance(
                    column_header.begin(),
                    std::find_if(column_header.begin(), column_header.end(),
                                 [](char c) { return (c != ' '); }));
                out[icol] = column_header.substr(first_non_whitespace_char);
            }
            else
            {
                out[icol] = "";
            }
        }
    }
    if (a_broadcast_to_all_ranks)
    {
        MayDay::Error("SmallDataIOReader::get_header_strings: broadcast to all "
                      "ranks not implemented");
    }

    return out;
}
