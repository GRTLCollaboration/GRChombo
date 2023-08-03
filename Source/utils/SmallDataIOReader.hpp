/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SMALLDATAIOREADER_HPP
#define SMALLDATAIOREADER_HPP

#include <fstream>
#include <string>
#include <vector>

// A class to read files written using SmallDataIO.
// This class will only work on rank 0 so make sure to guard with
// if (procID() == 0) { <SmallDataIOReader code> }

class SmallDataIOReader
{
  public:
    using column_t = std::vector<double>;
    // A struct for information about the structure of a SmallDataIO file
    struct file_structure_t
    {
        int num_blocks; // a block is separated by 2 blank lines
        std::vector<unsigned long long int>
            block_starts; // position offsets from the beginning of the file
        std::vector<int> num_data_rows;      // the number of data rows in each
                                             // block
        std::vector<int> num_header_rows;    // the number of header rows in
                                             // each block
        std::vector<int> num_coords_columns; // number of coord columns in
                                             // each block - assume constant
                                             // in each block
        int coords_width;                    // assume the same throughout file
        std::vector<int> num_data_columns;   // number of data columns in each
                                             // block

        int data_width;               // assume the same throughout file
        std::vector<int> num_columns; // sum of coord and data columns in each
                                      // block
        void clear();
    };

  protected:
    std::string m_filename;
    std::ifstream m_file;
    file_structure_t m_file_structure;
    bool m_structure_defined;

    void abort_if_not_rank_zero() const;

  public:
    // Constructor
    SmallDataIOReader();

    // Destructor
    ~SmallDataIOReader();

    // Opens the file and sets m_filename. Note that this does not determine the
    // file structure
    void open(const std::string &a_filename);

    // Closes the file
    void close();

    // Parses the file and determines its structure
    void determine_file_structure();

    // Set structure if known already (e.g. same as another file already
    // determined)
    void set_file_structure(const file_structure_t &a_file_structure);

    // File struture getter
    const file_structure_t &
    get_file_structure(bool a_broadcast_to_all_ranks = true);

    // Returns true if there are some valid non-header rows and false otherwise
    bool contains_data();

    // Get an interval of columns (inclusive) from a block
    std::vector<column_t> get_columns(int a_min_column, int a_max_column,
                                      int a_block = 0,
                                      bool a_broadcast_to_all_ranks = true);

    // Get all data columns from a block
    std::vector<column_t>
    get_all_data_columns(int a_block = 0, bool a_broadcast_to_all_ranks = true);

    // Get all columns from a block
    std::vector<column_t> get_all_columns(int a_block = 0,
                                          bool a_broadcast_to_all_ranks = true);

    // Get a single column from a block
    column_t get_column(int a_column, int a_block = 0,
                        bool a_broadcast_to_all_ranks = true);

    // Returns a vector of numeric values from a header row
    std::vector<double>
    get_data_from_header(int a_header_row_number, int a_block = 0,
                         bool a_broadcast_to_all_ranks = true);

    // Returns a vector of strings from a header row
    std::vector<std::string>
    get_header_strings(int a_header_row_number, int a_block = 0,
                       bool a_broadcast_to_all_ranks = false);
};

#endif /* SMALLDATAIOREADER_HPP */
