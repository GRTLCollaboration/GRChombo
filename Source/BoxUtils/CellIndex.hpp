/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CELLINDEX_HPP_
#define CELLINDEX_HPP_

#include "AlwaysInline.hpp"

/// Represents an index in the flattened Chombo array
/** CellIndex only stores one integer (the position in the flattened Chombo
 *array of the cell the CellIndex represents). A CellIndex expects a bool as
 *template arugment to check whether it represent "in" or "out" data of the
 *FABDriver. For better readibility  use ::CellIndexIn and ::CellIndexOut
 *instead of CellIndex<false> and CellIndex<true> respectively.
 * @tparam is_out_data A bool that makes sure that the indices in the "in" and
 *"out" arrays are never confused.
 **/
template <bool is_out_data> struct CellIndex
{
    int m_index;

    ALWAYS_INLINE
    CellIndex(int index) : m_index(index) {}

    ALWAYS_INLINE
    operator int() const { return m_index; }
};
using CellIndexIn =
    CellIndex<false>; //!< Short hand for a CellIndex in the "in"-box.
using CellIndexOut =
    CellIndex<true>; //!< Short hand for a CellIndex in the "out"-box.

#endif /* CELLINDEX_HPP_ */
