/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOXLOOPS_HPP_
#define BOXLOOPS_HPP_

// Chombo includes
#include "FArrayBox.H"
#include "LevelData.H"

// Our includes
#include "BoxPointers.hpp"
#include "ComputePack.hpp"

// Chombo namespace
#include "UsingNamespace.H"

enum
{
    EXCLUDE_GHOST_CELLS,
    SKIP_GHOST_CELLS = EXCLUDE_GHOST_CELLS,
    INCLUDE_GHOST_CELLS,
    FILL_GHOST_CELLS = INCLUDE_GHOST_CELLS
};

namespace BoxLoops
{
/// Perform the innermost (i.e. the x) loop, standard version with simd
template <typename... compute_ts>
ALWAYS_INLINE void
innermost_loop(const ComputePack<compute_ts...> &compute_pack,
               const BoxPointers &box_pointers, const int iy, const int iz,
               const int loop_lo_x, const int loop_hi_x);

/// Perform the innermost (i.e. the x) loop, simd switched off
template <typename... compute_ts>
ALWAYS_INLINE void
innermost_loop(const ComputePack<compute_ts...> &compute_pack,
               const BoxPointers &box_pointers, const int iy, const int iz,
               const int loop_lo_x, const int loop_hi_x, disable_simd);

/// Performs loop insde the box loop_box and calls compute(...) for all compute
/// classes in the compute_pack  with input FArrayBox 'in' and output FArrayBox
/// 'out'.
template <typename... compute_ts, typename... simd_info>
void loop(const ComputePack<compute_ts...> &compute_pack, const FArrayBox &in,
          FArrayBox &out, const Box &loop_box, simd_info... info);

/// Same as above but for only one compute class (rather than a pack of them)
template <typename compute_t, typename... simd_info>
std::enable_if_t<!is_compute_pack<compute_t>::value, void>
loop(compute_t compute_class, const FArrayBox &in, FArrayBox &out,
     const Box &loop_box, simd_info... info);

/// Performs loop insde the whole box of 'out' and calls compute(...) for all
/// compute classes in the compute_pack  with input FArrayBox 'in' and output
/// FArrayBox 'out'.
template <typename... compute_ts, typename... simd_info>
void loop(const ComputePack<compute_ts...> &compute_pack, const FArrayBox &in,
          FArrayBox &out, simd_info... info); // Uses out.box() as loop_box

/// Same as above but for only one compute class (rather than a pack of them)
template <typename compute_t, typename... simd_info>
std::enable_if_t<!is_compute_pack<compute_t>::value, void>
loop(compute_t compute_class, const FArrayBox &in, FArrayBox &out,
     simd_info... info); // Uses out.box() as loop_box

/// Performs loop over all boxes and inside all boxes of the LevelData 'out' and
/// calls compute(...) for all compute
// classes in the compute_pack with input data taken from 'in' and output
// written to 'out'  MK: Could give the ghost treatment a default argument but I
// think it's better to force the user to make a concious decision  Wrong
// fill_ghosts can give errors that are very hard to debug
template <typename... compute_ts, typename... simd_info>
void loop(const ComputePack<compute_ts...> &compute_pack,
          const LevelData<FArrayBox> &in, LevelData<FArrayBox> &out,
          bool fill_ghosts, simd_info... info);

/// Same as above but for only one compute class (rather than a pack of them)
template <typename compute_t, typename... simd_info>
std::enable_if_t<!is_compute_pack<compute_t>::value, void>
loop(compute_t compute_class, const LevelData<FArrayBox> &in,
     LevelData<FArrayBox> &out, bool fill_ghosts, simd_info... info);
} // namespace BoxLoops

#include "BoxLoops.impl.hpp"

#endif /* BOXLOOPS_HPP_ */
