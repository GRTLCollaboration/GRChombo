/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOXLOOPS_HPP_)
#error "This file should only be included through BoxLoops.hpp"
#endif

#ifndef BOXLOOPS_IMPL_HPP_
#define BOXLOOPS_IMPL_HPP_

#include "BoxPointers.hpp"
#include "Cell.hpp"
#include "DebuggingTools.hpp"
#include "simd.hpp"

template <typename... compute_ts>
void BoxLoops::innermost_loop(const ComputePack<compute_ts...> &compute_pack,
                              const BoxPointers &box_pointers, const int iy,
                              const int iz, const int loop_lo_x,
                              const int loop_hi_x)
{
    int simd_width = simd<double>::simd_len;
    int x_simd_max =
        loop_lo_x +
        simd_width * (((loop_hi_x - loop_lo_x + 1) / simd_width) - 1);
// SIMD LOOP
#ifdef __INTEL_COMPILER
#pragma novector
#elif !defined(__clang__)
#pragma omp simd safelen(1)
#endif /* __INTEL_COMPILER */
    for (int ix = loop_lo_x; ix <= x_simd_max; ix += simd_width)
    {
        compute_pack.call_compute(
            Cell<simd<double>>(IntVect(ix, iy, iz), box_pointers));
    }

    // REMAINDER LOOP
    for (int ix = x_simd_max + simd<double>::simd_len; ix <= loop_hi_x; ++ix)
    {
        compute_pack.call_compute(
            Cell<double>(IntVect(ix, iy, iz), box_pointers));
    }
}

template <typename... compute_ts>
void BoxLoops::innermost_loop(const ComputePack<compute_ts...> &compute_pack,
                              const BoxPointers &box_pointers, const int iy,
                              const int iz, const int loop_lo_x,
                              const int loop_hi_x, disable_simd)
{
    for (int ix = loop_lo_x; ix <= loop_hi_x; ++ix)
    {
        compute_pack.call_compute(
            Cell<double>(IntVect(ix, iy, iz), box_pointers));
    }
}

template <typename... compute_ts, typename... simd_info>
void BoxLoops::loop(const ComputePack<compute_ts...> &compute_pack,
                    const FArrayBox &in, FArrayBox &out, const Box &loop_box,
                    simd_info... info)
{
    // Makes sure we are not requesting data outside the box of 'out'
    CH_assert(out.box().contains(loop_box));

    BoxPointers box_pointers(in, out);

    const int *loop_lo = loop_box.loVect();
    const int *loop_hi = loop_box.hiVect();

#pragma omp parallel for default(shared) collapse(CH_SPACEDIM - 1)
#if CH_SPACEDIM >= 3
    for (int iz = loop_lo[2]; iz <= loop_hi[2]; ++iz)
#endif
        for (int iy = loop_lo[1]; iy <= loop_hi[1]; ++iy)
        {
#ifdef EQUATION_DEBUG_MODE
            innermost_loop(compute_pack, box_pointers, iy, iz, loop_lo[0],
                           loop_hi[0], disable_simd());
#else
        innermost_loop(compute_pack, box_pointers, iy, iz, loop_lo[0],
                       loop_hi[0], std::forward<simd_info>(info)...);
#endif
        }
}

template <typename compute_t, typename... simd_info>
std::enable_if_t<!is_compute_pack<compute_t>::value, void>
BoxLoops::loop(compute_t compute_class, const FArrayBox &in, FArrayBox &out,
               const Box &loop_box, simd_info... info)
{
    loop(make_compute_pack(compute_class), in, out, loop_box,
         std::forward<simd_info>(info)...);
}

template <typename... compute_ts, typename... simd_info>
void BoxLoops::loop(const ComputePack<compute_ts...> &compute_pack,
                    const FArrayBox &in, FArrayBox &out, simd_info... info)
{
    loop(compute_pack, in, out, out.box(), std::forward<simd_info>(info)...);
}

template <typename compute_t, typename... simd_info>
std::enable_if_t<!is_compute_pack<compute_t>::value, void>
BoxLoops::loop(compute_t compute_class, const FArrayBox &in, FArrayBox &out,
               simd_info... info)
{
    loop(make_compute_pack(compute_class), in, out,
         std::forward<simd_info>(info)...);
}

template <typename... compute_ts, typename... simd_info>
void BoxLoops::loop(const ComputePack<compute_ts...> &compute_pack,
                    const LevelData<FArrayBox> &in, LevelData<FArrayBox> &out,
                    bool fill_ghosts, simd_info... info)
{
    DataIterator dit0 = in.dataIterator();
    int nbox = dit0.size();
    for (int ibox = 0; ibox < nbox; ++ibox)
    {
        DataIndex di = dit0[ibox];
        const FArrayBox &in_fab = in[di];
        FArrayBox &out_fab = out[di];

        Box out_box;
        if (fill_ghosts)
            out_box = out_fab.box();
        else
            out_box = in.disjointBoxLayout()[di];

        loop(compute_pack, in_fab, out_fab, out_box,
             std::forward<simd_info>(info)...);
    }
}

template <typename compute_t, typename... simd_info>
std::enable_if_t<!is_compute_pack<compute_t>::value, void>
BoxLoops::loop(compute_t compute_class, const LevelData<FArrayBox> &in,
               LevelData<FArrayBox> &out, bool fill_ghosts, simd_info... info)
{
    // TODO think about perfect forwarding of compute_class
    loop(make_compute_pack(compute_class), in, out, fill_ghosts,
         std::forward<simd_info>(info)...);
}

#endif /* BOXLOOPS_IMPL_HPP_ */
