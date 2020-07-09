/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"
#include "computeNorm.H"
#include "computeSum.H"
#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>

// This is used by computeSum, computeNorm, etc.
Vector<LevelData<FArrayBox> *> GRAMR::getLevelDataPtrs()
{
    // First get a std::vector of AMRLevel pointers
    const std::vector<AMRLevel *> level_ptrs{getAMRLevels().stdVector()};

    // Instatiate a Vector of LevelData<FArrayBox> pointers as this is the
    // the format that Chombo's computeSum requires
    Vector<LevelData<FArrayBox> *> level_data_ptrs(level_ptrs.size());

    // Now get the level data pointers
    std::transform(
        level_ptrs.begin(), level_ptrs.end(),
        level_data_ptrs.stdVector().begin(), [](AMRLevel *level_ptr) {
            return const_cast<GRLevelData *>(
                &(GRAMRLevel::gr_cast(level_ptr)->getDiagnosticsLevelData()));
            //                           &(GRAMRLevel::gr_cast(level_ptr)->getLevelData()));
        });

    return level_data_ptrs;
}
// Returns the volume-weighted sum of a grid variable
Real GRAMR::compute_sum(const int a_comp, const Real a_dx_coarse)
{
    CH_TIME("GRAMR::compute_sum");
    const Vector<LevelData<FArrayBox> *> level_data_ptrs{getLevelDataPtrs()};
    return computeSum(level_data_ptrs, m_ref_ratios, a_dx_coarse,
                      Interval(a_comp, a_comp), 0);
}

// Returns the volume-weighted p-norm of an interval of grid variables
Real GRAMR::compute_norm(const Interval a_comps, const double a_p,
                         const Real a_dx_coarse)
{
    CH_TIME("GRAMR::compute_norm");
    const Vector<LevelData<FArrayBox> *> level_data_ptrs{getLevelDataPtrs()};
    return computeNorm(level_data_ptrs, m_ref_ratios, a_dx_coarse, a_comps, a_p,
                       0);
}

// Returns the max value of an interval of grid variables
Real GRAMR::compute_max(const Interval a_comps)
{
    CH_TIME("GRAMR::compute_max");
    const Vector<LevelData<FArrayBox> *> level_data_ptrs{getLevelDataPtrs()};
    return computeMax(level_data_ptrs, m_ref_ratios, a_comps, 0);
}

// Returns the min value of an interval of grid variables
Real GRAMR::compute_min(const Interval a_comps)
{
    CH_TIME("GRAMR::compute_min");
    const Vector<LevelData<FArrayBox> *> level_data_ptrs{getLevelDataPtrs()};
    return computeMin(level_data_ptrs, m_ref_ratios, a_comps, 0);
}

// Returns the Infinity norm of an interval of grid variables
// This function is a bit pointless because a_p = 0 in compute_norm does the
// same thing
Real GRAMR::compute_inf_norm(const Interval a_comps)
{
    CH_TIME("GRAMR::compute_inf_norm");
    return std::max(std::abs(compute_max(a_comps)),
                    std::abs(compute_min(a_comps)));
}