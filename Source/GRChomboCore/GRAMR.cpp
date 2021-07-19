/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"

GRAMR::GRAMR() : m_interpolator(nullptr) {}

// Called after AMR object set up
void GRAMR::set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator)
{
    m_interpolator = a_interpolator;
}

// returs a std::vector of GRAMRLevel pointers
// similar to AMR::getAMRLevels()
std::vector<GRAMRLevel *> GRAMR::get_gramrlevels()
{
    std::vector<GRAMRLevel *> out(m_amrlevels.size());
    std::transform(m_amrlevels.stdVector().cbegin(),
                   m_amrlevels.stdVector().cend(), out.begin(),
                   [](AMRLevel *amrlevel_ptr) {
                       return GRAMRLevel::gr_cast(amrlevel_ptr);
                   });
    return out;
}

// const version of above
std::vector<const GRAMRLevel *> GRAMR::get_gramrlevels() const
{
    std::vector<const GRAMRLevel *> out(m_amrlevels.size());
    std::transform(m_amrlevels.constStdVector().cbegin(),
                   m_amrlevels.constStdVector().cend(), out.begin(),
                   [](const AMRLevel *amrlevel_ptr) {
                       return GRAMRLevel::gr_cast(amrlevel_ptr);
                   });

    return out;
}

void GRAMR::fill_multilevel_ghosts(const VariableType a_var_type,
                                   const Interval &a_comps,
                                   const int a_min_level,
                                   const int a_max_level) const
{
    int max_level = std::min(m_finest_level, a_max_level);

    for (int level_idx = a_min_level; level_idx <= max_level; ++level_idx)
    {
        GRAMRLevel &level = *GRAMRLevel::gr_cast(m_amrlevels[level_idx]);
        level.fillAllGhosts(a_var_type, a_comps);
    }
}

bool GRAMR::need_to_regrid(const int a_level) const
{
    bool out = false;
    std::vector<const GRAMRLevel *> gramrlevels = get_gramrlevels();
    for (int level = a_level; level >= 0; --level)
    {
        // this value is irrelevant, just needs to be != 0
        constexpr int fake_steps_left = 1;
        out |= (needToRegrid(level, fake_steps_left) &&
                gramrlevels[a_level]->at_level_timestep_multiple(level));
    }
    return out;
}

void GRAMR::defer_regridding()
{
    for (int level = 0; level <= m_max_level; ++level)
    {
        if (m_steps_since_regrid[level] = 0)
        {
            m_steps_since_regrid[level] += m_regrid_intervals[level];
        }
    }
}
