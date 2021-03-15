/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if defined(USE_SLURM_INTEGRATION) && defined(CH_NAMESPACE)
#include <slurm/slurm.h>
#endif
#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"

GRAMR::GRAMR() : m_interpolator(nullptr)
{
#if defined(USE_SLURM_INTEGRATION) && defined(CH_NAMESPACE)
    set_end_walltime();
#endif
}

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

#if defined(USE_SLURM_INTEGRATION) && defined(CH_NAMESPACE)
void GRAMR::set_end_walltime()
{
    uint32_t jobid;
    int success = slurm_pid2jobid(getpid(), &jobid);
    if (success == 0)
    {
        m_in_slurm_job = true;
        long remaining_secs = slurm_get_rem_time(jobid);
        end_walltime = Clock::now() + std::chrono::seconds(remaining_secs);
    }
    else
    {
        m_in_slurm_job = false;
    }
}

bool GRAMR::in_slurm_job() { return m_in_slurm_job; }
#endif /* USE_SLURM_INTEGRATION */
