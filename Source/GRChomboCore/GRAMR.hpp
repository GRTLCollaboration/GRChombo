/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRAMR_HPP_
#define GRAMR_HPP_

// Chombo includes
#include "AMR.H"
#include "Interval.H"

// Other includes
#include "Lagrange.hpp"
#include "VariableType.hpp"
#include <algorithm>
#include <chrono>
#include <ratio>
#include <vector>

// Chombo namespace
#include "UsingNamespace.H"

/// A child of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy (such as the AMRInterpolator)
/**
 *It is necessary for many experimental features and allows us to
 *add said features later without breaking any user code.
 */

// Forward declaration for get_gramrlevels function declarations
class GRAMRLevel;

// Forward declaration for AMRInterpolator
template <typename InterpAlgo> class AMRInterpolator;

class GRAMR : public AMR
{
  public:
    using Clock = std::chrono::steady_clock;
    using Hours = std::chrono::duration<double, std::ratio<3600, 1>>;
    const std::chrono::time_point<Clock> start_walltime = Clock::now();

#if defined(USE_SLURM_INTEGRATION) && defined(CH_NAMESPACE)
  private:
    bool m_in_slurm_job;
    std::chrono::time_point<Clock> end_walltime;

    void set_end_walltime();

  public:
    bool in_slurm_job();

    // defined here due to auto return types
    auto get_remaining_duration()
    {
        if (!m_in_slurm_job)
        {
            return Clock::duration::max();
        }

        return end_walltime - Clock::now();
    }
#endif

  public:
    AMRInterpolator<Lagrange<4>> *m_interpolator; //!< The interpolator pointer

    GRAMR();

    // defined here due to auto return type
    auto get_walltime()
    {
        auto now = Clock::now();
        auto duration = std::chrono::duration_cast<Hours>(now - start_walltime);

        return duration.count();
    }

    // Called after AMR object set up
    virtual void set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator);

    // returs a std::vector of GRAMRLevel pointers
    // similar to AMR::getAMRLevels()
    std::vector<GRAMRLevel *> get_gramrlevels();

    // const version of above
    std::vector<const GRAMRLevel *> get_gramrlevels() const;

    // Fill ghosts on multiple levels
    void fill_multilevel_ghosts(
        const VariableType a_var_type,
        const Interval &a_comps = Interval(0, std::numeric_limits<int>::max()),
        const int a_min_level = 0,
        const int a_max_level = std::numeric_limits<int>::max()) const;
};

#endif /* GRAMR_HPP_ */
