/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRAMR_HPP_
#define GRAMR_HPP_

// Chombo includes
#include "AMR.H"

// Other includes
#include "AMRInterpolator.hpp"
#include "Lagrange.hpp"
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

class GRAMR : public AMR
{
  private:
    using Clock = std::chrono::steady_clock;
    using Hours = std::chrono::duration<double, std::ratio<3600, 1>>;
    std::chrono::time_point<Clock> start_time = Clock::now();

  public:
    AMRInterpolator<Lagrange<4>> *m_interpolator; //!< The interpolator pointer

    GRAMR();

    // defined here due to auto return type
    auto get_walltime()
    {
        auto now = Clock::now();
        auto duration = std::chrono::duration_cast<Hours>(now - start_time);

        return duration.count();
    }

    // Called after AMR object set up
    virtual void set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator);

    // returs a std::vector of GRAMRLevel pointers
    // similar to AMR::getAMRLevels()
    std::vector<GRAMRLevel *> get_gramrlevels();

    // const version of above
    std::vector<const GRAMRLevel *> get_gramrlevels() const;

    int get_max_level() const { return m_max_level; }
    int get_restart_step() const { return m_restart_step; }
};

#endif /* GRAMR_HPP_ */
