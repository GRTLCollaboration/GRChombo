/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRAMR_HPP_
#define GRAMR_HPP_

#include "AMR.H"
#include "AMRInterpolator.hpp"
#include "Lagrange.hpp"
#include <chrono>
#include <ratio>

/// A child of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy (such as the AMRInterpolator)
/**
 *It is necessary for many experimental features and allows us to
 *add said features later without breaking any user code.
 */
class GRAMR : public AMR
{
  private:
    using Clock = std::chrono::steady_clock;
    using Hours = std::chrono::duration<double, std::ratio<3600, 1>>;
    std::chrono::time_point<Clock> start_time = Clock::now();

    // This is used by computeSum, computeNorm, etc.
    Vector<LevelData<FArrayBox> *> getLevelDataPtrs();

  public:
    AMRInterpolator<Lagrange<4>> *m_interpolator; //!< The interpolator pointer

    GRAMR() { m_interpolator = nullptr; }

    auto get_walltime()
    {
        auto now = Clock::now();
        auto duration = std::chrono::duration_cast<Hours>(now - start_time);

        return duration.count();
    }

    // Called after AMR object set up
    void set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        m_interpolator = a_interpolator;
    }

    // Returns the volume-weighted sum of a grid variable
    Real compute_sum(const int a_comp, const Real a_dx_coarse);

    // Returns the volume-weighted p-norm of an interval of grid variables
    Real compute_norm(const Interval a_comps, const double a_p,
                      const Real a_dx_coarse);

    // Returns the max value of an interval of grid variables
    Real compute_max(const Interval a_comps);

    // Returns the min value of an interval of grid variables
    Real compute_min(const Interval a_comps);

    // Returns the Infinity norm of an interval of grid variables
    // This function is a bit pointless because a_p = 0 in compute_norm does the
    // same
    Real compute_inf_norm(const Interval a_comps);
};

#endif /* GRAMR_HPP_ */
