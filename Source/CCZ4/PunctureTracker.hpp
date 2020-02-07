/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PUNCTURETRACKER_HPP_
#define PUNCTURETRACKER_HPP_

#include "AMRInterpolator.hpp"
#include "BHAMR.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "SimulationParametersBase.hpp"
#include "SmallDataIO.hpp"   // for writing data
#include "UserVariables.hpp" // Needs c_shift etc

//!  The class tracks the puncture locations by integrating the shift at
//!  The puncture position
class PunctureTracker
{
  private:
    //! Params for puncture tracking
    const int m_num_punctures;
    const double m_time;
    const double m_restart_time;
    const double m_dt;
    const std::string m_punctures_filename;

  public:
    //! The constructor
    PunctureTracker(const double a_time, const double a_restart_time,
                    const double a_dt, const std::string a_checkpoint_prefix,
                    const int a_num_punctures = 2)
        : m_num_punctures(a_num_punctures), m_time(a_time),
          m_punctures_filename(a_checkpoint_prefix + "Punctures"),
          m_restart_time(a_restart_time), m_dt(a_dt)
    {
    }

    //! set puncture locations on restart
    void restart_punctures(BHAMR &a_gramr,
                           std::vector<std::array<double, CH_SPACEDIM>>
                               initial_puncture_coords) const;

    //! set and write initial puncture locations
    void set_initial_punctures(BHAMR &a_gramr,
                               std::vector<std::array<double, CH_SPACEDIM>>
                                   initial_puncture_coords) const;

    //! Set punctures post restart if m_time > 0
    void read_in_punctures(BHAMR &a_gramr) const;

    //! Execute the tracking and write out
    void execute_tracking(BHAMR &a_gramr,
                          const bool write_punctures = true) const;

    //! Use the interpolator to get the value of the shift at
    //! given coords
    void get_interp_shift(
        std::vector<std::array<double, CH_SPACEDIM>> &interp_shift,
        BHAMR &a_gramr,
        std::vector<std::array<double, CH_SPACEDIM>> puncture_coords) const;

    //! Get a vector of the puncture coords - used for write out
    std::vector<double> get_puncture_vector(
        std::vector<std::array<double, CH_SPACEDIM>> puncture_coords) const;
};

#include "PunctureTracker.impl.hpp"

#endif /* PUNCTURETRACKER_HPP_ */
