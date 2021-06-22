/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PUNCTURETRACKER_HPP_
#define PUNCTURETRACKER_HPP_

#include "AMRInterpolator.hpp"
#include "AlwaysInline.hpp"
#include "Lagrange.hpp"

//!  The class tracks the puncture locations by integrating the shift at
//!  The puncture position
class PunctureTracker
{
  private:
    //! Params for puncture tracking
    int m_num_punctures;
    std::vector<std::array<double, CH_SPACEDIM>> m_puncture_coords;
    std::vector<std::array<double, CH_SPACEDIM>> m_puncture_shift;
    int m_min_level; //!< the min level on which punctures will be
                     //!< (to fill ghosts)

    std::string m_punctures_filename;

    // saved pointer to external interpolator
    AMRInterpolator<Lagrange<4>> *m_interpolator;

  public:
    //! The constructor
    PunctureTracker() : m_num_punctures(0), m_interpolator(nullptr) {}

    //! set puncture locations on start (or restart)
    //! this needs to be done before 'setupAMRObject'
    //! if the puncture locations are required for Tagging Criteria
    void initial_setup(const std::vector<std::array<double, CH_SPACEDIM>>
                           &initial_puncture_coords,
                       const std::string &a_filename = "punctures",
                       const std::string &a_output_path = "",
                       const int a_min_level = 0);

    //! set puncture locations on start (or restart)
    void restart_punctures();

    //! Execute the tracking and write out
    void execute_tracking(double a_time, double a_restart_time, double a_dt,
                          const bool write_punctures = true);

    ALWAYS_INLINE void
    set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        m_interpolator = a_interpolator;
    }

    // function to get punctures
    ALWAYS_INLINE const std::vector<std::array<double, CH_SPACEDIM>> &
    get_puncture_coords() const
    {
        return m_puncture_coords;
    }

  private:
    //! set and write initial puncture locations
    void set_initial_punctures();

    //! Set punctures post restart if m_time > 0
    void read_in_punctures(int a_int_step, double a_restart_time);

    //! Use the interpolator to get the value of the shift at
    //! given coords
    void interp_shift();

    //! Get a vector of the puncture coords - used for write out
    std::vector<double> get_puncture_vector() const;
};

#endif /* PUNCTURETRACKER_HPP_ */
