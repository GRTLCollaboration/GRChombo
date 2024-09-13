/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CosmoAMR_HPP_
#define CosmoAMR_HPP_

#include "GRAMR.hpp"
#include "PunctureTracker.hpp"

/// A child of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy, and those of GRAMR
/**
 * This object inherits from GRAMR and adds tools required for BH spacetimes
 */
class CosmoAMR : public GRAMR
{
  private:
    double m_K_mean;
    double m_rho_mean;
    double m_S_mean;

  public:
    PunctureTracker m_puncture_tracker;

    CosmoAMR() {}

    void set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator) override
    {
        GRAMR::set_interpolator(a_interpolator);
        m_puncture_tracker.set_interpolator(a_interpolator);
    }

    void set_K_mean(double a_K_mean) { m_K_mean = a_K_mean; }
    double get_K_mean() { return m_K_mean; }

    void set_rho_mean(double a_rho_mean) { m_rho_mean = a_rho_mean; }
    double get_rho_mean() { return m_rho_mean; }

    void set_S_mean(double a_S_mean) { m_S_mean = a_S_mean; }
    double get_S_mean() { return m_S_mean; }
};

#endif /* CosmoAMR_HPP_ */
