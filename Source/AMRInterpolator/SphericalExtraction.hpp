/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPHERICALEXTRACTION_HPP_
#define SPHERICALEXTRACTION_HPP_

#include "SphericalGeometry.hpp"
#include "SphericalHarmonics.hpp"
#include "SurfaceExtraction.hpp"

//! A child class of SurfaceExtraction for extraction on spherical shells
class SphericalExtraction : public SurfaceExtraction<SphericalGeometry>
{
  public:
    struct params_t : SurfaceExtraction::params_t
    {
        int &num_extraction_radii = num_surfaces;
        std::vector<double> &extraction_radii = surface_param_values;
        int &num_points_theta = num_points_u;
        int &num_points_phi = num_points_v;
        std::array<double, CH_SPACEDIM> center; //!< the center of the spherical
                                                //!< shells
        std::array<double, CH_SPACEDIM> &extraction_center = center;
        int num_modes; //!< the number of modes to extract
        std::vector<std::pair<int, int>> modes; //!< the modes to extract
                                                //!< l = first, m = second
    };
    const std::array<double, CH_SPACEDIM> m_center;
    const int m_num_modes;
    const std::vector<std::pair<int, int>> m_modes;

    SphericalExtraction(const params_t &a_params, double a_dt, double a_time,
                        bool a_first_step, double a_restart_time = 0.0)
        : SurfaceExtraction(a_params.center, a_params,
                            a_dt, a_time, a_first_step, a_restart_time),
          m_center(a_params.center), m_num_modes(a_params.num_modes),
          m_modes(a_params.modes)
    {
    }

    SphericalExtraction(const params_t &a_params,
                        const std::vector<std::pair<int, Derivative>> &a_vars,
                        double a_dt, double a_time, bool a_first_step,
                        double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time)
    {
        add_vars(a_vars);
    }

    SphericalExtraction(const params_t &a_params,
                        const std::vector<int> &a_vars, double a_dt,
                        double a_time, bool a_first_step,
                        double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time)
    {
        add_vars(a_vars);
    }

    //! Computes the amplitude of a mode in the spin weighted spherical harmonic
    //! decomposition of a complex valued function, a_function
    //! (normalised by 1/r^2), over each spherical shell
    std::pair<std::vector<double>, std::vector<double>> integrate_mode(
        int es, int el, int em,
        std::function<std::pair<double, double>(std::vector<double>, double,
                                                double, double)>
            a_function)
    {
        auto integrand_re = [center = m_center, &geom = m_geom, es, el,
                             em, &a_function](std::vector<double> a_data_here,
                                              double r, double theta,
                                              double phi) {
            // note that spin_Y_lm requires the coordinates with the center
            // at the origin
            double x = geom.get_grid_coord(0, r, theta, phi) - center[0];
            double y = geom.get_grid_coord(1, r, theta, phi) - center[1];
            double z = geom.get_grid_coord(2, r, theta, phi) - center[2];
            SphericalHarmonics::Y_lm_t<double> Y_lm =
                SphericalHarmonics::spin_Y_lm(x, y, z, es, el, em);
            auto function_here = a_function(a_data_here, r, theta, phi);
            return (function_here.first * Y_lm.Real +
                    function_here.second * Y_lm.Im) /
                   (r * r);
        };

        auto integrand_im = [center = m_center, &geom = m_geom, es, el,
                             em, &a_function](std::vector<double> a_data_here,
                                              double r, double theta,
                                              double phi) {
            // note that spin_Y_lm requires the coordinates with the center
            // at the origin
            double x = geom.get_grid_coord(0, r, theta, phi) - center[0];
            double y = geom.get_grid_coord(1, r, theta, phi) - center[1];
            double z = geom.get_grid_coord(2, r, theta, phi) - center[2];
            SphericalHarmonics::Y_lm_t<double> Y_lm =
                SphericalHarmonics::spin_Y_lm(x, y, z, es, el, em);
            auto function_here = a_function(a_data_here, r, theta, phi);
            return (function_here.second * Y_lm.Real -
                    function_here.first * Y_lm.Im) /
                   (r * r);
        };

        auto integral_re = integrate(integrand_re, IntegrationMethod::simpson);
        auto integral_im = integrate(integrand_im, IntegrationMethod::simpson);

        return std::make_pair(integral_re, integral_im);
    }
};

#endif /* SPHERICALEXTRACTION_HPP_ */
