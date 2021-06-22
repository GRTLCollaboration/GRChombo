/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPHEROIDALEXTRACTION_HPP_
#define SPHEROIDALEXTRACTION_HPP_

#include "SpheroidalGeometry.hpp"
#include "SurfaceExtraction.hpp"

//! A child class of SurfaceExtraction for extraction on spheroidal shells
class SpheroidalExtraction : public SurfaceExtraction<SpheroidalGeometry>
{
  public:
    struct params_t : SurfaceExtraction::params_t
    {
        int &num_extraction_radii = num_surfaces;
        std::vector<double> &extraction_radii = surface_param_values;
        int &num_points_t = num_points_u;
        int &num_points_phi = num_points_v;
        std::array<double, CH_SPACEDIM> center; //!< the center of the
                                                //!< spheroidal shells
        std::array<double, CH_SPACEDIM> &extraction_center = center;
        double zaxis_over_xaxis;

        // constructor
        params_t() = default;

        // copy constructor defined due to references pointing to the wrong
        // things with the default copy constructor
        params_t(const params_t &params)
            : center(params.center), zaxis_over_xaxis(params.zaxis_over_xaxis),
              SurfaceExtraction::params_t(params)
        {
        }
    };

    SpheroidalExtraction(const params_t &a_params, double a_dt, double a_time,
                         bool a_first_step, double a_restart_time = 0.0)
        : SurfaceExtraction(
              SpheroidalGeometry(a_params.center, a_params.zaxis_over_xaxis),
              a_params, a_dt, a_time, a_first_step, a_restart_time)
    {
    }

    // Allows a general form to be passed for derivs
    // vars_t = std::tuple<int, VariableType, Derivative>
    SpheroidalExtraction(const params_t &a_params,
                         const std::vector<vars_t> &a_vars, double a_dt,
                         double a_time, bool a_first_step,
                         double a_restart_time = 0.0)
        : SpheroidalExtraction(a_params, a_dt, a_time, a_first_step,
                               a_restart_time)
    {
        add_vars(a_vars);
    }

    // if only variables passed, assume they are evolution ones
    SpheroidalExtraction(const params_t &a_params,
                         const std::vector<int> &a_vars, double a_dt,
                         double a_time, bool a_first_step,
                         double a_restart_time = 0.0)
        : SpheroidalExtraction(a_params, a_dt, a_time, a_first_step,
                               a_restart_time)
    {
        add_evolution_vars(a_vars);
    }

    // Add specific methods?
};

#endif /* SPHEROIDALEXTRACTION_HPP_ */
