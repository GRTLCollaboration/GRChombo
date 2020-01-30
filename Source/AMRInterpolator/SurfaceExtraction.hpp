/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SURFACEEXTRACTION_HPP_
#define SURFACEEXTRACTION_HPP_

#include "AMRInterpolator.hpp"
#include "CH_assert.H"
#include "DimensionDefinitions.hpp"
#include "IntegrationMethod.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "SmallDataIO.hpp" // for writing data
#include "UserVariables.hpp"
#include <algorithm>
#include <array>
#include <utility>
#include <vector>

//! This class extracts grid variables on 2 dimensional surfaces each
//! parameterised by u and v with different surfaces given by level sets of
//! another parameter
template <class SurfaceGeometry>
class SurfaceExtraction
{
  public:
    struct params_t
    {
        int num_surfaces; //!< number of surfaces over which to extraction
        std::vector<double>
            surface_param_values; //!< the values of the
                                  //!< parameter that gives the required
                                  //!< surfaces with SurfaceGeom geometry (e.g.
                                  //!< radii for spherical shells)
        int num_points_u; //!< the number of points for the first parameter
                          //!< that parameterises each surface
        int num_points_v; //!< the number of points for the second parameter
                          //!< that parameterises each surfaces
        std::vector<int> extraction_levels; //!< the level on which to do the
                                            //!< extraction for each surface
        bool write_extraction; //!< whether or not to write the extracted data

        int min_extraction_level()
        {
            return *(std::min_element(extraction_levels.begin(),
                                      extraction_levels.end()));
        }
    };

  protected:
    const SurfaceGeometry m_geom; //!< the geometry class which knows about
                                  //!< the particular surface
    const params_t m_params;
    std::vector<std::pair<int, Derivative>> m_vars; //!< the vector of pairs of
    //!< variables and derivatives to extract
    const double m_dt;
    const double m_time;
    const bool m_first_step;
    const double m_restart_time;
    const int m_num_points; //!< the total number of points per surface
    const double m_du;      //!< the grid spacing in u (used in integrate)
    const double m_dv;      //!< the grid spacing in v (used in integrate)

    std::vector<std::vector<double>> m_interp_data;
    std::array<std::vector<double>, CH_SPACEDIM> m_interp_coords;

    bool m_done_extraction; //!< whether or not the extract function has been
                            //!< called for this object

    //! returns the flattened index for m_interp_data and m_interp_coords
    //! associated to given surface, u and v indices
    int index(int a_isurface, int a_iu, int a_iv) const
    {
        return a_isurface * m_num_points + a_iu * m_params.num_points_v + a_iv;
    }

  public:
    //! Normal constructor which requires vars to be added after construction
    //! using add_var or add_vars
    SurfaceExtraction(const SurfaceGeometry &a_geom,
                      const params_t &a_params, double a_dt, double a_time,
                      bool a_first_step, double a_restart_time = 0.0);

    //! add a single variable or derivative of variable
    void add_var(int a_var, const Derivative &a_deriv = Derivative::LOCAL);

    //! add a vector of variables/derivatives of variables
    void add_vars(const std::vector<std::pair<int, Derivative>> &a_vars);

    //! add a vector of variables (no derivatives)
    void add_vars(const std::vector<int> &a_vars);

    //! Alternative constructor with a predefined vector of variables and
    //! derivatives
    SurfaceExtraction(const SurfaceGeometry &a_geom,
                      const params_t &a_params,
                      const std::vector<std::pair<int, Derivative>> &a_vars,
                      double a_dt, double a_time, bool a_first_step,
                      double a_restart_time = 0.0);

    //! Another alternative constructor with a predefined vector of variables
    //! no derivatives
    SurfaceExtraction(const SurfaceGeometry &a_geom,
                      const params_t &a_params, const std::vector<int> &a_vars,
                      double a_dt, double a_time, bool a_first_step,
                      double a_restart_time = 0.0);


//! Do the extraction
template <typename InterpAlgo>
void extract(AMRInterpolator<InterpAlgo> *a_interpolator);

    //! Integrate some integrand dependent on the interpolated data over the
    //! surface. The integrand function should be of the signature
    //! double integrand(std::vector<double> data_here,
    //!     double a_surface_param_value, double a_u, double a_v)
    //! where data_here is a vector of all the interpolated variables at the
    //! point specified by the other arguments.
    std::vector<double> integrate(
        std::function<double(std::vector<double>, double, double, double)>
            a_integrand,
        const IntegrationMethod &a_method_u = IntegrationMethod::trapezium,
        const IntegrationMethod &a_method_v = IntegrationMethod::trapezium);

    //! Write the interpolated data to a file with a block for each surface
    void write_extraction(std::string a_file_prefix) const;

    //! write some integrals to a file at this timestep
    void write_integrals(const std::string &a_filename,
                         const std::vector<std::vector<double>> &a_integrals,
                         const std::vector<std::string> &a_labels = {}) const;

    //! convenience caller for write_integrals in the case of just integral per
    //! surface
    void write_integral(const std::string &a_filename,
                        const std::vector<double> a_integrals,
                        const std::string a_label = "") const;
};

#include "SurfaceExtraction.impl.hpp"

#endif /* SURFACEEXTRACTION_HPP_ */
