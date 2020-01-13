/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SURFACEEXTRACTION_HPP_
#define SURFACEEXTRACTION_HPP_

#include "AMRInterpolator.hpp"
#include "CH_assert.H"
#include "DimensionDefinitions.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "SmallDataIO.hpp" // for writing data
#include "SurfaceGeometry.hpp"
#include "UserVariables.hpp"
#include <algorithm>
#include <array>
#include <utility>
#include <vector>

//! A class to store and return the weights associated to a Newton-Cotes formula
//! for numerical integration/quadrature which can be closed (i.e. includes the
//! endpoints) or open (does not include the end points)
class IntegrationMethod
{
  private:
    std::vector<double> m_weights;
    int m_num_weights;
    bool m_is_closed;

  public:
    //! Constructor
    IntegrationMethod(const std::vector<double> &a_weights,
                      bool a_is_closed = true)
        : m_weights(a_weights), m_num_weights(a_weights.size()),
          m_is_closed(a_is_closed)
    {
    }

    //! Checks that this integration method is suitable given the number of
    //! points and periodicity
    inline bool is_valid(int a_num_points, bool a_is_periodic) const
    {
        if (m_is_closed && !a_is_periodic)
        {
            return (a_num_points % m_num_weights == 1 % m_num_weights);
        }
        else
        {
            return (a_num_points % m_num_weights == 0);
        }
    }

    //! Returns whether this IntegrationMethod is closed or not
    inline bool is_closed() const { return m_is_closed; }
    
    //! Returns the weight for a point with given index
    inline double weight(int a_index, int a_num_points,
                         bool a_is_periodic) const
    {
        const int weight_index = a_index % m_num_weights;
        const bool endpoint =
            (a_index == 0 || a_index == a_num_points - 1) && !a_is_periodic;
        if (m_is_closed && !endpoint && weight_index == 0)
            return 2.0 * m_weights[weight_index];
        else
            return m_weights[weight_index];
    }

    static const IntegrationMethod trapezium;
    static const IntegrationMethod midpoint;
    static const IntegrationMethod simpson;
    static const IntegrationMethod boole;
};

//! This class extracts grid variables on 2 dimensional surfaces each
//! parameterised by u and v with different surfaces given by level sets of
//! another parameter
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
    const SurfaceGeometry *m_geom_ptr; //!< the pointer to the geometry class
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
    SurfaceExtraction(SurfaceGeometry *a_geom_ptr, params_t &a_params,
                      double a_dt, double a_time, bool a_first_step,
                      double a_restart_time = 0.0);

    //! add a single variable or derivative of variable
    void add_var(const int a_var, const Derivative a_deriv = Derivative::LOCAL);

    //! add a vector of variables/derivatives of variables
    void add_vars(const std::vector<std::pair<int, Derivative>> &a_vars);

    //! Alternative constructor with a predefined vector of variables
    SurfaceExtraction(SurfaceGeometry *a_geom_ptr, params_t &a_params,
                      const std::vector<std::pair<int, Derivative>> &a_vars,
                      double a_dt, double a_time, bool a_first_step,
                      double a_restart_time = 0.0);

    //! Do the extraction
    // defined here as this is a templated function
    template <typename InterpAlgo>
    void extract(AMRInterpolator<InterpAlgo> *a_interpolator)
    {
        if (a_interpolator == nullptr)
        {
            MayDay::Error("SurfaceExtraction: invalid AMRInterpolator pointer");
        }
        // set up the interpolation query
        InterpolationQuery query(m_num_points * m_params.num_surfaces);
        FOR1(idir) { query.setCoords(idir, m_interp_coords[idir].data()); }
        for (int ivar = 0; ivar < m_vars.size(); ++ivar)
        {
            query.addComp(m_vars[ivar].first, m_interp_data[ivar].data(),
                          m_vars[ivar].second);
        }

        // submit the query
        a_interpolator->interp(query);
        m_done_extraction = true;
    }

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

#endif /* SURFACEEXTRACTION_HPP_ */
