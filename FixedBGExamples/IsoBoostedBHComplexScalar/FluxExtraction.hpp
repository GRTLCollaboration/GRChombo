/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FLUXEXTRACTION_HPP_
#define FLUXEXTRACTION_HPP_

#include "SphericalExtraction.hpp"
//!  The class allows extraction of the values of the Flux components on
//!  Spherical shells at specified radii, and integration over those shells
/*!
   The class allows the user to extract data from the grid for the Flux
   components over Spherical shells at specified radii. The values may then be
   written to an output file, or integrated across the surfaces.
*/
class FluxExtraction : public SphericalExtraction
{
  public:
    //! The constructor
    FluxExtraction(SphericalExtraction::params_t &a_params, double a_dt,
                   double a_time, bool a_first_step,
                   double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time)
    {
        add_var(c_Edot, VariableType::diagnostic);
        add_var(c_Jdot, VariableType::diagnostic);
        add_var(c_Mdot, VariableType::diagnostic);
    }

    //! The old constructor which assumes it is called in specificPostTimeStep
    //! so the first time step is when m_time == m_dt
    FluxExtraction(SphericalExtraction::params_t a_params, double a_dt,
                   double a_time, double a_restart_time = 0.0)
        : FluxExtraction(a_params, a_dt, a_time, (a_dt == a_time),
                         a_restart_time)
    {
    }

    // the references of the vars as used in the integrator
    enum M_VARS
    {
        m_Edot,
        m_Jdot,
        m_Mdot,
        NUM_COMPS
    };

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        // extract the values of the Flux scalars on the spheres
        extract(a_interpolator);

        // this would write out the values at every point on the sphere
        if (m_params.write_extraction)
        {
            write_extraction("FluxExtractionOut_");
        }

        // Setup to integrate Stress and Edot
        std::vector<std::vector<double>> flux_integrals(NUM_COMPS);
        add_var_integrand(m_Edot, flux_integrals[m_Edot],
                          IntegrationMethod::simpson);
        add_var_integrand(m_Jdot, flux_integrals[m_Jdot],
                          IntegrationMethod::simpson);
        add_var_integrand(m_Mdot, flux_integrals[m_Mdot],
                          IntegrationMethod::simpson);

        // do the integration over the surface
        integrate();

        // write the integrals
        std::vector<std::string> labels(NUM_COMPS);
        labels[m_Edot] = "Edot";
        labels[m_Jdot] = "Jdot";
        labels[m_Mdot] = "Mdot";
        std::string filename = "SurfaceIntegrals";
        write_integrals(filename, flux_integrals, labels);
    }
};

#endif /* FLUXEXTRACTION_HPP_ */
