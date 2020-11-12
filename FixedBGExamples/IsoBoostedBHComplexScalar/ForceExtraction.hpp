/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FORCEEXTRACTION_HPP_
#define FORCEEXTRACTION_HPP_

#include "SpheroidalExtraction.hpp"
//!  The class allows extraction of the values of the force components on
//!  spheroidal shells at specified radii, and integration over those shells
/*!
   The class allows the user to extract data from the grid for the force
   components over spheroidal shells at specified radii. The values may then be
   written to an output file, or integrated across the surfaces.
*/
class ForceExtraction : public SpheroidalExtraction
{
  public:
    //! The constructor
    ForceExtraction(SpheroidalExtraction::params_t &a_params, double a_dt,
                    double a_time, bool a_first_step,
                    double a_restart_time = 0.0)
        : SpheroidalExtraction(a_params, a_dt, a_time, a_first_step,
                               a_restart_time)
    {
        add_var(c_Stress, VariableType::diagnostic);
        add_var(c_BHMom, VariableType::diagnostic);
    }

    //! The old constructor which assumes it is called in specificPostTimeStep
    //! so the first time step is when m_time == m_dt
    ForceExtraction(SpheroidalExtraction::params_t a_params, double a_dt,
                    double a_time, double a_restart_time = 0.0)
        : ForceExtraction(a_params, a_dt, a_time, (a_dt == a_time),
                          a_restart_time)
    {
    }

    // the references of the vars as used in the integrator
    enum M_VARS
    {
        m_Stress,
        m_BHMom
    };

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        // extract the values of the Force scalars on the spheres
        extract(a_interpolator);

        // this would write out the values at every point on the sphere
        if (m_params.write_extraction)
        {
            write_extraction("Force4ExtractionOut_");
        }

        // Setup to integrate Stress and BHMom
        std::vector<std::vector<double>> force_integrals(2);
        add_var_integrand(m_Stress, force_integrals[m_Stress],
                          IntegrationMethod::simpson);
        add_var_integrand(m_BHMom, force_integrals[m_BHMom],
                          IntegrationMethod::simpson);

        // do the integration over the surface
        integrate();

        // write the integrals
        std::vector<std::string> labels(2);
        labels[m_Stress] = "Force";
        labels[m_BHMom] = "BHMom";
        std::string filename = "Force_integrals";
        write_integrals(filename, force_integrals, labels);
    }
};

#endif /* FORCEEXTRACTION_HPP_ */
