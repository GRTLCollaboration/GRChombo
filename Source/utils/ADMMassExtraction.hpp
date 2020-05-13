/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ADMMASSEXTRACTION_HPP_
#define ADMMASSEXTRACTION_HPP_

#include "SphericalExtraction.hpp"
//!  The class allows extraction of the values of the ADM mass and angular
//!  momentum over spherical shells at specified radii, and integration over
//!  those shells
/*!
   The class allows the user to extract data from the grid for the ADM mass and
   angular momentum over spherical shells at specified radii. The values may
   then be written to an output file, or integrated across the surfaces.
*/
class ADMMassExtraction : public SphericalExtraction
{
  public:
    //! The constructor
    ADMMassExtraction(SphericalExtraction::params_t &a_params, double a_dt,
                      double a_time, bool a_first_step,
                      double a_restart_time = 0.0, bool a_extract_J = true)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time),
          extract_J(a_extract_J)
    {
        add_var(c_Madm);
        if (extract_J)
            add_var(c_Jadm);
    }

    //! The old constructor which assumes it is called in specificPostTimeStep
    //! so the first time step is when m_time == m_dt
    ADMMassExtraction(SphericalExtraction::params_t a_params, double a_dt,
                      double a_time, double a_restart_time = 0.0,
                      bool a_extract_J = true)
        : ADMMassExtraction(a_params, a_dt, a_time, (a_dt == a_time),
                            a_restart_time, a_extract_J)
    {
    }

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        // extract the values of the ADM mass and spin on the spheres
        extract(a_interpolator);

        if (m_params.write_extraction)
        {
            write_extraction("MadmExtractionOut_");
        }

        // now calculate and write the requested spherical harmonic modes
        std::vector<std::vector<double>> out_integrals(1 + extract_J);

        add_var_integrand(0, out_integrals[0], IntegrationMethod::simpson);
        if (extract_J)
            add_var_integrand(1, out_integrals[1], IntegrationMethod::simpson);

        // do the integration over the surface
        integrate();
        int extrapolation_order = 2;
        std::vector<std::string> labels(1 + extract_J);
        labels[0] = "M_adm";
        if (extract_J)
            labels[1] = "J_adm";
        write_integrals("IntegralADMmass", out_integrals, labels,
                        extrapolation_order);
    }

  private:
    bool extract_J;
};

#endif /* ADMMASSEXTRACTION_HPP_ */
