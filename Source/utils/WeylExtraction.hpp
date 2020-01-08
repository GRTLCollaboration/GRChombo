/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WEYLEXTRACTION_HPP_
#define WEYLEXTRACTION_HPP_

#include "SphericalExtraction.hpp"
//!  The class allows extraction of the values of the Weyl scalar components on
//!  spherical shells at specified radii, and integration over those shells
/*!
   The class allows the user to extract data from the grid for the Weyl
   components over spherical shells at specified radii. The values may then be
   written to an output file, or integrated across the surfaces.
*/
class WeylExtraction : public SphericalExtraction
{
  public:
    //! The constructor
    WeylExtraction(SphericalExtraction::params_t &a_params, double a_dt,
                   double a_time, bool a_first_step,
                   double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time)
    {
        add_var(c_Weyl4_Re);
        add_var(c_Weyl4_Im);
    }

    //! The old constructor which assumes it is called in specificPostTimeStep
    //! so the first time step is when m_time == m_dt
    WeylExtraction(SphericalExtraction::params_t a_params, double a_dt,
                   double a_time, double a_restart_time = 0.0)
        : WeylExtraction(a_params, a_dt, a_time, (a_dt == a_time),
                         a_restart_time)
    {
    }

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        // extract the values of the Weyl scalars on the spheres
        extract(a_interpolator);

        if (m_params.write_extraction)
        {
            write_extraction("Weyl4ExtractionOut_");
        }

        // now calculate and write the requested spherical harmonic modes
        for (int imode = 0; imode < m_num_modes; ++imode)
        {
            auto normalised_Weyl4_complex =
                [](std::vector<double> Weyl4_reim_parts, double r, double,
                   double) {
                    // here the std::vector<double> passed will just have
                    // the real and imaginary parts of the Weyl4 scalar as its
                    // only components
                    return std::make_pair(r * Weyl4_reim_parts[0],
                                          r * Weyl4_reim_parts[1]);
                };

            const auto &mode = m_modes[imode];
            constexpr int es = -2;
            auto integrals = integrate_mode(es, mode.first, mode.second,
                                            normalised_Weyl4_complex);
            std::string integrals_filename = "Weyl_integral_" +
                                             std::to_string(mode.first) +
                                             std::to_string(mode.second);
            std::vector<std::vector<double>> integrals_for_writing = {
                std::move(integrals.first), std::move(integrals.second)};
            std::vector<std::string> labels = {"integral Re", "integral Im"};
            write_integrals(integrals_filename, integrals_for_writing, labels);
        }
    }
};

#endif /* WEYLEXTRACTION_HPP_ */
