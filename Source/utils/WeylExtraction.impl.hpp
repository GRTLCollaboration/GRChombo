/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(WEYLEXTRACTION_HPP_)
#error "This file should only be included through WeylExtraction.hpp"
#endif

#ifndef WEYLEXTRACTION_IMPL_HPP_
#define WEYLEXTRACTION_IMPL_HPP_

//! Set up and execute the interpolation query
inline void WeylExtraction::execute_query(
    AMRInterpolator<Lagrange<4>> *a_interpolator) const
{
    CH_TIME("WeylExtraction::execute_query");
    if (a_interpolator == nullptr)
    {
        MayDay::Error("Interpolator has not been initialised in GRAMR class.");
    }

    std::vector<double> interp_re_part(m_num_points *
                                       m_params.num_extraction_radii);
    std::vector<double> interp_im_part(m_num_points *
                                       m_params.num_extraction_radii);
    std::vector<double> interp_x(m_num_points * m_params.num_extraction_radii);
    std::vector<double> interp_y(m_num_points * m_params.num_extraction_radii);
    std::vector<double> interp_z(m_num_points * m_params.num_extraction_radii);

    // Work out the coordinates
    for (int iradius = 0; iradius < m_params.num_extraction_radii; ++iradius)
    {
        for (int idx = 0; idx < m_num_points; ++idx)
        {
            int itheta = idx / m_params.num_points_phi;
            int iphi = idx % m_params.num_points_phi;
            // don't put a point at z = 0
            double theta = (itheta + 0.5) * m_dtheta;
            double phi = iphi * m_dphi;
            interp_x[iradius * m_num_points + idx] =
                m_params.extraction_center[0] +
                m_params.extraction_radii[iradius] * sin(theta) * cos(phi);
            interp_y[iradius * m_num_points + idx] =
                m_params.extraction_center[1] +
                m_params.extraction_radii[iradius] * sin(theta) * sin(phi);
            interp_z[iradius * m_num_points + idx] =
                m_params.extraction_center[2] +
                m_params.extraction_radii[iradius] * cos(theta);
        }
    }
    // set up the query
    InterpolationQuery query(m_num_points * m_params.num_extraction_radii);
    query.setCoords(0, interp_x.data())
        .setCoords(1, interp_y.data())
        .setCoords(2, interp_z.data())
        .addComp(m_re_comp, interp_re_part.data())
        .addComp(m_im_comp, interp_im_part.data());

    // submit the query
    a_interpolator->interp(query);

    for (int imode = 0; imode < m_params.num_modes; ++imode)
    {
        const std::pair<int, int> &mode = m_params.modes[imode];
        auto integral = integrate_surface(-2, mode.first, mode.second,
                                          interp_re_part, interp_im_part);
        std::string integral_filename = "Weyl_integral_" +
                                        std::to_string(mode.first) +
                                        std::to_string(mode.second);
        write_integral(integral.first, integral.second, integral_filename);
    }

    if (m_params.write_extraction)
    {
        write_extraction("ExtractionOut_", interp_re_part, interp_im_part);
    }
}

//! integrate over a spherical shell with given harmonics for each extraction
//! radius and normalise by multiplying by radius
inline std::pair<std::vector<double>, std::vector<double>>
WeylExtraction::integrate_surface(int es, int el, int em,
                                  const std::vector<double> a_re_part,
                                  const std::vector<double> a_im_part) const
{
    CH_TIME("WeylExtraction::integrate_surface");
    int rank;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
    rank = 0;
#endif
    std::vector<double> integral_re(m_params.num_extraction_radii, 0.);
    std::vector<double> integral_im(m_params.num_extraction_radii, 0.);

    // only rank 0 does the integral, but use OMP threads if available
    if (rank == 0)
    {
        // integrate the values over the sphere (normalised by r^2)
        // for each radius
        // assumes spacings constant, uses trapezium rule for phi and rectangles
        // for theta  note we don't have to fudge the end points for phi because
        // the function is periodic  and so the last point (implied but not part
        // of vector) is equal to the first point
#ifdef _OPENMP
#if __GNUC__ > 8
#define OPENMP_CONST_SHARED shared(a_re_part, a_im_part)
#else
#define OPENMP_CONST_SHARED
#endif
#pragma omp parallel for collapse(2) default(none)                             \
    shared(es, el, em, integral_re, integral_im) OPENMP_CONST_SHARED
#undef OPENMP_CONST_SHARED
#endif
        for (int iradius = 0; iradius < m_params.num_extraction_radii;
             ++iradius)
        {
            for (int iphi = 0; iphi < m_params.num_points_phi; ++iphi)
            {
                double phi = iphi * 2 * M_PI / m_params.num_points_phi;
                double inner_integral_re = 0.;
                double inner_integral_im = 0.;
                for (int itheta = 0; itheta < m_params.num_points_theta;
                     itheta++)
                {
                    using namespace SphericalHarmonics;
                    double theta = (itheta + 0.5) * m_dtheta;
                    int idx = iradius * m_num_points +
                              itheta * m_params.num_points_phi + iphi;
                    double x = m_params.extraction_radii[iradius] * sin(theta) *
                               cos(phi);
                    double y = m_params.extraction_radii[iradius] * sin(theta) *
                               sin(phi);
                    double z = m_params.extraction_radii[iradius] * cos(theta);
                    Y_lm_t<double> Y_lm = spin_Y_lm(x, y, z, es, el, em);
                    double integrand_re =
                        a_re_part[idx] * Y_lm.Real + a_im_part[idx] * Y_lm.Im;
                    double integrand_im =
                        a_im_part[idx] * Y_lm.Real - a_re_part[idx] * Y_lm.Im;

                    // note the multiplication by radius here
                    double f_theta_phi_re = m_params.extraction_radii[iradius] *
                                            integrand_re * sin(theta);
                    double f_theta_phi_im = m_params.extraction_radii[iradius] *
                                            integrand_im * sin(theta);
                    inner_integral_re += m_dtheta * f_theta_phi_re;
                    inner_integral_im += m_dtheta * f_theta_phi_im;
                }
#ifdef _OPENMP
#pragma omp atomic
#endif
                integral_re[iradius] += m_dphi * inner_integral_re;
#ifdef _OPENMP
#pragma omp atomic
#endif
                integral_im[iradius] += m_dphi * inner_integral_im;
            }
        }
    }
    return std::make_pair(integral_re, integral_im);
}

//! Write out calculated value of integral for each extraction radius
inline void
WeylExtraction::write_integral(const std::vector<double> a_integral_re,
                               const std::vector<double> a_integral_im,
                               std::string a_filename) const
{
    CH_TIME("WeylExtraction::write_integral");
    // open file for writing
    SmallDataIO integral_file(a_filename, m_dt, m_time, m_restart_time,
                              SmallDataIO::APPEND, m_first_step);

    // remove any duplicate data if this is a restart
    integral_file.remove_duplicate_time_data();

    // need to write headers if this is the first timestep
    if (m_first_step)
    {
        // make header strings
        std::vector<std::string> header1_strings(2 *
                                                 m_params.num_extraction_radii);
        std::vector<std::string> header2_strings(2 *
                                                 m_params.num_extraction_radii);
        for (int iintegral = 0; iintegral < 2 * m_params.num_extraction_radii;
             iintegral += 2)
        {
            int iintegral1 = iintegral + 1;
            int iradius = iintegral / 2;
            header1_strings[iintegral] = "integral Re";
            header1_strings[iintegral1] = "integral Im";
            header2_strings[iintegral1] = header2_strings[iintegral] =
                std::to_string(m_params.extraction_radii[iradius]);
        }

        // write headers
        integral_file.write_header_line(header1_strings);
        integral_file.write_header_line(header2_strings, "r = ");
    }

    // make vector of data for writing
    std::vector<double> data_for_writing(2 * m_params.num_extraction_radii);
    for (int iintegral = 0; iintegral < 2 * m_params.num_extraction_radii;
         iintegral += 2)
    {
        int iintegral1 = iintegral + 1;
        int iradius = iintegral / 2;
        data_for_writing[iintegral] = a_integral_re[iradius];
        data_for_writing[iintegral1] = a_integral_im[iradius];
    }

    // write data
    integral_file.write_time_data_line(data_for_writing);
}

//! Write out the result of the extraction in phi and theta at each timestep for
//! each extraction radius
inline void
WeylExtraction::write_extraction(std::string a_file_prefix,
                                 const std::vector<double> a_re_part,
                                 const std::vector<double> a_im_part) const
{
    CH_TIME("WeylExtraction::write_extraction");
    SmallDataIO extraction_file(a_file_prefix, m_dt, m_time, m_restart_time,
                                SmallDataIO::NEW, m_first_step);

    for (int iradius = 0; iradius < m_params.num_extraction_radii; ++iradius)
    {
        // Write headers
        std::vector<std::string> header1_strings = {
            "time = " + std::to_string(m_time) + ",",
            "r = " + std::to_string(m_params.extraction_radii[iradius])};
        extraction_file.write_header_line(header1_strings, "");
        std::vector<std::string> components = {
            UserVariables::variable_names[m_re_comp],
            UserVariables::variable_names[m_im_comp]};
        std::vector<std::string> coords = {"theta", "phi"};
        extraction_file.write_header_line(components, coords);

        // Now the data
        for (int idx = iradius * m_num_points;
             idx < (iradius + 1) * m_num_points; ++idx)
        {
            int itheta =
                (idx - iradius * m_num_points) / m_params.num_points_phi;
            int iphi = idx % m_params.num_points_phi;
            // don't put a point at z = 0
            double theta = (itheta + 0.5) * m_dtheta;
            double phi = iphi * m_dphi;

            extraction_file.write_data_line({a_re_part[idx], a_im_part[idx]},
                                            {theta, phi});
        }
        extraction_file.line_break();
    }
}

#endif /* WEYLEXTRACTION_IMPL_HPP_ */
