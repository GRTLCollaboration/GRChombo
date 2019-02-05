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

    // calculate the integral over the surface -currently only 20, 21 and 22
    // modes implemented here, but easily extended if more required
    auto integrals20 =
        integrate_surface(-2, 2, 0, interp_re_part, interp_im_part);
    auto integrals21 =
        integrate_surface(-2, 2, 1, interp_re_part, interp_im_part);
    auto integrals22 =
        integrate_surface(-2, 2, 2, interp_re_part, interp_im_part);

    // Output the result to a single file over the whole run
    write_integral(integrals20.first, integrals20.second, "Weyl_integral_20");
    write_integral(integrals21.first, integrals21.second, "Weyl_integral_21");
    write_integral(integrals22.first, integrals22.second, "Weyl_integral_22");

    // This generates a lot of output, and it is usually better to output plot
    // files so for now it is commented out
    // write_extraction("ExtractionOut_", interp_re_part,
    //                  interp_im_part);
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
#pragma omp parallel for collapse(2) default(none)                             \
    shared(es, el, em, integral_re, integral_im)
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
    int rank;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
    rank = 0;
#endif
    // only rank 0 does the write out
    if (rank == 0)
    {
        std::ofstream outfile;
        // overwrite file if this is the first timestep, otherwise append.
        if (m_time == m_dt)
        {
            outfile.open(a_filename);
        }
        else
        {
            outfile.open(a_filename, std::ios_base::app);
        }
        if (!outfile)
        {
            MayDay::Error(
                "WeylExtraction::write_integral: error opening output file");
        }

        // Header data at first timestep
        if (m_time == m_dt)
        {
            outfile << "#" << std::setw(9) << "time";
            for (int iradius = 0; iradius < m_params.num_extraction_radii;
                 ++iradius)
            {
                outfile << std::setw(20) << "integral Re";
                outfile << std::setw(20) << "integral Im";
            }
            outfile << "\n";
            outfile << "#" << std::setw(9) << "r =";
            for (int iintegral = 0;
                 iintegral < 2 * m_params.num_extraction_radii; ++iintegral)
            {
                outfile << std::setw(20)
                        << m_params.extraction_radii[iintegral / 2];
            }
            outfile << "\n";
        }

        // Now the data
        outfile << std::fixed << std::setw(10) << m_time;
        outfile << std::scientific << std::setprecision(10);
        for (int iradius = 0; iradius < m_params.num_extraction_radii;
             ++iradius)
        {
            outfile << std::setw(20) << a_integral_re[iradius];
            outfile << std::setw(20) << a_integral_im[iradius];
        }
        outfile << std::endl;
        outfile.close();
    }
}

//! Write out the result of the extraction in phi and theta at each timestep for
//! each extraction radius
inline void
WeylExtraction::write_extraction(std::string a_file_prefix,
                                 const std::vector<double> a_re_part,
                                 const std::vector<double> a_im_part) const
{
    CH_TIME("WeylExtraction::write_extraction");
    int rank;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
    rank = 0;
#endif
    // only rank 0 does the write out
    if (rank == 0)
    {
        // set up file names and component names
        int step = std::round(m_time / m_dt);
        std::string file_str = a_file_prefix + std::to_string(step);
        std::string comp_str_re = UserVariables::variable_names[m_re_comp];
        std::string comp_str_im = UserVariables::variable_names[m_im_comp];

        // write out extraction data to a separate file at each step
        std::ofstream outfile;
        outfile.open(file_str);
        if (!outfile)
        {
            MayDay::Error(
                "WeylExtraction::write_extraction: error opening output file");
        }

        // header data
        for (int iradius = 0; iradius < m_params.num_extraction_radii;
             ++iradius)
        {
            outfile << "# time : " << m_time << ", r = ";
            outfile << m_params.extraction_radii[iradius] << "\n";
            outfile << "#" << std::setw(11) << "theta";
            outfile << std::setw(12) << "phi";
            outfile << std::setw(20) << comp_str_re;
            outfile << std::setw(20) << comp_str_im << "\n";

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
                outfile << std::fixed << std::setprecision(7);
                outfile << std::setw(12) << theta;
                outfile << std::setw(12) << phi;
                outfile << std::scientific << std::setprecision(10);
                outfile << std::setw(20) << a_re_part[idx];
                outfile << std::setw(20) << a_im_part[idx] << "\n";
            }
            outfile << "\n\n";
        }
        outfile.close();
    }
}

#endif /* WEYLEXTRACTION_IMPL_HPP_ */
