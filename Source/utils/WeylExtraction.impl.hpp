/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(WEYLEXTRACTION_HPP_)
#error "This file should only be included through WeylExtraction.hpp"
#endif

#ifndef WEYLEXTRACTION_IMPL_HPP_
#define WEYLEXTRACTION_IMPL_HPP_

//! Set up tnd execute the interpolation query
inline void WeylExtraction::execute_query(
    AMRInterpolator<Lagrange<4>> *m_interpolator) const
{
    if (m_interpolator == nullptr)
    {
        MayDay::Error("Interpolator has not been initialised in GRAMR class.");
    }

    std::unique_ptr<double[]> m_state_ptr_re{new double[m_num_points]};
    std::unique_ptr<double[]> m_state_ptr_im{new double[m_num_points]};
    std::unique_ptr<double[]> m_interp_x{new double[m_num_points]};
    std::unique_ptr<double[]> m_interp_y{new double[m_num_points]};
    std::unique_ptr<double[]> m_interp_z{new double[m_num_points]};

    // Work out the coordinates
    for (int idx = 0; idx < m_num_points; ++idx)
    {
        int itheta = idx / m_params.num_points_phi;
        int iphi = idx % m_params.num_points_phi;
        // don't put a point at z = 0
        double theta = (itheta + 0.5) * m_dtheta;
        double phi = iphi * m_dphi;
        m_interp_x[idx] = m_params.extraction_center[0] +
                          m_params.extraction_radius * sin(theta) * cos(phi);
        m_interp_y[idx] = m_params.extraction_center[1] +
                          m_params.extraction_radius * sin(theta) * sin(phi);
        m_interp_z[idx] = m_params.extraction_center[2] +
                          m_params.extraction_radius * cos(theta);
    }

    // set up the query
    InterpolationQuery query(m_num_points);
    query.setCoords(0, m_interp_x.get())
        .setCoords(1, m_interp_y.get())
        .setCoords(2, m_interp_z.get())
        .addComp(m_re_comp, m_state_ptr_re.get())
        .addComp(m_im_comp, m_state_ptr_im.get());

    // submit the query
    m_interpolator->interp(query);

    // clear the memory for the coordinates
    m_interp_x.reset();
    m_interp_y.reset();
    m_interp_z.reset();

    // calculate the integral over the surface -currently only 20, 21 and 22
    // modes implemented here, but easily extended if more required
    std::array<double, 2> integral20 = {0, 0};
    std::array<double, 2> integral21 = {0, 0};
    std::array<double, 2> integral22 = {0, 0};
    integral20 =
        integrate_surface(-2, 2, 0, m_state_ptr_re.get(), m_state_ptr_im.get());
    integral21 =
        integrate_surface(-2, 2, 1, m_state_ptr_re.get(), m_state_ptr_im.get());
    integral22 =
        integrate_surface(-2, 2, 2, m_state_ptr_re.get(), m_state_ptr_im.get());

    // Output the result to a single file over the whole run
    write_integral(integral20, "Weyl_integral_20");
    write_integral(integral21, "Weyl_integral_21");
    write_integral(integral22, "Weyl_integral_22");

    // This generates a lot of output, and it is usually better to output plot
    // files so for now it is commented out
    // write_extraction("ExtractionOut_", m_state_ptr_re.get(),
    //                  m_state_ptr_im.get());
}

//! integrate over a spherical shell with given harmonics
inline std::array<double, 2>
WeylExtraction::integrate_surface(int es, int el, int em,
                                  const double *m_state_ptr_re,
                                  const double *m_state_ptr_im) const
{
    int rank;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
    rank = 0;
#endif
    std::array<double, 2> integral = {0.0, 0.0};

    // only rank 0 does the integral, but use OMP threads if available
    if (rank == 0)
    {
        // integrate the values over the sphere (normalised by r^2)
        // assumes spacings constant, uses trapezium rule for phi and rectangles
        // for theta  note we don't have to fudge the end points for phi because
        // the function is periodic  and so the last point (implied but not part
        // of vector) is equal to the first point
#pragma omp parallel for
        for (int iphi = 0; iphi < m_params.num_points_phi; ++iphi)
        {
            double phi = iphi * 2 * M_PI / m_params.num_points_phi;
            double inner_integral_re = 0;
            double inner_integral_im = 0;
            for (int itheta = 0; itheta < m_params.num_points_theta; itheta++)
            {
                using namespace SphericalHarmonics;
                double theta = (itheta + 0.5) * m_dtheta;
                int idx = itheta * m_params.num_points_phi + iphi;
                double x = m_params.extraction_radius * sin(theta) * cos(phi);
                double y = m_params.extraction_radius * sin(theta) * sin(phi);
                double z = m_params.extraction_radius * cos(theta);
                Y_lm_t<double> Y_lm = spin_Y_lm(x, y, z, es, el, em);
                double integrand_re = m_state_ptr_re[idx] * Y_lm.Real +
                                      m_state_ptr_im[idx] * Y_lm.Im;
                double integrand_im = m_state_ptr_im[idx] * Y_lm.Real -
                                      m_state_ptr_re[idx] * Y_lm.Im;
                double f_theta_phi_re = integrand_re * sin(theta);
                double f_theta_phi_im = integrand_im * sin(theta);
                inner_integral_re += m_dtheta * f_theta_phi_re;
                inner_integral_im += m_dtheta * f_theta_phi_im;
            }
            integral[0] += m_dphi * inner_integral_re;
            integral[1] += m_dphi * inner_integral_im;
        }
    }
    return integral;
}

//! Write out calculated value of integral
inline void WeylExtraction::write_integral(std::array<double, 2> integral,
                                           std::string filename) const
{
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
            outfile.open(filename);
        }
        else
        {
            outfile.open(filename, std::ios_base::app);
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
            outfile << std::setw(20) << "integral Re";
            outfile << std::setw(20) << "integral Im" << std::endl;
        }

        // Now the data
        outfile << std::fixed << std::setw(10) << m_time;
        outfile << std::scientific << std::setw(20) << std::setprecision(10);
        outfile << integral[0] << std::setw(20) << integral[1] << std::endl;
        outfile.close();
    }
}

//! Write out the result of the extraction in phi and theta at each timestep
inline void WeylExtraction::write_extraction(std::string file_prefix,
                                             const double *m_state_ptr_re,
                                             const double *m_state_ptr_im) const
{
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
        std::string file_str = file_prefix + std::to_string(step);
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
        outfile << "# time : " << m_time << "\n";
        outfile << "#" << std::setw(11) << "theta";
        outfile << std::setw(12) << "phi";
        outfile << std::setw(20) << comp_str_re;
        outfile << std::setw(20) << comp_str_im << "\n";

        // Now the data
        for (int idx = 0; idx < m_num_points; ++idx)
        {
            int itheta = idx / m_params.num_points_phi;
            int iphi = idx % m_params.num_points_phi;
            // don't put a point at z = 0
            double theta = (itheta + 0.5) * m_dtheta;
            double phi = iphi * m_dphi;
            outfile << std::fixed << std::setprecision(7);
            outfile << std::setw(12) << theta;
            outfile << std::setw(12) << phi;
            outfile << std::scientific << std::setprecision(10);
            outfile << std::setw(20) << m_state_ptr_re[idx];
            outfile << std::setw(20) << m_state_ptr_im[idx] << "\n";
        }
        outfile.close();
    }
}

#endif /* WEYLEXTRACTION_IMPL_HPP_ */
