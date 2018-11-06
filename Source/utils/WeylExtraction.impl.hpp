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
    // files so for now it is commented out write_extraction("ExtractionOut_",
    // m_state_ptr_re, m_state_ptr_im);
}

//! integrate over a spherical shell with given harmonics
inline std::array<double, 2>
WeylExtraction::integrate_surface(int es, int el, int em,
                                  const double *m_state_ptr_re,
                                  const double *m_state_ptr_im) const
{
    int rank;
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    std::array<double, 2> integral = {0.0, 0.0};

    // only rank 0 does the integral, but use OMP threads if available
    if (rank == 0)
    {
    // integrate the values over the sphere (normalised by r^2)
    // assumes spacings constant, uses trapezium rule for phi and rectangles for
    // theta  note we don't have to fudge the end points for phi because the
    // function is periodic  and so the last point (implied but not part of
    // vector) is equal to the first point
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
                double x = m_params.extraction_center[0] +
                           m_params.extraction_radius * sin(theta) * cos(phi);
                double y = m_params.extraction_center[1] +
                           m_params.extraction_radius * sin(theta) * sin(phi);
                double z = m_params.extraction_center[2] +
                           m_params.extraction_radius * cos(theta);
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
                                           const char *filename) const
{
    int rank;
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_File mpi_file;
    MPI_File_open(Chombo_MPI::comm, filename,
                  MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,
                  MPI_INFO_NULL, &mpi_file);

    // only rank 0 does the write out
    if (rank == 0)
    {
        int char_length = 60;
        // Header data at first timestep
        if (m_time == m_dt)
        {
            char header[char_length];
            snprintf(header, char_length,
                     "# m_time      integral Re   integral Im \n");
            MPI_File_write(mpi_file, header, strlen(header), MPI_CHAR,
                           MPI_STATUS_IGNORE);
        }
        // Now the data
        char data[char_length];
        snprintf(data, char_length, "%f     %e     %e \n", m_time, integral[0],
                 integral[1]);
        MPI_File_write(mpi_file, data, strlen(data), MPI_CHAR,
                       MPI_STATUS_IGNORE);
    }
    MPI_File_close(&mpi_file);
}

//! Write out the result of the extraction in phi and theta at each timestep
inline void WeylExtraction::write_extraction(char *file_prefix,
                                             const double *m_state_ptr_re,
                                             const double *m_state_ptr_im) const
{
    int rank;
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_File mpi_file;
    char step_str[10];
    sprintf(step_str, "%d", (int)round(m_time / m_dt));
    char filename[20] = "";
    strcat(filename, file_prefix);
    strcat(filename, step_str);
    MPI_File_open(Chombo_MPI::comm, filename,
                  MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,
                  MPI_INFO_NULL, &mpi_file);

    // only rank 0 does the write out
    if (rank == 0)
    {
        // Header data
        char comp_str_re[20];
        char comp_str_im[20];
        sprintf(comp_str_re, UserVariables::variable_names[m_re_comp]);
        sprintf(comp_str_im, UserVariables::variable_names[m_im_comp]);
        char header1[20];
        snprintf(header1, 20, "# m_time : %f \n", m_time);
        MPI_File_write(mpi_file, header1, strlen(header1), MPI_CHAR,
                       MPI_STATUS_IGNORE);
        char header2[70] = "";
        strcat(header2, "# theta     phi     ");
        strcat(header2, comp_str_re);
        strcat(header2, comp_str_im);
        strcat(header2, "/n");
        MPI_File_write(mpi_file, header2, strlen(header2), MPI_CHAR,
                       MPI_STATUS_IGNORE);

        // Now the data
        for (int idx = 0; idx < m_num_points; ++idx)
        {
            int itheta = idx / m_params.num_points_phi;
            int iphi = idx % m_params.num_points_phi;
            // don't put a point at z = 0
            double theta = (itheta + 0.5) * m_dtheta;
            double phi = iphi * m_dphi;
            char data[70];
            snprintf(data, 70, "%f     %f     %e     %e \n", theta, phi,
                     m_state_ptr_re[idx], m_state_ptr_im[idx]);
            MPI_File_write(mpi_file, data, strlen(data), MPI_CHAR,
                           MPI_STATUS_IGNORE);
        }
    }
    MPI_File_close(&mpi_file);
}

#endif /* WEYLEXTRACTION_IMPL_HPP_ */
