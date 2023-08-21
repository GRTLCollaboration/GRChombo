/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _APPARENTHORIZON_HPP_
#define _APPARENTHORIZON_HPP_

#include "AHInterpolation.hpp"
#include "AHParams.hpp"
#include "IntegrationMethod.hpp"
#include "PETScAHSolver.hpp"

// Chombo namespace
#include "UsingNamespace.H"

// Class to manage ApparentHorizon for 2+1D and 3+1D simulations
// (has only been implemented for 3+1 spacetimes, 3+1 cartoon-reduced to 2+1 and
// 4+1 cartoon-reduced to 2+1)
//! AHFunction defines the optimizing function (see AHFunction.hpp for
//! expansion example calculation)
template <class SurfaceGeometry, class AHFunction> class ApparentHorizon
{
    using AHInterpolation = AHInterpolation_t<SurfaceGeometry, AHFunction>;
    using AHParams = AHParams_t<AHFunction>;

  public:
    //! AH that finds the zero of AHFunction
    ApparentHorizon(
        const AHInterpolation &a_interp, //!< Geometry class to exchange data
        const AHInitialGuessPtr
            &a_initial_guess,     //!< Initial guess for radius (or whatever
                                  //!< coordinate you're solving for)
        const AHParams &a_params, //!< set of AH parameters
        const std::string &a_stats =
            "stats", //!< name for output file with area, spin and AH origin
        const std::string &a_coords =
            "coords_", //!< name for output file with AH coordinates at each
                       //!< time step
        bool solve_first_step = true //!< whether or not to solve if t=0
    );

    void solve(double a_dt, double a_time, double a_restart_time);

    bool good_to_go(double a_dt, double a_time) const;
    bool get_converged()
        const; //!< PETSc didn't converge last time solve() was called
    int get_failed_convergences()
        const; //!< PETSc didn't converge last time solve() was called
    bool has_been_found() const; //!< PETSc converged once and then
                                 //!< stopped converging (AH collapsed)

    const std::array<double, CH_SPACEDIM> &get_origin() const;
    const std::array<double, CH_SPACEDIM> &get_center() const;
    const AHInterpolation &get_ah_interp() const;
    PETScAHSolver<SurfaceGeometry, AHFunction> &get_petsc_solver();

    // set origin to whatever you want (e.g. punctures) before solving if you
    // want, otherwise we use the last center or the estimate next center
    void set_origin(const std::array<double, CH_SPACEDIM> &a_origin);

    double get_max_F() const;
    double get_min_F() const;
    double get_ave_F() const;
    double get_std_F() const;

    bool do_solve(double a_dt, double a_time)
        const; //!< decide (based times passed to 'solve') whether or not
               //!< to print (uses params::solve_interval)
    bool do_print(double a_dt, double a_time)
        const; //!< decide when to print (only params::print_interval
               //!< out of all 'solve's)

    // variables
  public:
    const AHParams m_params; //!< set of AH parameters
    const std::string m_stats,
        m_coords; //!< public base names for output files (no need for a set as
                  //!< they are const)

  private:
    void write_outputs(double a_dt, double a_time, double a_restart_time);

    void
    restart(bool solve_first_step); //!< restart AH, updating the coordinates
                                    //!< based on the last output file and the
                                    //!< origin based on the stats file

    // compute area, linear momentum and angular momentum of BH  (P only for 3D,
    // J only for 3D without cartoon)
    void calculate_ah_quantities(double &area
#if GR_SPACEDIM == 3 // GR_SPACEDIM, not CH_SPACEDIM !!!
                                 ,
                                 Tensor<1, double> &P
#if CH_SPACEDIM == 3
                                 ,
                                 Tensor<1, double> &J
#endif
#endif
    );

#if GR_SPACEDIM == 3 // GR_SPACEDIM, not CH_SPACEDIM !!!
    // estimate based on area and angular momentum J
    ALWAYS_INLINE static double calculate_mass(double area, double J_norm)
    {
        return sqrt(area / (16. * M_PI) + J_norm * J_norm * 4. * M_PI / area);
    }

#if CH_SPACEDIM == 3
    ALWAYS_INLINE static double calculate_irreducible_mass(double area)
    {
        return calculate_mass(area, 0.);
    }

    double calculate_spin_dimensionless(
        double a_area); //!< calculate spin with 'z' direction, ONLY FOR 3D
#endif
#endif

    void calculate_minmax_F() const;
    void calculate_average_F() const;

    void update_old_centers(std::array<double, CH_SPACEDIM>);
    void predict_next_origin();
    std::array<double, CH_SPACEDIM>
    calculate_center(); //!< update location of center by calculating the
                        //!< centroid of the AH

    void write_coords_file(double a_dt, double a_time, double a_restart_time,
                           const std::string &filename,
                           bool write_geometry_data = false)
        const; //!< write coords in (m_u, m_v, m_F) to 'filename'

    //! write metric and extrinsic curvature for each point of AH
    void write_geometry_data(const std::string &a_file_prefix, double a_dt,
                             double a_time, double a_restart_time) const;

    void check_convergence(); //!< check if PETSc has converged and share that
                              //!< info across all Chombo processes

    // default to simpson rule and change if invalid
    void check_integration_methods();

    // variables
  private:
    bool m_printed_once, m_printed_after_restart;

    int m_converged; //!< flag saying if PETSc has converged the last 'N' times
                     //!<(read using 'get_converged()')
    bool m_has_been_found; //!< flag saying if an horizon has ever been found
                           //!< (read using 'has_been_found()')
    int m_num_failed_convergences; //!< the number of failed consecutive
                                   //!< convergences

    //! used to estimate where the new center will be before solving
    //! (currently last 2 centers saved to make a parabolic regression)
    //! (elements ordered in reverse time - older with higher index)
    std::vector<std::array<double, CH_SPACEDIM>> m_old_centers;

    // mutable so that 'const' objects can still change them
    //! stores the global maximum and minimum of F - calculate with
    //! calculate_minmax_F
    mutable double m_max_F, m_min_F, m_ave_F, m_std_F;

    // Integration method to calculate area and spin of AH. Defaults to Simpson
    // method
    std::array<IntegrationMethod, CH_SPACEDIM - 1> m_integration_methods;

    // just to save the result temporarily at each iteration
    double m_area;
#if GR_SPACEDIM == 3
    double m_mass, m_linear_momentum_P_norm;
#if CH_SPACEDIM == 3
    double m_irreducible_mass, m_spin, m_spin_z_alt;
    Tensor<1, double> m_dimensionless_spin_vector, m_linear_momentum_P;
#endif
#endif

    // prevents resetting the origin when the user externally did 'set_origin'
    // before 'solve'
    bool origin_already_updated;

    PETScAHSolver<SurfaceGeometry, AHFunction> solver;
};

#include "ApparentHorizon.impl.hpp"

#endif /* _APPARENTHORIZON_HPP_ */
