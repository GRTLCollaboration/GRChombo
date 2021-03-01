/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _APPARENTHORIZON_HPP_
#define _APPARENTHORIZON_HPP_

#include <petsc.h>
// #include <petscviewerhdf5.h>

#include "AHData.hpp"
#include "AHDeriv.hpp"
#include "AHFinder.hpp"
#include "AHGeometryData.hpp"
#include "AHInterpolation.hpp"
#include "IntegrationMethod.hpp"

// Class to manage ApparentHorizon for 2+1D and 3+1D simulations
// Class AHFinder manages it
//! AHFunction defines the optimizing function (see AHFunction.hpp for
//! expansion example calculation)
template <class SurfaceGeometry, class AHFunction> class ApparentHorizon
{
    using Interpolation = AHInterpolation<SurfaceGeometry, AHFunction>;

  public:
    //! AH that finds the zero of expansion
    ApparentHorizon(
        const Interpolation &a_interp, //!< Geometry class to exchange data
        double a_initial_guess, //!< Initial guess for radius (or whatever
                                //!< coordinate you're solving for)
        const AHFinder::params &a_params, //!< set of AH parameters
        const std::string &a_stats =
            "stats.dat", //!< name for output file with area, spin and AH origin
        const std::string &a_coords =
            "coords_", //!< name for output file with AH coordinates at each
                       //!< time step
        bool solve_first_step = true //!< whether or not to solve if t=0
    );
    //! personalized optimizer that finds zero of function
    //! 'a_function_to_optimize' (a void* 'a_function_to_optimize_params' can be
    //! passed for auxiliary parameters passed to 'a_function_to_optimize')
    ApparentHorizon(const Interpolation &a_interp, double a_initial_guess,
                    const AHFinder::params &a_params,
                    const typename AHFunction::params &a_func_params,
                    const std::string &a_stats = "stats",
                    const std::string &a_coords = "coords",
                    bool solve_first_step = true);
    ~ApparentHorizon();

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
    double get_initial_guess() const;
    void set_initial_guess(double a_initial_guess);
    const Interpolation &get_ah_interp() const;

    // set origin to whatever you want (e.g. punctures) before solving if you
    // want, otherwise we use the last center or the estimate next center
    void set_origin(const std::array<double, CH_SPACEDIM> &a_origin);

    double get_max_F() const;
    double get_min_F() const;
    double get_ave_F() const;
    double get_std_F() const;

    bool do_solve(double a_dt, double a_time)
        const; //!< decide (based times passed to 'solve') whether or not
               //!< to print (uses AHFinder::params::solve_interval)
    bool do_print(double a_dt, double a_time)
        const; //!< decide when to print (only AHFinder::params::print_interval
               //!< out of all 'solve's)

    const AHFinder::params &m_params; //!< set of AH parameters
    const std::string m_stats,
        m_coords; //!< public base names for output files (no need for a set as
                  //!< they are const)

    // any parameters that want to be saved to be passed to optimization
    // function
    typename AHFunction::params m_func_params;

  private:
    void write_outputs(double a_dt, double a_time, double a_restart_time);

    void
    restart(bool solve_first_step); //!< restart AH, updating the coordinates
                                    //!< based on the last output file and the
                                    //!< origin based on the stats file

    double calculate_area(); //!< calculate AH area

#if CH_SPACEDIM == 3
    Tensor<1, double> calculate_angular_momentum_J(); //!< calculate spin, ONLY
                                                      //!< FOR 3D
    double calculate_spin_dimensionless(
        double a_area); //!< calculate spin with 'z' direction, ONLY FOR 3D
#endif
    // estimate based on area and angular momentum J
    ALWAYS_INLINE double calculate_mass(double area, double J_norm)
    {
        return sqrt(area / (16. * M_PI) + J_norm * J_norm * 4. * M_PI / area);
    }
    ALWAYS_INLINE double calculate_irreducible_mass(double area)
    {
        return calculate_mass(area, 0.);
    }

    void calculate_minmax_F() const;
    void calculate_average_F() const;

    void reset_initial_guess();

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

    void initialise_PETSc(); //!< initialise automatically done in constructor
    void finalise_PETSc();   //!< finalise automatically done in destructor

    // variables
  private:
    bool m_printed_once;

    int m_converged; //!< flag saying if PETSc has converged the last 'N' times
                     //!<(read using 'get_converged()')
    bool m_has_been_found; //!< flag saying if an horizon has ever been found
                           //!< (read using 'has_been_found()')
    int m_num_failed_convergences; //!< the number of failed consecutive
                                   //!< convergences

    double m_initial_guess; //!< initial guess for AH (saved so that it can be
                            //!< re-used when atempting to solve again)

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
    double m_area, m_spin, m_mass, m_irreducible_mass, m_spin_z_alt;
    Tensor<1, double> m_dimensionless_spin_vector;

    // prevents resetting the origin when the user externally did 'set_origin'
    // before 'solve'
    bool origin_already_updated;

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////// PETSc stuff below ///////////////////////////
    /////////////////////////////////////////////////////////////////////////

#if CH_SPACEDIM == 3
    typedef PetscScalar **dmda_arr_t;
#elif CH_SPACEDIM == 2
    typedef PetscScalar *dmda_arr_t;
#endif

    //! Geometries of the AH
    //! 'm_geom_plus' and 'm_geom_minus' are used to calculate the
    //! jacobian of the expansion using a 'delta' numerical differentiation
    Interpolation m_interp;
    Interpolation m_interp_plus;
    Interpolation m_interp_minus;

    //!< used to compute jacobian of expansion (numerical differentiation)
    static constexpr double eps = 1e-7;

    const bool m_periodic_u; //!< is 'u' periodic?
    PetscInt m_num_global_u; //!< total number of grid points in 'u' coordinate
    double m_du;             //!< physical 'delta' in 'u' coordinate

#if CH_SPACEDIM == 3
    const bool m_periodic_v; //!< is 'v' periodic?
    PetscInt m_num_global_v; //!< total number of grid points in 'u' coordinate
    double m_dv;             //!< physical 'delta' in 'v' coordinate
#endif

    //! vectors to store and manipulate 'F', u', and 'v'
    //! internally, in interaction with the 'Interpolation's
    std::vector<double> m_F;
    std::vector<double> m_u;
#if CH_SPACEDIM == 3
    std::vector<double> m_v;
#endif

    //! minimums and maximums of coordinates 'u' and 'v'
    //! of the PETSc grid specific to the current rank
    PetscInt m_umin;
    PetscInt m_umax;

#if CH_SPACEDIM == 3
    PetscInt m_vmin;
    PetscInt m_vmax;
#endif

    //! number of points in 'u' and 'v' direction
    //! (m_nu = m_umax - m_umin)
    PetscInt m_nu;
#if CH_SPACEDIM == 3
    PetscInt m_nv;
#endif

    // PETSc main object
    DM m_dmda;
    //! Scalable Nonlinear Equations Solvers
    SNES m_snes;

    Vec m_snes_soln;
    Vec m_snes_rhs;
    Mat m_snes_jac;

    //! interpolate (u,v) 2D grid points at restart if number of points in
    //! either direction changed
    //! This only interpolates the points that the PETSc rank that called it has
    //! returns whether or not points were interpolated
    bool interpolate_ah(dmda_arr_t f,
                        const std::vector<std::vector<double>> &old_coords);

    //! set the default stencils of AHDeriv at position {u,v}
    void set_stencils(AHDeriv &out, int u
#if CH_SPACEDIM == 3
                      ,
                      int v
#endif
    );

    //! function to calculate 1st and 2nd derivatives of 'in'
    //! (tipically corresponds to our 'f' function)
    //! in the 'u' and 'v' directions
    AHDeriv diff(const dmda_arr_t in, int u
#if CH_SPACEDIM == 3
                 ,
                 int v
#endif
    );

    //! private functions used to compute the RHS (the expansion) and it's
    //! jacobian
    void form_function(Vec F, Vec Rhs);
    void form_jacobian(Vec F, Mat J);

    //! helper for 'form_jacobian'
    double point_jacobian(int u, int u_stencil,
#if CH_SPACEDIM == 3
                          int v, int v_stencil,
#endif
                          dmda_arr_t in, int idx,
                          const Interpolation &interp_plus,
                          const Interpolation &interp_minus);

    //! functions used by PETSc based on 'form_function' and 'form_jacobian'
    static PetscErrorCode Petsc_form_function(SNES snes, Vec F, Vec Rhs,
                                              void *ptr);

    static PetscErrorCode
#if PETSC_VERSION_GE(3, 5, 0)
    Petsc_form_jacobian(SNES snes, Vec F, Mat Amat, Mat Pmat, void *ptr);
#else
    Petsc_form_jacobian(SNES snes, Vec F, Mat *Amat, Mat *Pmat,
                        MatStructure *flag, void *ptr);
#endif

    // monitor function required for SNES
    static PetscErrorCode Petsc_SNES_monitor(SNES snes, PetscInt its,
                                             PetscReal norm, void *ptr);
};

#include "ApparentHorizon.impl.hpp"
#include "ApparentHorizon_petsc.impl.hpp"

#endif /* _APPARENTHORIZON_HPP_ */
