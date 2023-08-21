/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _PETSCAHSOLVER_HPP_
#define _PETSCAHSOLVER_HPP_

#include <petsc.h>
// #include <petscviewerhdf5.h>

#include "AHDerivData.hpp"
#include "AHInterpolation.hpp"
#include "AHParams.hpp"

#if CH_SPACEDIM == 3
typedef PetscScalar **dmda_arr_t;
#elif CH_SPACEDIM == 2
typedef PetscScalar *dmda_arr_t;
#endif

//! Helper class for ApparentHorizon class that manages all weird PETSc stuff
template <class SurfaceGeometry, class AHFunction> class PETScAHSolver
{
  public:
    using AHInterpolation = AHInterpolation_t<SurfaceGeometry, AHFunction>;
    using AHParams = AHParams_t<AHFunction>;

  public:
    //! AH that finds the zero of expansion
    PETScAHSolver(const AHInterpolation &a_interp,
                  const AHInitialGuessPtr &a_initial_guess,
                  const AHParams &a_params);
    ~PETScAHSolver();

    //! function to calculate 1st and 2nd derivatives of 'in'
    //! (tipically corresponds to our 'f' function)
    //! in the 'u' and 'v' directions
    AHDerivData diff(D_DECL(const dmda_arr_t in, int u, int v));

    //! interpolate (u,v) 2D grid points at restart if number of points in
    //! either direction changed
    //! This only interpolates the points that the PETSc rank that called it has
    //! returns whether or not points were interpolated
    bool interpolate_ah(const std::vector<std::vector<double>> &old_coords);

    void solve();
    SNESConvergedReason getConvergedReason() const;

    void get_dmda_arr_t(Vec &localF, dmda_arr_t &in);
    // must be called at the end of the function that called 'get_dmda_arr_t'
    void restore_dmda_arr_t(Vec &localF, dmda_arr_t &in);

    const AHInitialGuessPtr &get_initial_guess() const;
    void reset_initial_guess();

    const std::array<double, CH_SPACEDIM> &get_origin() const;
    void set_origin(const std::array<double, CH_SPACEDIM> &a_origin);

    void initialise(); //!< initialise automatically done in constructor
    void finalise();   //!< finalise automatically done in destructor

    // variables
  public:
    //! Geometries of the AH
    //! 'm_geom_plus' and 'm_geom_minus' are used to calculate the
    //! jacobian of the expansion using a 'delta' numerical differentiation
    AHInterpolation m_interp;
    AHInterpolation m_interp_plus;
    AHInterpolation m_interp_minus;

    const bool m_periodic_u; //!< is 'u' periodic?
    PetscInt m_num_global_u; //!< total number of grid points in 'u' coordinate
    double m_du;             //!< physical 'delta' in 'u' coordinate

#if CH_SPACEDIM == 3
    const bool m_periodic_v; //!< is 'v' periodic?
    PetscInt m_num_global_v; //!< total number of grid points in 'u' coordinate
    double m_dv;             //!< physical 'delta' in 'v' coordinate
#endif

    //! vectors to store and manipulate 'F', u', and 'v'
    //! internally, in interaction with the 'AHInterpolation's
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

    // variables
  private:
    const AHInitialGuessPtr
        m_initial_guess; //!< initial guess for AH (saved so that it can be
                         //!< re-used when atempting to solve again)

    const AHParams &m_params; //!< set of AH parameters

    //!< used to compute jacobian of expansion (numerical differentiation)
    static constexpr double eps = 1e-7;

    // PETSc main object
    DM m_dmda;
    //! Scalable Nonlinear Equations Solvers
    SNES m_snes;

    Vec m_snes_soln;
    Vec m_snes_rhs;
    Mat m_snes_jac;

  private:
    //! set the default stencils of AHDerivData at position {u,v}
    void set_stencils(D_DECL(AHDerivData &out, int u, int v));

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
                          const AHInterpolation &interp_plus,
                          const AHInterpolation &interp_minus);

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

#include "PETScAHSolver.impl.hpp"

#endif /* _PETSCAHSOLVER_HPP_ */
