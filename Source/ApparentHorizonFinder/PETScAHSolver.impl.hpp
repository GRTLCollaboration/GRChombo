/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _PETSCAHSOLVER_HPP_
#error "This file should only be included through PETScAHSolver.hpp"
#endif

#ifndef _PETSCAHSOLVER_IMPL_HPP_
#define _PETSCAHSOLVER_IMPL_HPP_

// petsc libraries
#include "petscsys.h"
#include "petscviewerhdf5.h"

// for interpolation at restart:
#include "Lagrange.hpp"
#include "SimpleArrayBox.hpp"
#include "SimpleInterpSource.hpp"

template <class SurfaceGeometry, class AHFunction>
PETScAHSolver<SurfaceGeometry, AHFunction>::PETScAHSolver(
    const AHInterpolation &a_interp, const AHInitialGuessPtr &a_initial_guess,
    const AHParams &a_params)
    : m_interp(a_interp), m_interp_plus(a_interp), m_interp_minus(a_interp),

      m_periodic_u(a_interp.get_coord_system().is_u_periodic()),
      m_num_global_u(a_params.num_points_u),
#if CH_SPACEDIM == 3
      m_periodic_v(a_interp.get_coord_system().is_v_periodic()),
      m_num_global_v(a_params.num_points_v),
#endif

      m_initial_guess(a_initial_guess), m_params(a_params)
{
    initialise();
}

template <class SurfaceGeometry, class AHFunction>
PETScAHSolver<SurfaceGeometry, AHFunction>::~PETScAHSolver()
{
    finalise();
}

template <class SurfaceGeometry, class AHFunction>
void PETScAHSolver<SurfaceGeometry, AHFunction>::initialise()
{
    CH_TIME("PETScAHSolver::initialise_PETSc");

    if (!PETScCommunicator::is_rank_active())
        return;

#if PETSC_VERSION_LT(3, 5, 0)
#define DM_BOUNDARY_PERIODIC DMDA_BOUNDARY_PERIODIC
#define DM_BOUNDARY_GHOSTED DMDA_BOUNDARY_PERIODIC
#endif

#if CH_SPACEDIM == 3
    DMDACreate2d(PETSC_COMM_WORLD,
                 m_periodic_u ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                 m_periodic_v ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                 DMDA_STENCIL_BOX, m_num_global_u,
                 m_num_global_v, /* grid size (negative means that we can
                                    override it using command line) */
                 PETSC_DECIDE, PETSC_DECIDE, /* distribution */
                 1,                          /* number of degrees of freedom */
                 3, /* stencil width (each side from central point) */
                 NULL, NULL, &m_dmda);
#elif CH_SPACEDIM == 2
    DMDACreate1d(PETSC_COMM_WORLD,
                 m_periodic_u ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                 m_num_global_u, /* grid size (negative means that we can
                                    override it using command line) */
                 1,              /* number of degrees of freedom */
                 3, /* stencil width (each side from central point) */
                 NULL, &m_dmda);
#endif

    DMSetUp(m_dmda);

    DMDAGetInfo(m_dmda, NULL, &m_num_global_u,
                D_SELECT(, NULL, &m_num_global_v), NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL);

    m_du = m_interp.get_coord_system().du(m_num_global_u);

#if CH_SPACEDIM == 3
    m_dv = m_interp.get_coord_system().dv(m_num_global_v);
#endif

    DMDAGetCorners(m_dmda, &m_umin, D_SELECT(, NULL, &m_vmin), NULL, &m_nu,
                   D_SELECT(, NULL, &m_nv), NULL);
#if CH_SPACEDIM == 3
    m_vmax = m_vmin + m_nv;
#endif
    m_umax = m_umin + m_nu;

#if CH_SPACEDIM == 3
    int vec_size = m_nu * m_nv;
#elif CH_SPACEDIM == 2
    int vec_size = m_nu;
#endif

    m_F.resize(vec_size);
    m_u.resize(vec_size);
#if CH_SPACEDIM == 3
    m_v.resize(vec_size);
#endif

#if CH_SPACEDIM == 3
    for (int v = m_vmin; v < m_vmax; ++v)
#endif
    {
        for (int u = m_umin; u < m_umax; ++u)
        {
#if CH_SPACEDIM == 3
            m_u[(v - m_vmin) * m_nu + (u - m_umin)] =
                m_interp.get_coord_system().u(u, m_num_global_u);
            m_v[(v - m_vmin) * m_nu + (u - m_umin)] =
                m_interp.get_coord_system().v(v, m_num_global_v);
#elif CH_SPACEDIM == 2
            m_u[u - m_umin] = m_interp.get_coord_system().u(u, m_num_global_u);
#endif
        }
    }

    DMCreateGlobalVector(m_dmda, &m_snes_soln);
    PetscObjectSetName((PetscObject)m_snes_soln, "F");

    VecDuplicate(m_snes_soln, &m_snes_rhs);
    PetscObjectSetName((PetscObject)m_snes_rhs, "expansion");

#if PETSC_VERSION_GE(3, 5, 0)
    DMSetMatType(m_dmda, MATAIJ);
    DMCreateMatrix(m_dmda, &m_snes_jac);
    MatSetOption(m_snes_jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
#else
    DMCreateMatrix(m_dmda, MATAIJ, &m_snes_jac);
    MatSetOption(m_snes_jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
#endif

    SNESCreate(PETSC_COMM_WORLD, &m_snes);
    SNESSetFunction(m_snes, m_snes_rhs, &Petsc_form_function, this);
    SNESSetJacobian(m_snes, m_snes_jac, m_snes_jac, &Petsc_form_jacobian, this);
    SNESMonitorSet(m_snes, &Petsc_SNES_monitor, this, NULL);

    SNESSetFromOptions(m_snes);

    // setup options
    {
        const PETSc_params &pars = m_params.petsc_params;

        SNESType snes_type;
        SNESGetType(m_snes, &snes_type);
        PetscReal snes_atol, snes_rtol, snes_stol;
        PetscInt snes_maxit, snes_maxf;
        SNESGetTolerances(m_snes, &snes_atol, &snes_rtol, &snes_stol,
                          &snes_maxit, &snes_maxf);

        PetscReal snes_divtol;
        SNESGetDivergenceTolerance(m_snes, &snes_divtol);

        KSP snes_ksp;
        SNESGetKSP(m_snes, &snes_ksp);
        KSPType ksp_type;
        KSPGetType(snes_ksp, &ksp_type);
        PetscReal ksp_rtol, ksp_abstol, ksp_dtol;
        PetscInt ksp_maxits;
        KSPGetTolerances(snes_ksp, &ksp_rtol, &ksp_abstol, &ksp_dtol,
                         &ksp_maxits);

        bool any_changed = false; // SNES params
        if (pars.snes_atol > 0 && pars.snes_atol != snes_atol)
        {
            any_changed = true;
            snes_atol = pars.snes_atol;
        }
        if (pars.snes_rtol > 0 && pars.snes_rtol != snes_rtol)
        {
            any_changed = true;
            snes_rtol = pars.snes_rtol;
        }
        if (pars.snes_stol > 0 && pars.snes_stol != snes_stol)
        {
            any_changed = true;
            snes_stol = pars.snes_stol;
        }
        if (pars.snes_maxit > 0 && pars.snes_maxit != snes_maxit)
        {
            any_changed = true;
            snes_maxit = pars.snes_maxit;
        }
        if (pars.snes_maxf > 0 && pars.snes_maxf != snes_maxf)
        {
            any_changed = true;
            snes_maxf = pars.snes_maxf;
        }
        if (any_changed)
        {
            SNESSetTolerances(m_snes, snes_atol, snes_rtol, snes_stol,
                              snes_maxit, snes_maxf);
        }

        if (pars.snes_divtol > 0 && pars.snes_divtol != snes_divtol)
        {
            snes_divtol = pars.snes_divtol;
            SNESSetDivergenceTolerance(m_snes, snes_divtol);
        }

        any_changed = false; // KSP params
        if (pars.ksp_rtol > 0 && pars.ksp_rtol != ksp_rtol)
        {
            any_changed = true;
            ksp_rtol = pars.ksp_rtol;
        }
        if (pars.ksp_abstol > 0 && pars.ksp_abstol != ksp_abstol)
        {
            any_changed = true;
            ksp_abstol = pars.ksp_abstol;
        }
        if (pars.ksp_dtol > 0 && pars.ksp_dtol != ksp_dtol)
        {
            any_changed = true;
            ksp_dtol = pars.ksp_dtol;
        }
        if (pars.ksp_maxits > 0 && pars.ksp_maxits != ksp_maxits)
        {
            any_changed = true;
            ksp_maxits = pars.ksp_maxits;
        }
        if (any_changed)
        {
            KSPSetTolerances(snes_ksp, ksp_rtol, ksp_abstol, ksp_dtol,
                             ksp_maxits);
        }

        if (m_params.verbose > AHParams::MIN)
        {
            pout() << "-------------------------------------\n";
            pout() << "PETScAHSolver Options:\n";
            pout() << "-------------------------------------\n";
            pout() << "PETSc SNES Options:\n";
            pout() << "Type: " << snes_type << "\n";
            pout() << "atol = " << snes_atol << ", rtol = " << snes_rtol
                   << ", stol = " << snes_stol << ",\n";
            pout() << "maxit = " << snes_maxit << ", maxf = " << snes_maxf
                   << "\n";
            pout() << "divtol = " << snes_divtol << "\n";
            pout() << "-------------------------------------\n";
            pout() << "PETSc KSP Options:\n";
            pout() << "Type: " << ksp_type << "\n";
            pout() << "rtol = " << ksp_rtol << ", abstol = " << ksp_abstol
                   << ", dtol = " << ksp_dtol << ", maxits = " << ksp_maxits
                   << "\n";
            pout() << "-------------------------------------" << std::endl;
        }
    }
}

template <class SurfaceGeometry, class AHFunction>
void PETScAHSolver<SurfaceGeometry, AHFunction>::finalise()
{
    if (!PETScCommunicator::is_rank_active())
        return;

    SNESDestroy(&m_snes);
    VecDestroy(&m_snes_soln);
    VecDestroy(&m_snes_rhs);
    MatDestroy(&m_snes_jac);
    DMDestroy(&m_dmda);
}

template <class SurfaceGeometry, class AHFunction>
bool PETScAHSolver<SurfaceGeometry, AHFunction>::interpolate_ah(
    const std::vector<std::vector<double>> &old_coords)
{
    // read PETSc array to 'f'
    dmda_arr_t f;
    DMDAVecGetArray(m_dmda, m_snes_soln, &f);

    // check if number of points changed
    bool points_changed = false;
    const double du = old_coords[0][1] - old_coords[0][0];

#if CH_SPACEDIM == 3
    // assumes ordering in file is done first by fixed 'v', varying 'u'
    // such that number of equal 'v's will be the number of 'u' points
    int old_number_of_u = 1;
    while (abs(old_coords[1][old_number_of_u - 1] -
               old_coords[1][old_number_of_u]) < 1e-7)
        ++old_number_of_u;
    int old_number_of_v = old_coords[0].size() / old_number_of_u;
    const double dv =
        old_coords[1][old_number_of_u] - old_coords[1][old_number_of_u - 1];
    // int total_new_points = m_num_global_u * m_num_global_v;

    points_changed |= (old_number_of_v != m_num_global_v);
#elif CH_SPACEDIM == 2
    int old_number_of_u = old_coords[0].size();
    // int total_new_points = m_num_global_u;
#endif
    points_changed |= (old_number_of_u != m_num_global_u);

    bool force_restart = false;

    if (!points_changed)
    {
        // these should be equal if all is good
        // remember that only PETSc ranks have m_du and m_dv defined
        // but only those should enter here, so all good
#if CH_SPACEDIM == 3
        CH_assert(abs(m_dv - dv) < 1e-7);
        CH_assert(old_coords[0].size() == m_num_global_u * m_num_global_v);
#elif CH_SPACEDIM == 2
        CH_assert(old_coords[0].size() == m_num_global_u);
#endif
        CH_assert(abs(m_du - du) < 1e-7);

        // select that F's that matter to this ranking
#if CH_SPACEDIM == 3
        for (int v = m_vmin; v < m_vmax; ++v)
#endif
        {
            for (int u = m_umin; u < m_umax; ++u)
            {
#if CH_SPACEDIM == 3
                int idx_global = v * m_num_global_u + u;
                double val_v = m_interp.get_coord_system().v(v, m_num_global_v);
                int idx_local = (v - m_vmin) * m_nu + (u - m_umin);
                if (std::abs(val_v - old_coords[1][idx_global]) > 1.e-7)
                    MayDay::Error(
                        "'v' coordinates incompatible with restart data");
#elif CH_SPACEDIM == 2
                int idx_global = u;
                int idx_local = (u - m_umin);
#endif

                double val_u = m_interp.get_coord_system().u(u, m_num_global_u);
                if (std::abs(val_u - old_coords[0][idx_global]) > 1.e-7)
                    MayDay::Error(
                        "'u' coordinates incompatible with restart data");

                m_F[idx_local] = old_coords[CH_SPACEDIM - 1][idx_global];

#if CH_SPACEDIM == 3
                // as 'write_coords_file' preserves order
                f[v][u] = m_F[idx_local];
#elif CH_SPACEDIM == 2
                f[u] = m_F[idx_local];
#endif
            }
        }
    }
    else
    {
        force_restart = true;

        if (m_params.verbose > AHParams::NONE)
        {
            MayDay::Error("The same number of u or v points in the AH finder "
                          "has to be used on the restart!");
        }

        // start interpolating

#if CH_SPACEDIM == 3
        // local, no derivative
        std::array<int, CH_SPACEDIM - 1> derivs = {0, 0};
        std::array<double, CH_SPACEDIM - 1> dxs = {du, dv};

        SimpleArrayBox<CH_SPACEDIM - 1> box(
            {old_number_of_u, old_number_of_v}, old_coords[CH_SPACEDIM - 1],
            {m_interp.get_coord_system().is_u_periodic(),
             m_interp.get_coord_system().is_v_periodic()});
        SimpleInterpSource<CH_SPACEDIM - 1> source(
            {old_number_of_u, old_number_of_v}, dxs,
            {m_interp.get_coord_system().is_u_periodic(),
             m_interp.get_coord_system().is_v_periodic()});
#elif CH_SPACEDIM == 2
        std::array<int, CH_SPACEDIM - 1> derivs = {0};
        std::array<double, CH_SPACEDIM - 1> dxs = {du};

        SimpleArrayBox<CH_SPACEDIM - 1> box(
            {old_number_of_u}, old_coords[CH_SPACEDIM - 1],
            {m_interp.get_coord_system().is_u_periodic()});
        SimpleInterpSource<CH_SPACEDIM - 1> source(
            {old_number_of_u}, dxs,
            {m_interp.get_coord_system().is_u_periodic()});
#endif

        const int Order = 4;
        bool verbose = false;
        Lagrange<Order, CH_SPACEDIM - 1> interpolator4(source, verbose);

#if CH_SPACEDIM == 3
        for (int v = m_vmin; v < m_vmax; ++v)
        {
            double v_val = m_interp.get_coord_system().v(v, m_num_global_v);
            double v_old_idx = v_val / dv;
            // round such that we don't get AMRInterpolator errors
            // (e.g. the last point idx=N-1 should be exact, and not being exact
            // was generating this error in the Lagrange intrepolator trying to
            // access cells beyond the last cell, and hence error... this fixed
            // it)
            if (v_old_idx - int(v_old_idx) < 1.e-7)
                v_old_idx = int(v_old_idx);
#else
        {
#endif
            for (int u = m_umin; u < m_umax; ++u)
            {
                double u_val = m_interp.get_coord_system().u(u, m_num_global_u);
                double u_old_idx = u_val / du;
                // round such that we don't get AMRInterpolator errors (same as
                // for 'v')
                if (u_old_idx - int(u_old_idx) < 1.e-7)
                    u_old_idx = int(u_old_idx);

#if CH_SPACEDIM == 3
                // int idx_global = v * m_num_global_u + u;
                int idx_local = (v - m_vmin) * m_nu + (u - m_umin);

                std::array<double, CH_SPACEDIM - 1> evalCoord = {u_old_idx,
                                                                 v_old_idx};
#elif CH_SPACEDIM == 2
                // int idx_global = u;
                int idx_local = (u - m_umin);
                std::array<double, CH_SPACEDIM - 1> evalCoord = {u_old_idx};
#endif

                interpolator4.setup(derivs, evalCoord);
                m_F[idx_local] = interpolator4.interpData(box);

#if CH_SPACEDIM == 3
                // as 'write_coords_file' preserves order
                f[v][u] = m_F[idx_local];
#elif CH_SPACEDIM == 2
                f[u] = m_F[idx_local];
#endif
            }
        }

        if (m_params.verbose > AHParams::NONE)
        {
            pout() << "Interpolation Successfull." << std::endl;
        }
    }

    // write PETSc array back
    DMDAVecRestoreArray(m_dmda, m_snes_soln, &f);

    return force_restart;
}

template <class SurfaceGeometry, class AHFunction>
const std::array<double, CH_SPACEDIM> &
PETScAHSolver<SurfaceGeometry, AHFunction>::get_origin() const
{
    return m_interp.get_coord_system().get_origin();
}

template <class SurfaceGeometry, class AHFunction>
void PETScAHSolver<SurfaceGeometry, AHFunction>::set_origin(
    const std::array<double, CH_SPACEDIM> &a_origin)
{
    m_interp.set_origin(a_origin);
    m_interp_plus.set_origin(a_origin);
    m_interp_minus.set_origin(a_origin);

    if (m_params.verbose > AHParams::SOME)
    {
        pout() << "Setting origin to (" << a_origin[0] << "," << a_origin[1]
#if CH_SPACEDIM == 3
               << "," << a_origin[2]
#endif
               << ")" << std::endl;
    }
}

template <class SurfaceGeometry, class AHFunction>
const AHInitialGuessPtr &
PETScAHSolver<SurfaceGeometry, AHFunction>::get_initial_guess() const
{
    return m_initial_guess;
}

template <class SurfaceGeometry, class AHFunction>
void PETScAHSolver<SurfaceGeometry, AHFunction>::reset_initial_guess()
{
    CH_TIME("ApparentHorizon::reset_initial_guess");

    if (!PETScCommunicator::is_rank_active())
        return;

    auto origin = get_origin();

    // read PETSc array to 'f'
    dmda_arr_t f;
    DMDAVecGetArray(m_dmda, m_snes_soln, &f);

    bool out_of_grid = false;
    const SurfaceGeometry &coord_system = m_interp.get_coord_system();

#if CH_SPACEDIM == 3
    for (int v = m_vmin; v < m_vmax; ++v)
#endif
    {
        for (int u = m_umin; u < m_umax; ++u)
        {
            double u_val = m_interp.get_coord_system().u(u, m_num_global_u);
#if CH_SPACEDIM == 3
            double v_val = m_interp.get_coord_system().v(v, m_num_global_v);
            double &f_point = f[v][u];
            f_point = m_initial_guess->get(u_val, v_val);
            double z = coord_system.get_grid_coord(2, f_point, u_val, v_val);
#elif CH_SPACEDIM == 2
            double &f_point = f[u];
            f_point = m_initial_guess->get(u_val);
#endif
            double x =
                coord_system.get_grid_coord(0, D_DECL(f_point, u_val, v_val));
            double y =
                coord_system.get_grid_coord(1, D_DECL(f_point, u_val, v_val));
            out_of_grid |= m_interp.is_in_grid(D_DECL(x, y, z));
        }
    }

    CH_assert(!out_of_grid);

    // write PETSc array back
    DMDAVecRestoreArray(m_dmda, m_snes_soln, &f);
}

template <class SurfaceGeometry, class AHFunction>
void PETScAHSolver<SurfaceGeometry, AHFunction>::solve()
{
    // actual solve happens here!
    SNESSolve(m_snes, NULL, m_snes_soln);

    PetscInt its;
    SNESGetIterationNumber(m_snes, &its);
    if (m_params.verbose > AHParams::MIN)
    {
        pout() << "SNES Iteration number " << its << endl;
    }
    SNESGetLinearSolveIterations(m_snes, &its);
    if (m_params.verbose > AHParams::MIN)
    {
        pout() << "KSP Iteration number " << its << endl;
    }
}

template <class SurfaceGeometry, class AHFunction>
SNESConvergedReason
PETScAHSolver<SurfaceGeometry, AHFunction>::getConvergedReason() const
{
    SNESConvergedReason reason;
    SNESGetConvergedReason(m_snes, &reason);
    return reason;
}

template <class SurfaceGeometry, class AHFunction>
void PETScAHSolver<SurfaceGeometry, AHFunction>::get_dmda_arr_t(Vec &localF,
                                                                dmda_arr_t &in)
{
    DMGetLocalVector(m_dmda, &localF);
    DMGlobalToLocalBegin(m_dmda, m_snes_soln, INSERT_VALUES, localF);
    DMGlobalToLocalEnd(m_dmda, m_snes_soln, INSERT_VALUES, localF);
    DMDAVecGetArray(m_dmda, localF, &in);
}
template <class SurfaceGeometry, class AHFunction>
void PETScAHSolver<SurfaceGeometry, AHFunction>::restore_dmda_arr_t(
    Vec &localF, dmda_arr_t &in)
{
    DMDAVecRestoreArray(m_dmda, localF, &in);
    DMRestoreLocalVector(m_dmda, &localF);
}

template <class SurfaceGeometry, class AHFunction>
void PETScAHSolver<SurfaceGeometry, AHFunction>::set_stencils(AHDerivData &out,
                                                              int u
#if CH_SPACEDIM == 3
                                                              ,
                                                              int v
#endif
)
{
#if CH_SPACEDIM == 3
    const bool periodic[CH_SPACEDIM - 1] = {m_periodic_u, m_periodic_v};

    int coord[CH_SPACEDIM - 1] = {u, v};
    const int N[CH_SPACEDIM - 1] = {m_num_global_u, m_num_global_v};

    int d_start[CH_SPACEDIM - 1];
    double *const d_weights[CH_SPACEDIM - 1] = {out.du_weights, out.dv_weights};

    int dd_start[CH_SPACEDIM - 1];
    double *const dd_weights[CH_SPACEDIM - 1] = {out.dudu_weights,
                                                 out.dvdv_weights};
#elif CH_SPACEDIM == 2
    const bool periodic[CH_SPACEDIM - 1] = {m_periodic_u};

    int coord[CH_SPACEDIM - 1] = {u};
    const int N[CH_SPACEDIM - 1] = {m_num_global_u};

    int d_start[CH_SPACEDIM - 1];
    double *const d_weights[CH_SPACEDIM - 1] = {out.du_weights};

    int dd_start[CH_SPACEDIM - 1];
    double *const dd_weights[CH_SPACEDIM - 1] = {out.dudu_weights};
#endif

    for (int i = 0; i < CH_SPACEDIM - 1; ++i)
    {

        // Non-periodic boundaries
        if (!periodic[i] && (coord[i] == 0 || coord[i] == 1 ||
                             coord[i] == N[i] - 2 || coord[i] == N[i] - 1))
        {

            if (coord[i] == 0)
            {

                d_start[i] = 0;

                d_weights[i][0] = -2.083333333333333333333;
                d_weights[i][1] = +4.000000000000000000000;
                d_weights[i][2] = -3.000000000000000000000;
                d_weights[i][3] = +1.333333333333333333333;
                d_weights[i][4] = -0.250000000000000000000;

                dd_start[i] = 0;

                dd_weights[i][0] = +3.75000000000000000000;
                dd_weights[i][1] = -12.8333333333333333333;
                dd_weights[i][2] = +17.8333333333333333333;
                dd_weights[i][3] = -13.0000000000000000000;
                dd_weights[i][4] = +5.08333333333333333333;
                dd_weights[i][5] = -0.83333333333333333333;
            }
            else if (coord[i] == 1)
            {

                d_start[i] = -1;

                d_weights[i][0] = -0.250000000000000000000;
                d_weights[i][1] = -0.833333333333333333333;
                d_weights[i][2] = +1.500000000000000000000;
                d_weights[i][3] = -0.500000000000000000000;
                d_weights[i][4] = +0.083333333333333333333;

                dd_start[i] = -1;

                dd_weights[i][0] = +0.83333333333333333333;
                dd_weights[i][1] = -1.25000000000000000000;
                dd_weights[i][2] = -0.33333333333333333333;
                dd_weights[i][3] = +1.16666666666666666667;
                dd_weights[i][4] = -0.50000000000000000000;
                dd_weights[i][5] = +0.08333333333333333333;
            }
            else if (coord[i] == N[i] - 2)
            {

                d_start[i] = -3;

                d_weights[i][0] = -0.083333333333333333333;
                d_weights[i][1] = +0.500000000000000000000;
                d_weights[i][2] = -1.500000000000000000000;
                d_weights[i][3] = +0.833333333333333333333;
                d_weights[i][4] = +0.250000000000000000000;

                dd_start[i] = -4;

                dd_weights[i][0] = +0.08333333333333333333;
                dd_weights[i][1] = -0.50000000000000000000;
                dd_weights[i][2] = +1.16666666666666666667;
                dd_weights[i][3] = -0.33333333333333333333;
                dd_weights[i][4] = -1.25000000000000000000;
                dd_weights[i][5] = +0.83333333333333333333;
            }
            else if (coord[i] == N[i] - 1)
            {

                d_start[i] = -4;

                d_weights[i][0] = +0.250000000000000000000;
                d_weights[i][1] = -1.333333333333333333333;
                d_weights[i][2] = +3.000000000000000000000;
                d_weights[i][3] = -4.000000000000000000000;
                d_weights[i][4] = +2.083333333333333333333;

                dd_start[i] = -5;

                dd_weights[i][0] = -0.83333333333333333333;
                dd_weights[i][1] = +5.08333333333333333333;
                dd_weights[i][2] = -13.0000000000000000000;
                dd_weights[i][3] = +17.8333333333333333333;
                dd_weights[i][4] = -12.8333333333333333333;
                dd_weights[i][5] = +3.75000000000000000000;
            }
        }
        // Standard points
        else
        {

            d_start[i] = -2;

            d_weights[i][0] = +0.083333333333333333333;
            d_weights[i][1] = -0.666666666666666666666;
            d_weights[i][2] = 0;
            d_weights[i][3] = +0.666666666666666666666;
            d_weights[i][4] = -0.083333333333333333333;

            if (coord[i] == N[i] - 3)
            {
                dd_start[i] = -3;

                dd_weights[i][0] = 0;
                dd_weights[i][1] = -0.08333333333333333333;
                dd_weights[i][2] = +1.33333333333333333333;
                dd_weights[i][3] = -2.50000000000000000000;
                dd_weights[i][4] = +1.33333333333333333333;
                dd_weights[i][5] = -0.08333333333333333333;
            }
            else
            {
                dd_start[i] = -2;

                dd_weights[i][0] = -0.08333333333333333333;
                dd_weights[i][1] = +1.33333333333333333333;
                dd_weights[i][2] = -2.50000000000000000000;
                dd_weights[i][3] = +1.33333333333333333333;
                dd_weights[i][4] = -0.08333333333333333333;
                dd_weights[i][5] = 0;
            }
        }
    }

    out.du_stencil_start = d_start[0];
    out.dudu_stencil_start = dd_start[0];

#if CH_SPACEDIM == 3
    out.dv_stencil_start = d_start[1];
    out.dvdv_stencil_start = dd_start[1];
#endif
}

template <class SurfaceGeometry, class AHFunction>
AHDerivData PETScAHSolver<SurfaceGeometry, AHFunction>::diff(
    D_DECL(const dmda_arr_t in, int u, int v))
{
    CH_TIME("PETScAHSolver::diff");

    AHDerivData out;
    set_stencils(out, u
#if CH_SPACEDIM == 3
                 ,
                 v
#endif
    );

#if CH_SPACEDIM == 3
    // d/du
    for (int j = 0; j < DWIDTH; ++j)
    {
        out.duF +=
            out.du_weights[j] * in[v][u + out.du_stencil_start + j] / m_du;

        // d2/dudv
        for (int k = 0; k < DWIDTH; ++k)
        {
            out.dudvF +=
                out.du_weights[j] * out.dv_weights[k] *
                in[v + out.dv_stencil_start + k][u + out.du_stencil_start + j] /
                (m_du * m_dv);
        }
    }

    // d/dv
    for (int j = 0; j < DWIDTH; ++j)
    {
        out.dvF +=
            out.dv_weights[j] * in[v + out.dv_stencil_start + j][u] / m_dv;
    }

    // d2/du2
    for (int j = 0; j < DDWIDTH; ++j)
    {
        out.duduF += out.dudu_weights[j] *
                     in[v][u + out.dudu_stencil_start + j] / (m_du * m_du);
    }

    // d2/dv2
    for (int j = 0; j < DDWIDTH; ++j)
    {
        out.dvdvF += out.dvdv_weights[j] *
                     in[v + out.dvdv_stencil_start + j][u] / (m_dv * m_dv);
    }

#elif CH_SPACEDIM == 2

    // d/du
    for (int j = 0; j < DWIDTH; ++j)
    {
        out.duF += out.du_weights[j] * in[u + out.du_stencil_start + j] / m_du;
    }

    // d2/du2
    for (int j = 0; j < DDWIDTH; ++j)
    {
        out.duduF += out.dudu_weights[j] * in[u + out.dudu_stencil_start + j] /
                     (m_du * m_du);
    }

#endif

    return out;
}

template <class SurfaceGeometry, class AHFunction>
void PETScAHSolver<SurfaceGeometry, AHFunction>::form_function(Vec F, Vec Rhs)
{
    CH_TIME("PETScAHSolver::form_function");

    // Scatter ghost cells
    Vec localF;
    DMGetLocalVector(m_dmda, &localF);
    DMGlobalToLocalBegin(m_dmda, F, INSERT_VALUES, localF);
    DMGlobalToLocalEnd(m_dmda, F, INSERT_VALUES, localF);

    dmda_arr_t in;
    DMDAVecGetArray(m_dmda, localF, &in);

    dmda_arr_t out;
    DMDAVecGetArray(m_dmda, Rhs, &out);

#if CH_SPACEDIM == 3
    for (int v = m_vmin; v < m_vmax; ++v)
#endif
    {
        for (int u = m_umin; u < m_umax; ++u)
        {
#if CH_SPACEDIM == 3
            m_F[(v - m_vmin) * m_nu + (u - m_umin)] = in[v][u];
#elif CH_SPACEDIM == 2
            m_F[u - m_umin] = in[u];
#endif
        }
    }

    bool out_of_grid = m_interp.set_coordinates(D_DECL(m_F, m_u, m_v));
    // abort if out of grid - reduces the time in divergence dramatically!
    if (out_of_grid)
        SNESSetFunctionDomainError(m_snes);

    // must interpolate anyway, as other PETSc ranks do not know if this one is
    // "out_of_grid"
    m_interp.interpolate();

    int idx = 0;

#if CH_SPACEDIM == 3
    for (int v = m_vmin; v < m_vmax; ++v)
#endif
    {
        for (int u = m_umin; u < m_umax; ++u)
        {

#if CH_SPACEDIM == 3
            double &_out = out[v][u];
#elif CH_SPACEDIM == 2
            double &_out = out[u];
#endif

            AHDerivData deriv = diff(D_DECL(in, u, v));

            if (out_of_grid) // no need to calculate expansion
                _out = 0.;
            else if (!m_periodic_u && (u == 0 || u == m_num_global_u - 1))
                _out = deriv.duF;

#if CH_SPACEDIM == 3
            else if (!m_periodic_v && (v == 0 || v == m_num_global_v - 1))
                _out = deriv.dvF;
#endif

            else
            {
                const auto geometry_data = m_interp.get_geometry_data(idx);
                const auto data = m_interp.get_data(idx);
                const auto coords = m_interp.get_coords(idx);
                const auto coords_cart = m_interp.get_cartesian_coords(idx);
                AHFunction func(data, coords, coords_cart);
                _out = func.get(geometry_data, deriv, m_params.func_params);
            }

            ++idx;
        }
    }

    DMDAVecRestoreArray(m_dmda, localF, &in);
    DMDAVecRestoreArray(m_dmda, Rhs, &out);
    DMRestoreLocalVector(m_dmda, &localF);
}

template <class SurfaceGeometry, class AHFunction>
void PETScAHSolver<SurfaceGeometry, AHFunction>::form_jacobian(Vec F, Mat J)
{
    CH_TIME("PETScAHSolver::form_jacobian");

    // Scatter ghost cells
    Vec localF;
    DMGetLocalVector(m_dmda, &localF);
    DMGlobalToLocalBegin(m_dmda, F, INSERT_VALUES, localF);
    DMGlobalToLocalEnd(m_dmda, F, INSERT_VALUES, localF);

    dmda_arr_t in;
    DMDAVecGetArray(m_dmda, localF, &in);

    bool out_of_grid =
        m_interp_plus.set_coordinates(D_DECL(m_F, m_u, m_v), eps);
    out_of_grid |= m_interp_minus.set_coordinates(D_DECL(m_F, m_u, m_v), -eps);
    // abort if out of grid - reduces the time in divergence dramatically!
    if (out_of_grid)
        SNESSetFunctionDomainError(m_snes);

    // must interpolate anyway, as other PETSc ranks do not know if this one is
    // "out_of_grid"
    m_interp_plus.interpolate();
    m_interp_minus.interpolate();

    int idx = 0;

#if CH_SPACEDIM == 3
    for (int v = m_vmin; v < m_vmax; ++v)
#endif
    {
        for (int u = m_umin; u < m_umax; ++u)
        {
            MatStencil row[1] = {0};
            row[0].i = u;
#if CH_SPACEDIM == 3
            row[0].j = v;
#endif

            if (out_of_grid) // no need to calculate expansion
            {
                MatStencil col[DWIDTH] = {0};
                double val[DWIDTH] = {0};
                MatSetValuesStencil(J, 1, row, DWIDTH, col, val, INSERT_VALUES);
            }
            else if (!m_periodic_u && (u == 0 || u == m_num_global_u - 1))
            {
                MatStencil col[DWIDTH] = {0};
                double val[DWIDTH] = {0};

                const AHDerivData deriv = diff(D_DECL(in, u, v));

                for (int a = 0; a < DWIDTH; ++a)
                {
                    col[a].i = u + deriv.du_stencil_start + a;
#if CH_SPACEDIM == 3
                    col[a].j = v;
#endif
                    val[a] = deriv.du_weights[a];
                }

                MatSetValuesStencil(J, 1, row, DWIDTH, col, val, INSERT_VALUES);
            }

#if CH_SPACEDIM == 3
            else if (!m_periodic_v && (v == 0 || v == m_num_global_v - 1))
            {
                MatStencil col[DWIDTH] = {0};
                double val[DWIDTH] = {0};

                const AHDerivData deriv = diff(in, u, v);

                for (int b = 0; b < DWIDTH; ++b)
                {
                    col[b].i = u;
                    col[b].j = v + deriv.dv_stencil_start + b;
                    val[b] = deriv.dv_weights[b];
                }

                MatSetValuesStencil(J, 1, row, DWIDTH, col, val, INSERT_VALUES);
            }
#endif
            else
            {

#if CH_SPACEDIM == 3
                const int NVAL = DDWIDTH * DDWIDTH;
#elif CH_SPACEDIM == 2
                const int NVAL = DDWIDTH;
#endif
                MatStencil col[NVAL] = {0};
                double val[NVAL] = {0};

                const AHDerivData deriv_default = diff(D_DECL(in, u, v));

                // "local" and "stencil" jacobian terms
                // "local" (point {u,v}) corresponds to (a ==
                // -deriv_default.dudu_stencil_start
                // && b == -deriv_default.dvdv_stencil_start) and has to use
                // 'm_interp_plus' and 'm_interp_minus'
                {
#if CH_SPACEDIM == 3
                    for (int b = 0; b < DDWIDTH; ++b)
                    {
#elif CH_SPACEDIM == 2
                    {
                        int b = 0;
#endif
                        for (int a = 0; a < DDWIDTH; ++a)
                        {
                            col[b * DDWIDTH + a].i =
                                u + deriv_default.dudu_stencil_start + a;
                            const int uu = col[b * DDWIDTH + a].i;
                            const int u_start =
                                -deriv_default.dudu_stencil_start;
#if CH_SPACEDIM == 3
                            col[b * DDWIDTH + a].j =
                                v + deriv_default.dvdv_stencil_start + b;
                            const int vv = col[b * DDWIDTH + a].j;
                            const int v_start =
                                -deriv_default.dvdv_stencil_start;
#elif CH_SPACEDIM == 2
                            const int v_start = 0;
#endif

                            if (a == u_start && b == v_start) // <=> uu=u, vv=v
                            {
                                val[b * DDWIDTH + a] = point_jacobian(
                                    u, uu,
#if CH_SPACEDIM == 3
                                    v, vv,
#endif
                                    in, idx, m_interp_plus, m_interp_minus);
                            }
                            else
                            {
                                val[b * DDWIDTH + a] =
                                    point_jacobian(u, uu,
#if CH_SPACEDIM == 3
                                                   v, vv,
#endif
                                                   in, idx, m_interp, m_interp);
                            }
                        }
                    }
                }

                MatSetValuesStencil(J, 1, row, NVAL, col, val, INSERT_VALUES);
            }

            ++idx;
        }
    }

    DMDAVecRestoreArray(m_dmda, localF, &in);
    DMRestoreLocalVector(m_dmda, &localF);

    MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
}

template <class SurfaceGeometry, class AHFunction>
double PETScAHSolver<SurfaceGeometry, AHFunction>::point_jacobian(
    int u, int u_stencil,
#if CH_SPACEDIM == 3
    int v, int v_stencil,
#endif
    dmda_arr_t in, int idx, const AHInterpolation &interp_plus,
    const AHInterpolation &interp_minus)
{
#if CH_SPACEDIM == 3
    double &_in = in[v_stencil][u_stencil];
#elif CH_SPACEDIM == 2
    double &_in = in[u_stencil];
#endif

    // "plus" perturbation
    double expansionPlus;
    {
        double in_old = _in;
        _in += eps; // perturb just the point {u,v}

        const AHDerivData deriv = diff(D_DECL(in, u, v));

        const auto geometry_data = interp_plus.get_geometry_data(idx);
        const auto data = interp_plus.get_data(idx);
        const auto coords = interp_plus.get_coords(idx);
        const auto coords_cart = interp_plus.get_cartesian_coords(idx);
        AHFunction func(data, coords, coords_cart);
        expansionPlus = func.get(geometry_data, deriv, m_params.func_params);

        _in = in_old;
    }

    // "minus" perturbation
    double expansionMinus;
    {
        double in_old = _in;
        _in -= eps; // perturb just the point {u,v}

        const AHDerivData deriv = diff(D_DECL(in, u, v));

        const auto geometry_data = interp_minus.get_geometry_data(idx);
        const auto data = interp_minus.get_data(idx);
        const auto coords = interp_minus.get_coords(idx);
        const auto coords_cart = interp_minus.get_cartesian_coords(idx);
        AHFunction func(data, coords, coords_cart);
        expansionMinus = func.get(geometry_data, deriv, m_params.func_params);

        _in = in_old;
    }

    return (expansionPlus - expansionMinus) / (2. * eps);
}

//! functions used by PETSc based on 'form_function' and 'form_jacobian'
template <class SurfaceGeometry, class AHFunction>
PetscErrorCode PETScAHSolver<SurfaceGeometry, AHFunction>::Petsc_form_function(
    SNES snes, Vec F, Vec Rhs, void *ptr)
{
    PETScAHSolver &ah = *reinterpret_cast<PETScAHSolver *>(ptr);
    CH_assert(ah.m_snes == snes);
    ah.form_function(F, Rhs);
    return 0;
}

template <class SurfaceGeometry, class AHFunction>
PetscErrorCode
#if PETSC_VERSION_GE(3, 5, 0)
PETScAHSolver<SurfaceGeometry, AHFunction>::Petsc_form_jacobian(SNES snes,
                                                                Vec F, Mat Amat,
                                                                Mat Pmat,
                                                                void *ptr)
#else
PETScAHSolver<SurfaceGeometry, AHFunction>::Petsc_form_jacobian(
    SNES snes, Vec F, Mat *Amat, Mat *Pmat, MatStructure *flag, void *ptr)
#endif
{
    PETScAHSolver &ah = *reinterpret_cast<PETScAHSolver *>(ptr);

#if PETSC_VERSION_GE(3, 5, 0)
    CH_assert(ah.m_snes == snes && Amat == Pmat);
    ah.form_jacobian(F, Amat);
#else
    CH_assert(ah.m_snes == snes && *Amat == *Pmat);
    ah.form_jacobian(F, *Amat);
#endif

    return 0;
}

template <class SurfaceGeometry, class AHFunction>
PetscErrorCode PETScAHSolver<SurfaceGeometry, AHFunction>::Petsc_SNES_monitor(
    SNES snes, PetscInt its, PetscReal norm, void *ptr)
{
    PETScAHSolver &ah = *reinterpret_cast<PETScAHSolver *>(ptr);
    CH_assert(ah.m_snes == snes);
    return 0;
}

#endif // _PETSCAHSOLVER_IMPL_HPP_
