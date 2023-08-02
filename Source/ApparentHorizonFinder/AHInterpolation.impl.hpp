/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHINTERPOLATION_HPP_
#error "This file should only be included through AHInterpolation.hpp"
#endif

#ifndef _AHINTERPOLATION_IMPL_HPP_
#define _AHINTERPOLATION_IMPL_HPP_

#include "PETScCommunicator.hpp"

template <class SurfaceGeometry, class AHFunction>
AHInterpolation_t<SurfaceGeometry, AHFunction>::AHInterpolation_t(
    const SurfaceGeometry &a_coord_system,
    AMRInterpolator<Lagrange<4>> *a_interpolator)
    : m_coord_system(a_coord_system), m_interpolator(a_interpolator)
{
    // below:
    // determine maximum and minimum physical coordinate of the grid
    // (so that 'is_in_grid' knows when PETSc has diverged out of the grid and
    // doesn't let him do so,
    //  as it would cause an error in the AMRInterpolator)

    const Box &domainBox = const_cast<AMR &>(m_interpolator->getAMR())
                               .getAMRLevels()[0]
                               ->problemDomain()
                               .domainBox();
    const IntVect &small_end = domainBox.smallEnd();
    const IntVect &big_end = domainBox.bigEnd();

    const std::array<double, CH_SPACEDIM> &coarsest_dx =
        m_interpolator->get_coarsest_dx();
    const std::array<double, CH_SPACEDIM> &coarsest_origin =
        m_interpolator->get_coarsest_origin();

    // set coordinate minimum and maximum as ONE cells before the boundary
    for (unsigned i = 0; i < CH_SPACEDIM; ++i)
    {
        if (m_interpolator->get_boundary_reflective(Side::Lo, i))
        {
            m_coord_min[i] =
                -(big_end[i] - 1) * coarsest_dx[i] - coarsest_origin[i];
        }
        else if (m_interpolator->get_boundary_periodic(i))
        {
            m_coord_min[i] = small_end[i] * coarsest_dx[i];
        }
        else
        {
            m_coord_min[i] =
                (small_end[i] + 1) * coarsest_dx[i] + coarsest_origin[i];
        }

        if (m_interpolator->get_boundary_reflective(Side::Hi, i))
        {
            m_coord_max[i] =
                (2 * big_end[i] - 1) * coarsest_dx[i] + coarsest_origin[i];
        }
        else if (m_interpolator->get_boundary_periodic(i))
        {
            m_coord_max[i] = (big_end[i] + 1) * coarsest_dx[i];
        }
        else
        {
            m_coord_max[i] =
                (big_end[i] - 1) * coarsest_dx[i] + coarsest_origin[i];
        }
        // pout() << "m_coord = " << m_coord_min[i] << "\t" << m_coord_max[i]
        //        << " with coarse origin " << coarsest_origin[i] << " and dx "
        //        << coarsest_dx[i] << std::endl;
    }
}

template <class SurfaceGeometry, class AHFunction>
const AMRInterpolator<Lagrange<4>> *
AHInterpolation_t<SurfaceGeometry, AHFunction>::get_interpolator() const
{
    return m_interpolator;
}

template <class SurfaceGeometry, class AHFunction>
const SurfaceGeometry &
AHInterpolation_t<SurfaceGeometry, AHFunction>::get_coord_system() const
{
    return m_coord_system;
}

template <class SurfaceGeometry, class AHFunction>
std::vector<std::string>
AHInterpolation_t<SurfaceGeometry, AHFunction>::get_labels() const
{
    return
    {
        m_coord_system.u_name(),
#if CH_SPACEDIM == 3
            m_coord_system.v_name(),
#endif
            m_coord_system.param_name()
    };
}

template <class SurfaceGeometry, class AHFunction>
void AHInterpolation_t<SurfaceGeometry, AHFunction>::set_origin(
    const std::array<double, CH_SPACEDIM> &origin)
{
    m_coord_system.set_origin(origin);
}

template <class SurfaceGeometry, class AHFunction>
bool AHInterpolation_t<SurfaceGeometry, AHFunction>::is_in_grid(
    D_DECL(double &x, double &y, double &z))
{
    CH_TIME("AHInterpolation::is_in_grid");

    bool out_of_grid = false;

    // if out of bounds, put back in grid
    if (x < m_coord_min[0])
    {
        out_of_grid = true;
        x = m_coord_min[0];
    }
    if (x > m_coord_max[0])
    {
        out_of_grid = true;
        x = m_coord_max[0];
    }

    if (y < m_coord_min[1])
    {
        out_of_grid = true;
        y = m_coord_min[1];
    }
    if (y > m_coord_max[1])
    {
        out_of_grid = true;
        y = m_coord_max[1];
    }

#if CH_SPACEDIM == 3
    if (z < m_coord_min[2])
    {
        out_of_grid = true;
        z = m_coord_min[2];
    }
    if (z > m_coord_max[2])
    {
        out_of_grid = true;
        z = m_coord_max[2];
    }
#endif

    return out_of_grid;
}

//! triplet of functions to be used together in blocks of code that require
//! PETSc AND AMRInterpolator to interpolate 'keep_interpolating_if_inactive'
//! returns 'true' immediately for PETSc ranks for non-PETSc ranks, it loops in
//! a sequence of 'interpolate()' waiting for PETSc ranks to call
//! 'interpolate()'' as well (as this needs to be called by all Chombo
//! processes) can be aborted by calling 'break_interpolation_loop()' for PETSc
//! ranks
/** Example of usage:
 * if(m_geom.keep_interpolating_if_inactive())
 * {
 *     (... PETSc code that calls 'set_coordinates()'...)
 *     m_geom.break_interpolation_loop();
 * }
 */
template <class SurfaceGeometry, class AHFunction>
bool AHInterpolation_t<SurfaceGeometry,
                       AHFunction>::keep_interpolating_if_inactive()
{
    CH_TIME("AMRInterpolator::keep_interpolating_if_inactive");

    if (!PETScCommunicator::is_rank_active())
    {
        int keep_interpolating = 1;
        while (keep_interpolating)
            keep_interpolating = interpolate();
        return false;
    }
    else
        return true;
}

template <class SurfaceGeometry, class AHFunction>
void AHInterpolation_t<SurfaceGeometry, AHFunction>::break_interpolation_loop()
    const
{
    CH_TIME("AMRInterpolator::break_interpolation_loop");

#ifdef CH_MPI
    // break "keep_interpolating" loop
    int ZERO = 0;
    int keep_interpolating = 0;
    // this is also called in 'interpolate()', breaking the loop
    MPI_Allreduce(&ZERO, &keep_interpolating, 1, MPI_INT, MPI_LAND,
                  Chombo_MPI::comm);
#endif
}
template <class SurfaceGeometry, class AHFunction>
int AHInterpolation_t<SurfaceGeometry, AHFunction>::interpolate()
{
    CH_TIME("AHInterpolation::interpolate");

    // Code below used to allow PETSc to run in a sub-communicator.
    // For that, PETSc processes run 'set_coordinates' and 'interpolate', while
    // non-PETSc processes loop through 'interpolate'. If all cores run
    // interpolate, 'MPI_Allreduce' will return 'keep_interpolating' as true. To
    // exit the loop for non-PETSc cores, simply do on them a call to:
    // MPI_Allreduce(&ZERO, &keep_interpolating, 1, MPI_INT, MPI_LAND,
    // Chombo_MPI::comm); (this is done in 'break_interpolation_loop()')

    int keep_interpolating = 1;
#ifdef CH_MPI
    int ONE = 1;
    MPI_Allreduce(&ONE, &keep_interpolating, 1, MPI_INT, MPI_LAND,
                  Chombo_MPI::comm);
#endif

    if (!keep_interpolating)
        return 0;

    const int n = m_x.size();

    InterpolationQuery query(n);
    query.setCoords(0, m_x.data()).setCoords(1, m_y.data());
#if CH_SPACEDIM == 3
    query.setCoords(2, m_z.data());
#endif

    for (int i = AHFunction::vars_min(); i <= AHFunction::vars_max(); ++i)
        m_data.set_vars(query, i, i, VariableType::evolution, n);
    for (int i = AHFunction::d1_vars_min(); i <= AHFunction::d1_vars_max(); ++i)
        m_data.set_d1(query, i, i, VariableType::evolution, n);
    for (int i = AHFunction::d2_vars_min(); i <= AHFunction::d2_vars_max(); ++i)
        m_data.set_d2(query, i, i, VariableType::evolution, n);

    m_interpolator->interp(query);

    return 1;
}

template <class SurfaceGeometry, class AHFunction>
void AHInterpolation_t<SurfaceGeometry, AHFunction>::refresh_interpolator(
    bool printing_step,
    const std::map<std::string, std::tuple<int, VariableType, int>> &extra_vars)
{
    CH_TIME("AHInterpolation::refresh_interpolator");

    // brute force filling of all ghosts:
    // m_interpolator->refresh();

    // OR try to only interpolate the variables needed:

    int min_evolution_var =
        std::min(AHFunction::vars_min(), std::min(AHFunction::d1_vars_min(),
                                                  AHFunction::d2_vars_min()));
    int max_evolution_var =
        std::max(AHFunction::vars_max(), std::max(AHFunction::d1_vars_max(),
                                                  AHFunction::d2_vars_max()));

    int min_diagnostic_var = -1;
    int max_diagnostic_var = -1;

    for (auto &var : extra_vars)
    {
        int var_enum = std::get<0>(var.second);
        VariableType var_type = std::get<1>(var.second);

        if (var_type == VariableType::evolution)
        {
            if (var_enum < min_evolution_var)
                min_evolution_var = var_enum;
            else if (var_enum > max_evolution_var)
                max_evolution_var = var_enum;
        }
        else if (printing_step) // diagnostic
        {
            if (min_diagnostic_var == -1)
            {
                min_diagnostic_var = var_enum;
                max_diagnostic_var = var_enum;
            }
            else if (var_enum < min_diagnostic_var)
                min_diagnostic_var = var_enum;
            else if (var_enum > max_diagnostic_var)
                max_diagnostic_var = var_enum;
        }
    }

    // fill ghosts manually to minimise communication
    bool fill_ghosts = false;
    m_interpolator->refresh(fill_ghosts);
    m_interpolator->fill_multilevel_ghosts(
        VariableType::evolution,
        Interval(min_evolution_var, max_evolution_var));
    if (printing_step && min_diagnostic_var != -1)
    {
        m_interpolator->fill_multilevel_ghosts(
            VariableType::diagnostic,
            Interval(min_diagnostic_var, max_diagnostic_var));
    }
}

// 'set_coordinates' calls 'interpolate'. All Chombo_MPI:comm need to run
// 'interpolate' for it to work
template <class SurfaceGeometry, class AHFunction>
bool AHInterpolation_t<SurfaceGeometry, AHFunction>::set_coordinates(
    const vector<double> &f, const vector<double> &u,
#if CH_SPACEDIM == 3
    const vector<double> &v,
#endif
    double add_epsilon)
{
    CH_TIME("AHInterpolation::set_coordinates");

    CH_assert(f.size() == u.size());
#if CH_SPACEDIM == 3
    CH_assert(u.size() == v.size());
#endif

    const int n = f.size();

    // resize vectors
    m_f.resize(n);
    m_u.resize(n);
#if CH_SPACEDIM == 3
    m_v.resize(n);
#endif

    m_x.resize(n);
    m_y.resize(n);
#if CH_SPACEDIM == 3
    m_z.resize(n);
#endif

    bool out_of_grid = false;

    // Transform to Cartesian
    for (int i = 0; i < n; ++i)
    {
        m_f[i] = f[i] + add_epsilon;
        m_u[i] = u[i];
#if CH_SPACEDIM == 3
        m_v[i] = v[i];
#endif

        m_x[i] =
            m_coord_system.get_grid_coord(0, D_DECL(m_f[i], m_u[i], m_v[i]));
        m_y[i] =
            m_coord_system.get_grid_coord(1, D_DECL(m_f[i], m_u[i], m_v[i]));
#if CH_SPACEDIM == 3
        m_z[i] = m_coord_system.get_grid_coord(2, m_f[i], m_u[i], m_v[i]);
#endif

        // don't let PETSc diverge to outside of the grid (this can happen if
        // there is no BH)
        out_of_grid |= is_in_grid(m_x[i], m_y[i]
#if CH_SPACEDIM == 3
                                  ,
                                  m_z[i]
#endif
        );
    }

    return out_of_grid;
}

template <class SurfaceGeometry, class AHFunction>
const AHGeometryData
AHInterpolation_t<SurfaceGeometry, AHFunction>::get_geometry_data(int idx) const
{
    CH_TIME("AHInterpolation::get_geometry_data");

    return m_coord_system.get_geometry_data(
        D_DECL(m_f[idx], m_u[idx], m_v[idx]));
}

template <class SurfaceGeometry, class AHFunction>
const Tensor<1, double>
AHInterpolation_t<SurfaceGeometry, AHFunction>::get_cartesian_coords(
    int idx) const
{
    return {D_DECL(m_x[idx], m_y[idx], m_z[idx])};
}

template <class SurfaceGeometry, class AHFunction>
const Tensor<1, double>
AHInterpolation_t<SurfaceGeometry, AHFunction>::get_coords(int idx) const
{
    return
    {
        m_u[idx],
#if CH_SPACEDIM == 3
            m_v[idx],
#endif
            m_f[idx]
    };
}

template <class SurfaceGeometry, class AHFunction>
const AHVarsData<int, double>
AHInterpolation_t<SurfaceGeometry, AHFunction>::get_data(int idx) const
{
    return get_AHVarsData_idx(idx, m_data);
}

template <class SurfaceGeometry, class AHFunction>
void AHInterpolation_t<SurfaceGeometry, AHFunction>::interpolate_extra_vars(
    const std::map<std::string, std::tuple<int, VariableType, int>> &extra_vars)
{
    CH_TIME("AHInterpolation::interpolate_extra_vars");

    if (extra_vars.size() == 0)
        return;

    int n = m_x.size();

    InterpolationQuery query(n);
    query.setCoords(0, m_x.data()).setCoords(1, m_y.data());
#if CH_SPACEDIM == 3
    query.setCoords(2, m_z.data());
#endif

    for (auto &var : extra_vars)
    {
        int var_enum = std::get<0>(var.second);
        VariableType var_type = std::get<1>(var.second);
        int der_type = std::get<2>(var.second);

        if (der_type == 0)
            m_extra.set_vars(query, var_enum, var.first, var_type, n);
        else if (der_type == 1)
            m_extra.set_d1(query, var_enum, var.first, var_type, n);
        else // if(der_type == 2)
            m_extra.set_d2(query, var_enum, var.first, var_type, n);
    }

    // ghosts for extra evolution or diagnostic vars already filled in
    // 'refresh_interpolator'
    m_interpolator->interp(query);
}

template <class SurfaceGeometry, class AHFunction>
const AHVarsData<std::string, double>
AHInterpolation_t<SurfaceGeometry, AHFunction>::get_extra_data(int idx) const
{
    return get_AHVarsData_idx(idx, m_extra);
}

#endif // _AHINTERPOLATION_IMPL_HPP_