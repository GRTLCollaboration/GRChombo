/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOUNDARYCONDITIONS_HPP_
#define BOUNDARYCONDITIONS_HPP_

#include "BoxIterator.H"
#include "Copier.H"
#include "DimensionDefinitions.hpp"
#include "FourthOrderInterpStencil.H"
#include "GRLevelData.hpp"
#include "UserVariables.hpp"

/// Class which deals with the boundaries at the edge of the physical domain in
/// cases where they are not periodic. Currently only options are static BCs,
/// sommerfeld (outgoing radiation) and symmetric. The conditions can differ in
/// the high and low directions.
/// In cases where different variables/boundaries are required, the user should
/// (usually) write their own conditions class which inherits from this one.
/// Note that these conditions enforce a certain rhs based on the current values
/// of the grid variables. (Another option would be to enforce grid values, e.g.
/// by extrapolating from within the grid.)
class BoundaryConditions
{
  public:
    /// enum for possible boundary states
    enum
    {
        STATIC_BC,
        SOMMERFELD_BC,
        REFLECTIVE_BC
    };

    /// enum for possible parity states
    enum
    {
        EVEN,
        ODD_X,
        ODD_Y,
        ODD_Z,
        ODD_XY,
        ODD_YZ,
        ODD_XZ
    };

    /// Structure containing the boundary condition params
    struct params_t
    {
        std::array<int, CH_SPACEDIM> hi_boundary;
        std::array<int, CH_SPACEDIM> lo_boundary;
        std::array<bool, CH_SPACEDIM> is_periodic;
        std::array<int, NUM_VARS> vars_parity;
        std::array<double, NUM_VARS> vars_asymptotic_values;
    };

  protected:
    // Member values
    double m_dx;            // The grid spacing
    int m_num_ghosts;       // the number of ghosts (usually 3)
    params_t m_params;      // the boundary params
    RealVect m_center;      // the position of the center of the grid
    ProblemDomain m_domain; // the problem domain (excludes boundary cells)
    Box m_domain_box;       // The box representing the domain
    bool is_defined; // whether the BoundaryConditions class members are defined

  public:
    /// Default constructor - need to call define afterwards
    BoundaryConditions() { is_defined = false; }

    /// define function sets members and is_defined set to true
    void define(double a_dx, std::array<double, CH_SPACEDIM> a_center,
                params_t a_params, ProblemDomain a_domain, int a_num_ghosts)
    {
        m_dx = a_dx;
        m_params = a_params;
        m_domain = a_domain;
        m_domain_box = a_domain.domainBox();
        m_num_ghosts = a_num_ghosts;
        FOR1(i) { m_center[i] = a_center[i]; }
        is_defined = true;
    }

    // change the asymptotic values of the variables for the Sommerfeld BCs
    // this will allow them to evolve during a simulation if necessary
    void set_vars_asymptotic_values(
        std::array<double, NUM_VARS> &vars_asymptotic_values)
    {
        m_params.vars_asymptotic_values = vars_asymptotic_values;
    }

    // write out boundary params (used during setup for debugging)
    static void write_boundary_conditions(params_t a_params)
    {
        pout() << "You are using non periodic boundary conditions." << endl;
        pout() << "The boundary params chosen are:  " << endl;
        pout() << "---------------------------------" << endl;
        FOR1(idir)
        {
            if (!a_params.is_periodic[idir])
            {
                // high directions
                if (a_params.hi_boundary[idir] == STATIC_BC)
                {
                    pout() << "- Static boundaries in direction high " << idir
                           << endl;
                }
                else if (a_params.hi_boundary[idir] == REFLECTIVE_BC)
                {
                    pout() << "- Reflective boundaries in direction high "
                           << idir << endl;
                    pout() << "The variables that are parity odd in this "
                              "direction are : "
                           << endl;
                    for (int icomp = 0; icomp < NUM_VARS; icomp++)
                    {
                        int parity = get_vars_parity(icomp, idir, a_params);
                        if (parity == -1)
                        {
                            pout() << UserVariables::variable_names[icomp]
                                   << "    ";
                        }
                    }
                    pout() << endl;
                }
                else if (a_params.hi_boundary[idir] == SOMMERFELD_BC)
                {
                    pout() << "- Sommerfeld boundaries in direction high "
                           << idir << endl;
                    pout() << "The non zero asymptotic values of the variables "
                              "in this direction are : "
                           << endl;
                    for (int icomp = 0; icomp < NUM_VARS; icomp++)
                    {
                        if (a_params.vars_asymptotic_values[icomp] != 0)
                        {
                            pout()
                                << UserVariables::variable_names[icomp] << " = "
                                << a_params.vars_asymptotic_values[icomp]
                                << "    ";
                        }
                    }
                    pout() << endl;
                }

                // low directions
                if (a_params.lo_boundary[idir] == STATIC_BC)
                {
                    pout() << "- Static boundaries in direction low " << idir
                           << endl;
                }
                else if (a_params.lo_boundary[idir] == REFLECTIVE_BC)
                {
                    pout() << "- Reflective boundaries in direction low "
                           << idir << endl;
                    pout() << "The variables that are parity odd in this "
                              "direction are : "
                           << endl;
                    for (int icomp = 0; icomp < NUM_VARS; icomp++)
                    {
                        int parity = get_vars_parity(icomp, idir, a_params);
                        if (parity == -1)
                        {
                            pout() << UserVariables::variable_names[icomp]
                                   << "    ";
                        }
                    }
                    pout() << endl;
                }
                else if (a_params.lo_boundary[idir] == SOMMERFELD_BC)
                {
                    pout() << "- Sommerfeld boundaries in direction low "
                           << idir << endl;
                    pout() << "The non zero asymptotic values of the variables "
                              "in this direction are : "
                           << endl;
                    for (int icomp = 0; icomp < NUM_VARS; icomp++)
                    {
                        if (a_params.vars_asymptotic_values[icomp] != 0)
                        {
                            pout()
                                << UserVariables::variable_names[icomp] << " = "
                                << a_params.vars_asymptotic_values[icomp]
                                << "    ";
                        }
                    }
                    pout() << endl;
                }
            }
        }
        pout() << "---------------------------------" << endl;
    }

    // The function which returns the parity of each of the vars in
    // UserVariables.hpp The parity should be defined in the params file, and
    // will be output to the pout files for checking at start/restart of
    // simulation (It is only required for reflective boundary conditions.)
    int get_vars_parity(int a_comp, int a_dir)
    {
        int vars_parity = get_vars_parity(a_comp, a_dir, m_params);

        return vars_parity;
    }

    // static version used for initial output of boundary values
    static int get_vars_parity(int a_comp, int a_dir, params_t a_params)
    {
        int vars_parity = 1;
        if ((a_dir == 0) && (a_params.vars_parity[a_comp] == ODD_X ||
                             a_params.vars_parity[a_comp] == ODD_XY ||
                             a_params.vars_parity[a_comp] == ODD_XZ))
        {
            vars_parity = -1;
        }
        else if ((a_dir == 1) && (a_params.vars_parity[a_comp] == ODD_Y ||
                                  a_params.vars_parity[a_comp] == ODD_XY ||
                                  a_params.vars_parity[a_comp] == ODD_YZ))
        {
            vars_parity = -1;
        }
        else if ((a_dir == 2) && (a_params.vars_parity[a_comp] == ODD_Z ||
                                  a_params.vars_parity[a_comp] == ODD_XZ ||
                                  a_params.vars_parity[a_comp] == ODD_YZ))
        {
            vars_parity = -1;
        }
        return vars_parity;
    }

    /// Fill the rhs boundary values appropriately based on the params set
    void fill_boundary_rhs(const Side::LoHiSide a_side, GRLevelData &a_soln,
                           GRLevelData &a_rhs)
    {
        CH_assert(is_defined);
        CH_TIME("BoundaryConditions::fill_boundary_rhs");

        // cycle through the directions, filling the rhs
        FOR1(idir)
        {
            // only do something if this direction is not periodic
            if (!m_params.is_periodic[idir])
            {
                fill_boundary_rhs_dir(a_side, a_soln, a_rhs, idir);
            }
        }
    }

    /// Fill the boundary values appropriately based on the params set
    // in the direction dir
    void fill_boundary_rhs_dir(const Side::LoHiSide a_side, GRLevelData &a_soln,
                               GRLevelData &a_rhs, const int dir)
    {
        // iterate through the boxes, shared amongst threads
        DataIterator dit = a_rhs.dataIterator();
        int nbox = dit.size();
#pragma omp parallel for default(shared)
        for (int ibox = 0; ibox < nbox; ++ibox)
        {
            DataIndex dind = dit[ibox];
            FArrayBox &m_rhs_box = a_rhs[dind];
            FArrayBox &m_soln_box = a_soln[dind];
            Box this_box = m_rhs_box.box();
            IntVect offset_lo = -this_box.smallEnd() + m_domain_box.smallEnd();
            IntVect offset_hi = +this_box.bigEnd() - m_domain_box.bigEnd();

            // reduce box to the intersection of the box and the
            // problem domain ie remove all outer ghost cells
            this_box &= m_domain_box;
            // get the boundary box (may be Empty) and the condition on it
            int boundary_condition = get_boundary_condition(a_side, dir);
            Box boundary_box =
                get_boundary_box(a_side, dir, offset_lo, offset_hi, this_box);

            // now we have the appropriate box, fill it!
            BoxIterator bit(boundary_box);
            for (bit.begin(); bit.ok(); ++bit)
            {
                IntVect iv = bit();
                switch (boundary_condition)
                {
                // simplest case - boundary values are fixed to the initial ones
                case STATIC_BC:
                {
                    for (int icomp = 0; icomp < NUM_VARS; icomp++)
                    {
                        m_rhs_box(iv, icomp) = 0.0;
                    }
                    break;
                }
                // assumes an asymptotic value + radial waves and permits them
                // to exit grid with minimal reflections
                case SOMMERFELD_BC:
                {
                    // get real position on the grid
                    RealVect loc(iv + 0.5 * RealVect::Unit);
                    loc *= m_dx;
                    loc -= m_center;
                    double radius_squared = 0.0;
                    FOR1(i) { radius_squared += loc[i] * loc[i]; }
                    double radius = sqrt(radius_squared);
                    IntVect lo_local_offset = iv - m_soln_box.smallEnd();
                    IntVect hi_local_offset = m_soln_box.bigEnd() - iv;

                    // Apply Sommerfeld BCs to each variable
                    for (int icomp = 0; icomp < NUM_VARS; icomp++)
                    {
                        m_rhs_box(iv, icomp) = 0.0;
                        FOR1(idir2)
                        {
                            IntVect iv_offset1 = iv;
                            IntVect iv_offset2 = iv;
                            double d1;
                            // bit of work to get the right stencils for near
                            // the edges of the domain, only using second order
                            // stencils for now
                            if (lo_local_offset[idir2] < 1)
                            {
                                // near lo end
                                iv_offset1[idir2] += +1;
                                iv_offset2[idir2] += +2;
                                d1 = 1.0 / m_dx *
                                     (-1.5 * m_soln_box(iv, icomp) +
                                      2.0 * m_soln_box(iv_offset1, icomp) -
                                      0.5 * m_soln_box(iv_offset2, icomp));
                            }
                            else if (hi_local_offset[idir2] < 1)
                            {
                                // near hi end
                                iv_offset1[idir2] += -1;
                                iv_offset2[idir2] += -2;
                                d1 = 1.0 / m_dx *
                                     (+1.5 * m_soln_box(iv, icomp) -
                                      2.0 * m_soln_box(iv_offset1, icomp) +
                                      0.5 * m_soln_box(iv_offset2, icomp));
                            }
                            else
                            {
                                // normal case
                                iv_offset1[idir2] += +1;
                                iv_offset2[idir2] += -1;
                                d1 = 0.5 / m_dx *
                                     (m_soln_box(iv_offset1, icomp) -
                                      m_soln_box(iv_offset2, icomp));
                            }

                            // for each direction add dphidx * x^i / r
                            m_rhs_box(iv, icomp) += -d1 * loc[idir2] / radius;
                        }

                        // asymptotic values - these need to have been set in
                        // the params file
                        m_rhs_box(iv, icomp) +=
                            (m_params.vars_asymptotic_values[icomp] -
                             m_soln_box(iv, icomp)) /
                            radius;
                    }
                    break;
                }
                // assume boundary is a reflection of values within the grid
                // care must be taken with variable parity to maintain correct
                // values on reflection, e.g. x components of vectors are odd
                // parity in the x direction
                case REFLECTIVE_BC:
                {
                    IntVect iv_copy = iv;
                    /// where to copy the data from - mirror image in domain
                    if (a_side == Side::Lo)
                    {
                        iv_copy[dir] = -iv[dir] - 1;
                    }
                    else
                    {
                        iv_copy[dir] =
                            2 * m_domain_box.bigEnd(dir) - iv[dir] + 1;
                    }

                    // replace value at iv with value at iv_copy
                    for (int icomp = 0; icomp < NUM_VARS; icomp++)
                    {
                        int parity = get_vars_parity(icomp, dir);
                        m_rhs_box(iv, icomp) =
                            parity * m_rhs_box(iv_copy, icomp);
                    }
                    break;
                }
                default:
                    MayDay::Error(
                        "BoundaryCondition::Supplied boundary not supported.");
                } // end switch
            }     // end iterate over box
        }         // end iterate over boxes
    }

    /// Copy the boundary values from src to dest
    /// NB assumes same box layout of input and output data
    void copy_boundary_cells(const Side::LoHiSide a_side,
                             const GRLevelData &a_src, GRLevelData &a_dest)
    {
        CH_TIME("BoundaryConditions::copy_boundary_cells");

        CH_assert(is_defined);
        if (a_src.boxLayout() == a_dest.boxLayout())
        {
            // cycle through the directions
            FOR1(idir)
            {
                // only do something if this direction is not periodic
                if (!m_params.is_periodic[idir])
                {
                    // iterate through the boxes, shared amongst threads
                    DataIterator dit = a_dest.dataIterator();
                    int nbox = dit.size();
#pragma omp parallel for default(shared)
                    for (int ibox = 0; ibox < nbox; ++ibox)
                    {
                        DataIndex dind = dit[ibox];
                        FArrayBox &m_dest_box = a_dest[dind];
                        Box this_box = m_dest_box.box();
                        IntVect offset_lo =
                            -this_box.smallEnd() + m_domain_box.smallEnd();
                        IntVect offset_hi =
                            +this_box.bigEnd() - m_domain_box.bigEnd();

                        // reduce box to the intersection of the box and the
                        // problem domain ie remove all outer ghost cells
                        this_box &= m_domain_box;

                        // get the boundary box (may be Empty)
                        Box boundary_box = get_boundary_box(
                            a_side, idir, offset_lo, offset_hi, this_box);

                        BoxIterator bit(boundary_box);
                        for (bit.begin(); bit.ok(); ++bit)
                        {
                            IntVect iv = bit();
                            for (int icomp = 0; icomp < NUM_VARS; icomp++)
                            {
                                m_dest_box(iv, icomp) = a_src[dind](iv, icomp);
                            }
                        } // end iterate over box
                    }     // end iterate over boxes
                }         // end if(not periodic)
            }             // end iterate over spacedims
        }                 // end test for same box layout
    }

    /// enforce symmetric boundary conditions, e.g. after interpolation
    void enforce_symmetric_boundaries(const Side::LoHiSide a_side,
                                      GRLevelData &a_state)
    {
        CH_assert(is_defined);
        CH_TIME("BoundaryConditions::enforce_symmetric_boundaries");

        // cycle through the directions
        FOR1(idir)
        {
            // only do something if this direction is not periodic and symmetric
            if (!m_params.is_periodic[idir])
            {
                int boundary_condition = get_boundary_condition(a_side, idir);

                // as a bit of a hack, just use the rhs update, since it is the
                // same copying of cells which we require for the solution
                if (boundary_condition == REFLECTIVE_BC)
                {
                    fill_boundary_rhs_dir(a_side, a_state, a_state, idir);
                }
            }
        }
    }

    /// Fill the fine boundary values in a_state
    /// Required for interpolating onto finer levels at boundaries
    void interp_boundaries(GRLevelData &a_fine_state,
                           GRLevelData &a_coarse_state,
                           const Side::LoHiSide a_side)
    {
        CH_assert(is_defined);
        CH_TIME("BoundaryConditions::interp_boundaries");

        // cycle through the directions
        FOR1(idir)
        {
            // only do something if this direction is not periodic
            if (!m_params.is_periodic[idir])
            {
                // Ref ratio is always two
                int ref_ratio = 2;

                // create a coarsened fine layout and copy the coarse data onto
                // it
                DisjointBoxLayout coarsened_layout;
                coarsen(coarsened_layout, a_fine_state.disjointBoxLayout(),
                        ref_ratio * IntVect::Unit);
                GRLevelData coarsened_fine;
                coarsened_fine.define(coarsened_layout, NUM_VARS,
                                      m_num_ghosts * IntVect::Unit);
                Box coarse_domain_box = coarsen(m_domain_box, ref_ratio);

                // trick the copyTo into thinking the boundary cells are within
                // the domain by growing the domain
                Box grown_domain_box = coarse_domain_box;
                grown_domain_box.grow(m_num_ghosts * IntVect::Unit);
                Copier boundary_copier;
                boundary_copier.ghostDefine(
                    a_coarse_state.disjointBoxLayout(),
                    coarsened_fine.disjointBoxLayout(), grown_domain_box,
                    m_num_ghosts * IntVect::Unit, m_num_ghosts * IntVect::Unit);
                a_coarse_state.copyTo(a_coarse_state.interval(), coarsened_fine,
                                      coarsened_fine.interval(),
                                      boundary_copier);

                // iterate through the coarse boxes, shared amongst threads
                DataIterator dit = coarsened_layout.dataIterator();
                int nbox = dit.size();
#pragma omp parallel for default(shared)
                for (int ibox = 0; ibox < nbox; ++ibox)
                {
                    DataIndex dind = dit[ibox];
                    FArrayBox &m_fine_box = a_fine_state[dind];
                    FArrayBox &m_coarse_box = coarsened_fine[dind];
                    Box this_box = m_coarse_box.box();
                    Box fine_box = m_fine_box.box();
                    IntVect offset_lo =
                        -this_box.smallEnd() + coarse_domain_box.smallEnd();
                    IntVect offset_hi =
                        +this_box.bigEnd() - coarse_domain_box.bigEnd();

                    // reduce box to the intersection of the box and the
                    // problem domain ie remove all outer ghost cells
                    this_box &= coarse_domain_box;

                    // get the boundary box - remove one cell as we only want 2
                    // coarse cells filled in each direction, to fill the 3 fine
                    // cells on the level above
                    Box boundary_box = get_boundary_box(a_side, idir, offset_lo,
                                                        offset_hi, this_box, 1);

                    // define standard stencil for interp where not near
                    // boundaries in other dirs
                    IntVect default_offset =
                        IntVect::Zero + sign(a_side) * 2 * BASISV(idir);
                    FourthOrderInterpStencil default_stencil(default_offset,
                                                             ref_ratio);

                    // now interp the box from coarse to fine
                    BoxIterator bit(boundary_box);
                    for (bit.begin(); bit.ok(); ++bit)
                    {
                        IntVect iv = bit();
                        IntVect lo_local_offset = iv - m_coarse_box.smallEnd();
                        IntVect hi_local_offset = m_coarse_box.bigEnd() - iv;

                        // bit of work to get the right stencils for near the
                        // edges of the box
                        bool near_boundary = false;
                        IntVect local_boundary_offset = IntVect::Zero;
                        FOR1(idir2)
                        {
                            if (idir2 == idir)
                            {
                                local_boundary_offset[idir2] =
                                    default_offset[idir2];
                            }
                            else if ((idir2 != idir) &&
                                     (lo_local_offset[idir2] > 1) &&
                                     (hi_local_offset[idir2] > 1))
                            {
                                local_boundary_offset[idir2] = 0;
                            }
                            else if ((idir2 != idir) &&
                                     (lo_local_offset[idir2] == 1))
                            {
                                local_boundary_offset[idir2] = -2;
                                near_boundary = true;
                            }
                            else if ((idir2 != idir) &&
                                     (hi_local_offset[idir2] == 1))
                            {
                                local_boundary_offset[idir2] = +2;
                                near_boundary = true;
                            }
                            else
                            {
                                MayDay::Error(
                                    "BoundaryConditions::define bad boxes");
                            }
                        }

                        // if not near the boundary use the default stencil,
                        // otherwise use the one calculated locally
                        if (!near_boundary)
                        {
                            default_stencil.fillFine(m_fine_box, m_coarse_box,
                                                     iv);
                        }
                        else
                        {
                            FourthOrderInterpStencil local_stencil(
                                local_boundary_offset, ref_ratio);
                            local_stencil.fillFine(m_fine_box, m_coarse_box,
                                                   iv);
                        }
                    } // end loop box
                }     // end loop boxes
            }         // end if is_periodic
        }             // end loop idir
    }

    /// Get the boundary condition for a_dir and a_side
    int get_boundary_condition(const Side::LoHiSide a_side, const int a_dir)
    {
        int boundary_condition = 0;
        if (a_side == Side::Lo)
        {
            boundary_condition = m_params.lo_boundary[a_dir];
        }
        else
        {
            boundary_condition = m_params.hi_boundary[a_dir];
        }
        return boundary_condition;
    }

    /// get the boundary box to fill if we are at a boundary
    Box get_boundary_box(const Side::LoHiSide a_side, const int a_dir,
                         const IntVect &offset_lo, const IntVect &offset_hi,
                         Box &this_ghostless_box, int shrink_for_coarse = 0)
    {
        // default constructor gives empty box
        Box boundary_box;

        // check if we are over the edges of the domain - are we a boundary box?
        // if so create the box of the cells we want to fill
        if (((a_side == Side::Hi) && (offset_hi[a_dir] > 0)) ||
            ((a_side == Side::Lo) && (offset_lo[a_dir] > 0)))
        {
            // Get just the boundary box to iterate over, m_num_ghosts ghost
            // cells unless we are filling the coarse cells in the interp case
            // where we want to fill only two coarse ghost cells (to cover 3
            // fine ones)
            if (a_side == Side::Lo)
            {
                boundary_box = adjCellLo(this_ghostless_box, a_dir,
                                         m_num_ghosts - shrink_for_coarse);
            }
            else
            {
                boundary_box = adjCellHi(this_ghostless_box, a_dir,
                                         m_num_ghosts - shrink_for_coarse);
            }

            // adjust for any offsets - catches the corners etc
            // but only want to fill them once, so x fills y and z, y fills z
            // etc. Required in periodic direction corners in cases where there
            // are mixed boundaries, (otherwise these corners are full of nans)
            FOR1(idir)
            {
                if (offset_lo[idir] > 0) // this direction is a low end boundary
                {
                    if ((idir > a_dir) || (m_params.is_periodic[idir]))
                    {
                        // grow it to fill the corners
                        boundary_box.growLo(idir,
                                            m_num_ghosts - shrink_for_coarse);
                    }
                }
                else // cut off end ghost cell
                {
                    if (idir != a_dir)
                    {
                        boundary_box.growLo(idir, -shrink_for_coarse);
                    }
                }

                if (offset_hi[idir] > 0) // this direction is a high end
                                         // boundary
                {
                    if ((idir > a_dir) || (m_params.is_periodic[idir]))
                    {
                        // grow it to fill the corners
                        boundary_box.growHi(idir,
                                            m_num_ghosts - shrink_for_coarse);
                    }
                }
                else // cut off end ghost cell
                {
                    if (idir != a_dir)
                    {
                        boundary_box.growHi(idir, -shrink_for_coarse);
                    }
                }
            }
        }
        return boundary_box;
    }
};

#endif /* BOUNDARYCONDITIONS_HPP_ */
