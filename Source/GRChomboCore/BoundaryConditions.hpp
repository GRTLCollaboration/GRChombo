/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOUNDARYCONDITIONS_HPP_
#define BOUNDARYCONDITIONS_HPP_

#include "BoxIterator.H"
#include "Coordinates.hpp"
#include "Copier.H"
#include "DimensionDefinitions.hpp"
#include "FourthOrderInterpStencil.H"
#include "GRLevelData.hpp"
#include "Interval.H"
#include "RealVect.H"
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
        REFLECTIVE_BC,
        EXTRAPOLATING_BC,
        MIXED_BC
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
        std::vector<int> mixed_bc_extrapolating_vars;
        std::vector<int> mixed_bc_sommerfeld_vars;
        int extrapolation_order = 1;
    };

  protected:
    // Member values
    double m_dx;                 // The grid spacing
    int m_num_ghosts;            // the number of ghosts (usually 3)
    params_t m_params;           // the boundary params
    std::vector<int> m_all_vars; // a vector of c_nums for all the vars
    RealVect m_center;           // the position of the center of the grid
    ProblemDomain m_domain;      // the problem domain (excludes boundary cells)
    Box m_domain_box;            // The box representing the domain
    bool is_defined; // whether the BoundaryConditions class members are defined

  public:
    /// Default constructor - need to call define afterwards
    BoundaryConditions() { is_defined = false; }

    /// define function sets members and is_defined set to true
    void define(double a_dx, std::array<double, CH_SPACEDIM> a_center,
                params_t a_params, ProblemDomain a_domain, int a_num_ghosts);

    /// change the asymptotic values of the variables for the Sommerfeld BCs
    /// this will allow them to evolve during a simulation if necessary
    void set_vars_asymptotic_values(
        std::array<double, NUM_VARS> &vars_asymptotic_values);

    /// write out boundary params (used during setup for debugging)
    static void write_boundary_conditions(params_t a_params);

    /// The function which returns the parity of each of the vars in
    /// UserVariables.hpp The parity should be defined in the params file, and
    /// will be output to the pout files for checking at start/restart of
    /// simulation (It is only required for reflective boundary conditions.)
    int get_vars_parity(int a_comp, int a_dir) const;

    /// static version used for initial output of boundary values
    static int get_vars_parity(int a_comp, int a_dir, params_t a_params);

    /// Fill the rhs boundary values appropriately based on the params set
    void fill_boundary_rhs(const Side::LoHiSide a_side,
                           const GRLevelData &a_soln, GRLevelData &a_rhs);

    /// Fill the boundary values appropriately based on the params set
    /// in the direction dir
    void fill_boundary_cells_dir(const Side::LoHiSide a_side,
                                 const GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const int dir, const bool filling_rhs = true);

    /// Copy the boundary values from src to dest
    /// NB assumes same box layout of input and output data
    void copy_boundary_cells(const Side::LoHiSide a_side,
                             const GRLevelData &a_src, GRLevelData &a_dest);

    /// enforce solution boundary conditions, e.g. after interpolation
    /// for BCs where solution is specified, not rhs
    void enforce_solution_boundaries(const Side::LoHiSide a_side,
                                     GRLevelData &a_state);

    /// Fill the fine boundary values in a_state
    /// Required for interpolating onto finer levels at boundaries
    void interp_boundaries(GRLevelData &a_fine_state,
                           GRLevelData &a_coarse_state,
                           const Side::LoHiSide a_side);

    /// Get the boundary condition for a_dir and a_side
    int get_boundary_condition(const Side::LoHiSide a_side, const int a_dir);

    /// get the boundary box to fill if we are at a boundary
    Box get_boundary_box(const Side::LoHiSide a_side, const int a_dir,
                         const IntVect &offset_lo, const IntVect &offset_hi,
                         Box &this_ghostless_box, int shrink_for_coarse = 0);

    /// This function takes a default constructed open DisjointBoxLayout and
    /// grows the boxes lying along the boundary to include the boundaries if
    /// necessary (i.e. in the Sommerfeld BC case). It is used to define the
    /// correct DisjointBoxLayout for the exchange copier so that shared
    /// boundary ghosts are exchanged correctly.
    void expand_grids_to_boundaries(DisjointBoxLayout &a_out_grids,
                                    const DisjointBoxLayout &a_in_grids);

    friend class ExpandGridsToBoundaries;

  private:
    /// write out reflective conditions
    static void write_reflective_conditions(int idir, params_t a_params);

    /// write out sommerfeld conditions
    static void write_sommerfeld_conditions(int idir, params_t a_params);

    /// write out mixed conditions
    static void write_mixed_conditions(int idir, params_t a_params);

    void fill_sommerfeld_cell(FArrayBox &rhs_box, const FArrayBox &soln_box,
                              const IntVect iv,
                              const std::vector<int> &sommerfeld_comps) const;

    void fill_reflective_cell(FArrayBox &rhs_box, const IntVect iv,
                              const Side::LoHiSide a_side, const int dir,
                              const std::vector<int> &reflective_comps) const;

    void fill_extrapolating_cell(FArrayBox &rhs_box, const IntVect iv,
                                 const Side::LoHiSide a_side, const int dir,
                                 const std::vector<int> &extrapolating_comps,
                                 const int order = 1) const;
};

/// This derived class is used by expand_grids_to_boundaries to grow the
/// boxes along the Sommerfeld BC boundaries
class ExpandGridsToBoundaries : public BaseTransform
{
  public:
    ExpandGridsToBoundaries(BoundaryConditions &a_boundaries)
        : m_boundaries(a_boundaries)
    {
    }

    /// Operator called by transform to grow the boxes where required
    Box operator()(const Box &a_in_box) override;

  protected:
    BoundaryConditions &m_boundaries;
};

#endif /* BOUNDARYCONDITIONS_HPP_ */
