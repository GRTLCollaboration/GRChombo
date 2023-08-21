/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHINTERPOLATION_HPP_
#define _AHINTERPOLATION_HPP_

#include "AHGeometryData.hpp"
#include "AHVarsData.hpp"
#include "AMRInterpolator.hpp"
#include "Lagrange.hpp"

#include <map>
#include <string>

//! Class used for interpolation of the variables needed to calculate expansion
//! with the data from a given 'SurfaceGeometry'
template <class SurfaceGeometry, class AHFunction> class AHInterpolation_t
{
  private:
    SurfaceGeometry m_coord_system;
    AMRInterpolator<Lagrange<4>> *m_interpolator;

    // variables of AH in 'SurfaceGeometry'
    // (in Spherical Coordinates Class 'AHSphericalCoords',
    // they correspond to 'theta', 'phi' and 'sqrt(radius)')
    std::vector<double> m_u;
#if CH_SPACEDIM == 3
    std::vector<double> m_v;
#endif
    std::vector<double> m_f;

    // variables of AH in cartesian (for AMR interpolation)
    std::vector<double> m_x;
    std::vector<double> m_y;
#if CH_SPACEDIM == 3
    std::vector<double> m_z;
#endif

    AHVarsData<int, std::vector<double>> m_data;
    AHVarsData<std::string, std::vector<double>> m_extra;

    std::array<double, CH_SPACEDIM> m_coord_min,
        m_coord_max; //!< maximum and minimum of level 0 box, used in
                     //!< 'fit_in_grid'

  public:
    AHInterpolation_t(const SurfaceGeometry &a_coordSystem,
                      AMRInterpolator<Lagrange<4>> *a_interpolator);

    const AMRInterpolator<Lagrange<4>> *get_interpolator() const;
    const SurfaceGeometry &get_coord_system() const;

    std::vector<std::string> get_labels() const; //!< get all names (u, v, f)

    void set_origin(
        const std::array<double, CH_SPACEDIM> &); //!< set origin of CoordSystem

    void refresh_interpolator(
        bool printing_step,
        const std::map<std::string, std::tuple<int, VariableType, int>>
            &extra_vars); //!< refresh AMRInterpolator 'm_interpolator'

    //! returns whether any pointis outside of grid (==> diverging)
    bool set_coordinates(D_DECL(const std::vector<double> &f,
                                const std::vector<double> &u,
                                const std::vector<double> &v),
                         double add_epsilon = 0.);
    const AHGeometryData get_geometry_data(int idx) const;
    const Tensor<1, double> get_cartesian_coords(int idx) const;
    const Tensor<1, double> get_coords(int idx) const;
    const AHVarsData<int, double> get_data(int idx) const;

    //! verify point is inside the grid
    //! when PETSc tried to diverge out of the grid, this doesn't let him do so
    //! it forces him to stay on the grid, causing non-convergence
    bool is_in_grid(D_DECL(double &x, double &y, double &z));
    //! PETSc AND AMRInterpolator to interpolate
    //! 'keep_interpolating_if_inactive' returns 'true' immediately for PETSc
    //! ranks for non-PETSc ranks, it loops in a sequence of 'interpolate()'
    //! waiting for PETSc ranks to call 'interpolate()'' as well (as this needs
    //! to be called by all Chombo processes) can be aborted by calling
    //! 'break_interpolation_loop()' for PETSc ranks
    /** Example of usage:
     * if(m_geom.keep_interpolating_if_inactive())
     * {
     *     (... PETSc code that calls 'set_coordinates()'...)
     *     m_geom.break_interpolation_loop();
     * }
     */
    bool keep_interpolating_if_inactive();
    void
    break_interpolation_loop() const; //!< see 'keep_interpolating_if_inactive'
                                      //!< (this one could be static)
    //! returns whether or not to keep interpolating. See
    //! 'keep_interpolating_if_inactive'
    //! Again, all Chombo_MPI:comm need to run 'interpolate' for it to work (if
    //! some are not part of PETSc, they still need to run 'interpolate' on the
    //! side)
    int interpolate();

    void interpolate_extra_vars(
        const std::map<std::string, std::tuple<int, VariableType, int>>
            &extra_vars);
    const AHVarsData<std::string, double> get_extra_data(int idx) const;
};

#include "AHInterpolation.impl.hpp"

#endif /* _AHINTERPOLATION_HPP_ */
