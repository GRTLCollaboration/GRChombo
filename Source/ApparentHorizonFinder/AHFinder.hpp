/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHFINDER_HPP_
#define _AHFINDER_HPP_

// Chombo includes
#include "SPMD.H" //Chombo_MPI::comm

// Other includes:

#ifdef CH_MPI
#include <mpi.h>
#endif
#include <petsc.h>

#include "AHInterpolation.hpp"
#include "AHParams.hpp"
#include "ApparentHorizon.hpp"
#include "ChomboParameters.hpp"
#include "Lagrange.hpp"

// Chombo namespace
#include "UsingNamespace.H"

/////////////////////////////////////////////////////////////////////////////
// SurfaceGeometry and AHFunction defauls

// included so that user can re-define default AHSurfaceGeometry or AHFunction
#include "UserVariables.hpp"

// define SurfaceGeometry of AHFinder
#ifndef AHSurfaceGeometry
#include "AHSphericalGeometry.hpp"
#define AHSurfaceGeometry AHSphericalGeometry
#endif

// default to expansion
#ifndef AHFunction
#include "AHFunctions.hpp"
#define AHFunction ExpansionFunction
#endif

#include "AHInitialGuess.hpp"
// // default to constant initial guess
// #ifndef AHInitialGuess
// #include "AHInitialGuess.hpp"
// #define AHInitialGuess AHInitialGuessConstant
// #endif
/////////////////////////////////////////////////////////////////////////////

/*
TF:
A general note on the radius and the value of chi on the AH surface.
For a Kerr BH initial data with dimensionless spin 's' and mass 'M':
  r_AH   = M / 4 (1 + sqrt{1 - s^2})
  chi_AH = [(1-s^2)^(1/6) / 16 * ((1 + sqrt{1-s^2})/2)^{2/3} for theta=0,
            (1-s^2)^{1/6} / 16 * ((1 + sqrt{1-s^2})/2)^{1/3} for theta=pi/2]

After puncture gauge settled:
  r_AH   ~ M * (0.275 + 0.76 * sqrt{1-s^2})
  chi_AH ~ 0.26 * sqrt{1-s^2}

Take it with a grain of salt, but these results are from a numerical fit I did
for numerical Kerr BH runs, approximately settled at around t~50M, from spin=0
to spin=0.99. You can use it generically for other things, like realizing that
for spin=0 you can plot in VisIt the contour chi=0.26 to get the AH, or that for
spin=0 the AH goes from r=M/2 to r~M

Binaries:
After some inspiral cases I've been through, I realized that merger happened
when the distance between the punctures was about (roughtly!) ~[0.5,2]*(M1+M2)
and with a radius of ~[0.7,1.1]*(M1+M2) (wobbly of course, and bigger values for
smaller final spin). It makes sense to start looking for the merger before this,
but remember that for PETSc diverging is significantly slower than converging,
so the closer the better. Hence, as default, start looking for a merger when the
separation is about ~<2*(M1+M2) and with a initial guess of ~2*(M1+M2) (twice as
big as what we'll probably get to make sure it converges to the outer AH; even
then it might converge to the inner one, but it is much easer for the AH to
converge if the initial guess is bigger than the actual radius, but the opposite
is more senstitive).
*/

//! Class to manage AHs and its mergers + control PETSc MPI sub-communicator
template <class SurfaceGeometry = AHSurfaceGeometry,
          class AHFunction = AHFunction>
class AHFinder
{
    using AHInterpolation = AHInterpolation_t<SurfaceGeometry, AHFunction>;
    using AHParams = AHParams_t<AHFunction>;

  public:
    AHFinder(){};
    ~AHFinder();

    ALWAYS_INLINE ApparentHorizon<SurfaceGeometry, AHFunction> *
    get(unsigned AH_i)
    {
        CH_assert(AH_i < m_apparent_horizons.size());
        return m_apparent_horizons[AH_i];
    }

    //! returns the index of the AH in m_apparent_horizons
    template <class AHInitialGuess>
    int add_ah(const SurfaceGeometry &a_coord_system,
               AHInitialGuess
                   a_initial_guess, //!< Initial guess for radius (or whatever
                                    //!< coordinate you're solving for)
               const AHParams &a_params,    //!< set of AH parameters
               bool solve_first_step = true //!< whether or not to solve if t=0
    );
    //! returns the index of the AH in m_apparent_horizons
    int add_ah(const SurfaceGeometry &a_coord_system,
               const AHInitialGuessPtr
                   &a_initial_guess, //!< Initial guess for radius (or whatever
                                     //!< coordinate you're solving for)
               const AHParams &a_params,    //!< set of AH parameters
               bool solve_first_step = true //!< whether or not to solve if t=0
    );
    //! backward-compatibility: forces AHInitialGuess = AHInitialGuessConstant
    int
    add_ah(const SurfaceGeometry &a_coord_system,
           double a_initial_guess,   //!< Initial guess for radius (or whatever
                                     //!< coordinate you're solving for)
           const AHParams &a_params, //!< set of AH parameters
           bool solve_first_step = true //!< whether or not to solve if t=0
    );

    // returns the index of the AH in m_apparent_horizons
    int add_ah_merger(int ah1, int ah2, const AHParams &a_params);

    //! Find AH; Calculate area and spin; Update center; Print outputs
    void solve(double a_dt, double a_time, double a_restart_time);

    bool need_diagnostics(double a_dt, double a_time)
        const; //!< is any AH printing diagnostics? Useful if Diagnostics need
               //!< to be calculated

    ALWAYS_INLINE void
    set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        m_interpolator = a_interpolator;
    }

    void
    set_origins(const std::vector<std::array<double, CH_SPACEDIM>> &origins,
                bool includes_mergers = false);

  private:
    //! returns false if 'parent' AHs are too far
    //! sets the initial guess and the origin for the merger
    bool solve_merger(int ah1, int ah2, AHInitialGuessPtr &initial_guess_merger,
                      std::array<double, CH_SPACEDIM> &origin_merged);

  private:
    //! if this AH is supposed to track the formation of a merger
    //! this pair indicates the indices of the 2 AHs in the vector
    //! 'm_apparent_horizons'
    std::vector<std::pair<int, int>> m_merger_pairs;
    AMRInterpolator<Lagrange<4>> *m_interpolator; //!< The interpolator pointer

    std::vector<ApparentHorizon<SurfaceGeometry, AHFunction> *>
        m_apparent_horizons; //!< public in case user wants to solve by himself
};

#include "AHFinder.impl.hpp"

#endif /* _AHFINDER_HPP_ */
