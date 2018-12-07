/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CCZ4_HPP_
#define CCZ4_HPP_

#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"

#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

#include <array>

/// Compute class to calculate the CCZ4 right hand side
/**
 * This compute class implements the CCZ4 right hand side equations. Use it by
 *handing it to a loop in the BoxLoops namespace. CCZ4 includes two classes in
 *its scope: CCZ4::Vars (the CCZ4 variables like conformal factor, conformal
 *metric, extrinsic curvature, etc) and CCZ4::Params (parameters necessary for
 *CCZ4 like gauge and damping parameters).
 **/
class CCZ4
{
  public:
    enum
    {
        USE_CCZ4,
        USE_BSSN
    };

    /// CCZ4 variables
    template <class data_t> using Vars = CCZ4Vars::VarsWithGauge<data_t>;

    /// CCZ4 variables
    template <class data_t>
    using Diff2Vars = CCZ4Vars::Diff2VarsWithGauge<data_t>;

    /// Parameters for CCZ4
    /** This struct collects all parameters that are necessary for CCZ4 such as
     * gauge and damping parameters.  */
    struct params_t
    {
        double kappa1; //!< Damping parameter kappa1 as in arXiv:1106.2254
        double kappa2; //!< Damping parameter kappa2 as in arXiv:1106.2254
        double kappa3; //!< Damping parameter kappa3 as in arXiv:1106.2254
        double shift_Gamma_coeff = 0.75; //!< Gives the F in \f$\partial_t
                                         //!  \beta^i =  F B^i\f$
        double lapse_advec_coeff = 0.;   //!< Switches advection terms in
                                         //! the lapse condition on/off
        double shift_advec_coeff = 0.;   //!< Switches advection terms in the
                                         //! shift condition on/off
        double eta = 1.; //!< The eta in \f$\partial_t B^i = \partial_t \tilde
                         //!\Gamma - \eta B^i\f$
        double lapse_power = 1.; //!< The power p in \f$\partial_t \alpha = - c
                                 //!\alpha^p(K-2\Theta)\f$
        double lapse_coeff = 2.; //!< The coefficient c in \f$\partial_t \alpha
                                 //!= -c \alpha^p(K-2\Theta)\f$
    };

  protected:
    const params_t m_params; //!< CCZ4 parameters
    const double m_sigma;    //!< Coefficient for Kreiss-Oliger dissipation
    int m_formulation;
    double m_cosmological_constant;
    const FourthOrderDerivatives m_deriv;

  public:
    /// Constructor
    CCZ4(params_t params,            //!< The CCZ4 parameters
         double dx,                  //!< The grid spacing
         double sigma,               //!< Kreiss-Oliger dissipation coefficient
         int formulation = USE_CCZ4, //!< Switches between CCZ4, BSSN,...
         double cosmological_constant = 0 //!< Value of the cosmological const.
         );

    /// Compute function
    /** This function orchestrates the calculation of the rhs for one specific
     * grid cell. This function is called by the BoxLoops::loop for each grid
     * cell; there should rarely be a need to call it directly.
     */
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    /// Calculates the rhs for CCZ4
    /** Calculates the right hand side for CCZ4 with slicing \f$- n \alpha^m (K
     *- 2\Theta)\f$ and Gamma-Driver shift condition. The variables (the
     *template argument vars_t) must contain at least the members: chi,
     *h[i][j], Gamma[i], A[i][j], Theta, lapse and shift[i].
     **/
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    void rhs_equation(
        vars_t<data_t> &rhs, //!< Reference to the variables into which the
                             //! output right hand side is written
        const vars_t<data_t> &vars, //!< The values of the current variables
        const vars_t<Tensor<1, data_t>>
            &d1, //!< First derivative of the variables
        const diff2_vars_t<Tensor<2, data_t>>
            &d2, //!< The second derivative the variables
        const vars_t<data_t>
            &advec //!< The advection derivatives of the variables
        ) const;
};

#include "CCZ4.impl.hpp"

#endif /* CCZ4_HPP_ */
