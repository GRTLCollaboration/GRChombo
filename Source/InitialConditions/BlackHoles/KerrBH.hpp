/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef KERRBH_HPP_
#define KERRBH_HPP_

#include "ADMConformalVars.hpp"
#include "Cell.hpp"
#include "CoordinateTransformations.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which computes the Kerr initial conditions per arXiv 1401.1548
class KerrBH
{
    // Use the variable definition in CCZ4
    template <class data_t>
    using Vars = ADMConformalVars::VarsWithGauge<data_t>;

  public:
    //! Stuct for the params of the Kerr BH
    struct params_t
    {
        double mass;                            //!<< The mass of the Kerr BH
        std::array<double, CH_SPACEDIM> center; //!< The center of the Kerr BH
        double spin; //!< The spin param a = J/M, so 0 <= |a| <= M
        std::array<double, CH_SPACEDIM> spin_direction = {
            0., 0., 1.}; // default to 'z' axis; doesn't need to be normalized
    };

  protected:
    double m_dx;
    params_t m_params;

  public:
    KerrBH(params_t a_params, double a_dx) : m_dx(a_dx), m_params(a_params)

    {
        // check this spin param is sensible
        if (std::abs(m_params.spin) > m_params.mass)
        {
            MayDay::Error("The spin parameter must satisfy |a| <= M");
        }
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    //! Function which computes the components of the metric in spherical coords
    template <class data_t>
    void compute_kerr(
        Tensor<2, data_t>
            &spherical_g, //!<< The spatial metric in spherical coords
        Tensor<2, data_t>
            &spherical_K, //!<< The extrinsic curvature in spherical coords
        Tensor<1, data_t>
            &spherical_shift, //!<< The spherical components of the shift
        data_t &kerr_lapse,   //!<< The lapse for the kerr solution
        const Tensor<1, data_t> &coords //!<< Coords of current cell
    ) const;
};

#include "KerrBH.impl.hpp"

#endif /* KERRBH_HPP_ */
