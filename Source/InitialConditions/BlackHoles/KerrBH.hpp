/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef KERRBH_HPP_
#define KERRBH_HPP_

#include "ADMConformalVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "InitialDataTools.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

enum Lapse
{
    ONE,
    PRE_COLLAPSED,
    ANALYTIC
};

enum Kerr
{
    QUASI_ISOTROPIC,
    KERR_SCHILD,
    BOWEN_YORK
};

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
        double mass; //!<< The bare mass of the Kerr BH m_bare
        std::array<double, CH_SPACEDIM> center; //!< The center of the Kerr BH
        double spin; //!< The spin param a (NB for Bowen York this will set J_z)
        std::array<double, CH_SPACEDIM> boost = {
            0.0, 0.0, 0.0}; //!< Boost is only for Bowen York data
    };

  protected:
    const double m_dx;
    const params_t m_params;
    const int m_initial_lapse;
    const int m_kerr_solution;

  public:
    KerrBH(params_t a_params, double a_dx,
           int a_initial_lapse = Lapse::PRE_COLLAPSED,
           int a_kerr_solution = Kerr::QUASI_ISOTROPIC)
        : m_dx(a_dx), m_params(a_params), m_initial_lapse(a_initial_lapse),
          m_kerr_solution(a_kerr_solution)

    {
        // check this spin param is sensible
        if ((m_params.spin > 1.0) || (m_params.spin < 0.0))
        {
            MayDay::Error(
                "The spin parameter must be in the range 0 < a < 1.0");
        }
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    /// Quasi isotropic solution for high spin BHs per 1001.4077
    template <class data_t>
    void quasi_isotropic_kerr(Vars<data_t> &vars,
                              const Coordinates<data_t> &coords) const;

    /// Bowen York solution with puncture
    /// per http://www.livingreviews.org/Articles/Volume3/2000-5cook
    template <class data_t>
    void bowen_york(Vars<data_t> &vars,
                    const Coordinates<data_t> &coords) const;

    /// Kerr Schild solution, must be used with excision at centre
    /// vars per https://arxiv.org/abs/gr-qc/0109032
    template <class data_t>
    void kerr_schild(Vars<data_t> &vars,
                     const Coordinates<data_t> &coords) const;

    /// Work out the gradients of the quantities H and el appearing in the Kerr
    /// Schild solution
    template <class data_t>
    void get_KS_derivs(Tensor<1, data_t> &dHdx, Tensor<2, data_t> &dldx,
                       const data_t &r_BL, const data_t &H,
                       const Coordinates<data_t> &coords) const;
};

#include "KerrBH.impl.hpp"

#endif /* KERRBH_HPP_ */
