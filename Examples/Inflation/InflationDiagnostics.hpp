/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INFLATIONDIAGNOSTICS_HPP_
#define INFLATIONDIAGNOSTICS_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DiagnosticVariables.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "simd.hpp"
#include <array>
#include "CCZ4Vars.hpp"
#include "NewMatterConstraints.hpp"
#include "FourthOrderDerivatives.hpp"


//! Calculates all relevant variables, which
//! are stored as diagnostics

template <class matter_t> 
class InflationDiagnostics //public MatterConstraints<matter_t>
{
  public:
    // Inherit the variable definitions from CCZ4 + matter_t
    template <class data_t>
    using CCZ4Vars = typename CCZ4::template Vars<data_t>;

    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    template <class data_t>
    struct Vars : public CCZ4Vars<data_t>, public MatterVars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            CCZ4Vars<data_t>::enum_mapping(mapping_function);
            MatterVars<data_t>::enum_mapping(mapping_function);
        }
    };

    //! Constructor of class InflationDiagnostics
    InflationDiagnostics(const matter_t a_matter
                      , double dx, double G_Newton)
                      : m_matter(a_matter), m_deriv(dx)
                      {
                      }

    //! The compute member which calculates the constraints at each point in the
    //! box
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    matter_t m_matter; //!< The matter object, e.g. a scalar field
    const FourthOrderDerivatives m_deriv;

};

#include "InflationDiagnostics.impl.hpp"

#endif /* INFLATIONDIAGNOSTICS_HPP_ */
