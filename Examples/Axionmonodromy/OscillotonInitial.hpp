/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef OSCILLOTONINITIAL_HPP_
#define OSCILLOTONINITIAL_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"
#include <fstream>

//! Class which creates Oscillotons (A self gravitating scalar field) using
//! external data-files.
// Information taken form gr-qc/0301105 - Numerical studies of Phi^2 Oscillotons
// The read in method is a bit inefficient but it only gets done once so let's
// not be too fussy
class OscillotonInitial
{

  public:
    OscillotonInitial(double a_L, double a_dx, int sign_of_Pi,
                      std::array<double, CH_SPACEDIM> a_center, double spacing,
                      std::vector<double> lapse_values,
                      std::vector<double> psi_values,
                      std::vector<double> Pi_values)
        : m_L(a_L), m_dx(a_dx), m_spacing(spacing), m_center(a_center),
          m_lapse_values(lapse_values), m_psi_values(psi_values),
          m_Pi_values(Pi_values), m_sign_of_Pi(sign_of_Pi)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    void compute(Cell<double> current_cell) const;

  protected:
    double m_dx;
    double m_L;
    double m_spacing;
    int m_sign_of_Pi;
    std::array<double, CH_SPACEDIM> m_center;
    std::vector<double> m_lapse_values;
    std::vector<double> m_psi_values;
    std::vector<double> m_Pi_values;
};

#include "OscillotonInitial.impl.hpp"

#endif /* OSCILLOTONINITIAL_HPP_ */
