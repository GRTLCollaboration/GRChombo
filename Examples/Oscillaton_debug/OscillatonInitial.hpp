/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef OSCILLATONINITIAL_HPP_
#define OSCILLATONINITIAL_HPP_

#include "Cell.hpp"
#include "CoordinateTransformations.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"
#include <fstream>

//! Class which creates initial conditions for an oscilloton
// (A self gravitating scalar field) using an
// external data-file with radial data
// Information taken form gr-qc/0301105 - Numerical studies of Phi^2 Oscillotons
// The read in method is a bit inefficient but it only gets done once so let's
// not be too fussy
class OscillatonInitial
{
  public:
    OscillatonInitial(const double a_L, const double a_dx,
                      const std::array<double, CH_SPACEDIM> a_center,
                      const double spacing,
                      const std::vector<double> lapse_values,
                      const std::vector<double> grr_values,
                      const std::vector<double> Pi_values)
        : m_L(a_L), m_dx(a_dx), m_spacing(spacing), m_center(a_center),
          m_lapse_values(lapse_values), m_grr_values(grr_values),
          m_Pi_values(Pi_values)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    void compute(Cell<double> current_cell) const
    {
        // may want to add to existing values so load the vars
        auto vars =
            current_cell.template load_vars<MatterCCZ4<ScalarField<>>::Vars>();
        // assign them zero to be sure they are zeroed
        // including at boundaries where required
        VarsTools::assign(vars, 0.);

        // Define Coordinates
        const Coordinates<double> coords(current_cell, m_dx, m_center);
        const double r = coords.get_radius();
        const double x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const double rho = sqrt(x * x + y * y);

        // Interpolate data from read in values
        const int indxL = static_cast<int>(floor(r / m_spacing));
        const int indxH = static_cast<int>(ceil(r / m_spacing));

        // the lapse value
        vars.lapse = m_lapse_values[indxL] +
                     (r / m_spacing - indxL) *
                         (m_lapse_values[indxH] - m_lapse_values[indxL]);

        // the field values
        vars.Pi +=
            m_Pi_values[indxL] +
            (r / m_spacing - indxL) * (m_Pi_values[indxH] - m_Pi_values[indxL]);

        // the metric value dl^2 = grr dr^2 + r^2 dOmega^2
        const double grr = m_grr_values[indxL] +
                           (r / m_spacing - indxL) *
                               (m_grr_values[indxH] - m_grr_values[indxL]);
        Tensor<2, double> spherical_gamma;
        spherical_gamma[0][0] = grr;
        spherical_gamma[1][1] = r * r;
        spherical_gamma[2][2] = rho * rho; // r2 sin2theta
        vars.h = CoordinateTransformations::spherical_to_cartesian_LL(
            spherical_gamma, x, y, z);

        // decompose for BSSN form
        const double det_gamma = TensorAlgebra::compute_determinant(vars.h);
        vars.chi = pow(det_gamma, -1.0 / 3.0);
        FOR2(i, j) { vars.h[i][j] *= vars.chi; }

        // store the values
        current_cell.store_vars(vars);
    }

  protected:
    const double m_dx;
    const double m_L;
    const double m_spacing;
    const std::array<double, CH_SPACEDIM> m_center;
    const std::vector<double> m_lapse_values;
    const std::vector<double> m_grr_values;
    const std::vector<double> m_Pi_values;
};

#endif /* OSCILLOTONINITIAL_HPP_ */
