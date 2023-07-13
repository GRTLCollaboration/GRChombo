/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 #if !defined(MEANSVARS_HPP_)
 #error "Only include this file through MeansVars.hpp"
 #endif

 #ifndef MEANSVARS_IMPL_HPP
 #define MEANSVARS_IMPL_HPP

#include "Coordinates.hpp"
#include "VarsTools.hpp"
#include <cmath>

inline
 MeansVars::MeansVars(double dx, params_t a_params) : 
    m_dx (dx), m_params (a_params) {}

 template <class data_t>
 void MeansVars::compute(Cell<data_t> current_cell) const
 {
     CH_TIME("MeansVars::compute");

     Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

     const auto vars = current_cell.template load_vars<Vars>();

     data_t phisq = vars.phi*vars.phi;
     data_t kin = 0.5*vars.Pi*vars.Pi;

    //store class (Vars) variables as diagnostic variables on the grid
     current_cell.store_vars(vars.phi, c_sf);
     current_cell.store_vars(vars.Pi, c_sfd);
     current_cell.store_vars(vars.chi, c_a);
     current_cell.store_vars(vars.K, c_H);
     current_cell.store_vars(phisq, c_sf2);
     current_cell.store_vars(kin, c_kin);
 }

 template <class data_t>
 template <typename mapping_function_t>
 void MeansVars::Vars<data_t>::enum_mapping(mapping_function_t mapping_function)
 {
     using namespace VarsTools;
     define_enum_mapping(mapping_function, c_sf, phi);
     define_enum_mapping(mapping_function, c_sfd, Pi);
     define_enum_mapping(mapping_function, c_chi, chi);
     define_enum_mapping(mapping_function, c_K, K);
 }

 #endif /* MEANSVARS_IMPL_HPP_ */