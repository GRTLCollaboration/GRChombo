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
 MeansVars::MeansVars(double dx, params_t a_params, int a_slice, std::string a_data_path) : 
    m_dx (dx), m_params (a_params), m_slice (a_slice), m_data_path(a_data_path) {}

 template <class data_t>
 void MeansVars::compute(Cell<data_t> current_cell) const
 {
     CH_TIME("MeansVars::compute");

     fstream field_file;
     field_file.open(m_data_path+"field_step_"+to_string(m_slice)+".dat", std::fstream::app);

     Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

     const auto vars = current_cell.template load_vars<Vars>();

     data_t phisq = vars.phi*vars.phi;
     data_t chisq = vars.chi*vars.chi;
     data_t kin = vars.Pi*vars.Pi;

    //store class (Vars) variables as diagnostic variables on the grid
     current_cell.store_vars(vars.phi, c_sf);
     current_cell.store_vars(vars.Pi, c_sfd);
     current_cell.store_vars(vars.chi, c_a);
     current_cell.store_vars(vars.K, c_H);
     current_cell.store_vars(phisq, c_sf2);
     current_cell.store_vars(chisq, c_ch2);
     current_cell.store_vars(kin, c_kin);

     //current_cell.store_vars(vars.h11, c_h11);

     //field_file << vars.h11 << "\n";

     field_file.close();
 }

 template <class data_t>
 template <typename mapping_function_t>
 void MeansVars::Vars<data_t>::enum_mapping(mapping_function_t mapping_function)
 {
     using namespace VarsTools;
     define_enum_mapping(mapping_function, c_phi, phi);
     define_enum_mapping(mapping_function, c_Pi, Pi);
     define_enum_mapping(mapping_function, c_chi, chi);
     define_enum_mapping(mapping_function, c_K, K);

     define_enum_mapping(mapping_function, c_h11, h11);
 }

 #endif /* MEANSVARS_IMPL_HPP_ */