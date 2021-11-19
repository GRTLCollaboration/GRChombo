/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CatalystAdaptor.hpp"

#ifdef USE_CATALYST

CatalystAdaptor::CatalystAdaptor() {}

CatalystAdaptor::CatalystAdaptor(const GRAMR *a_gr_amr_ptr,
                                 std::string a_python_script_path)
{
    initialise(a_gr_amr_ptr, a_python_script_path);
}

CatalystAdaptor::~CatalystAdaptor() {}

void CatalystAdaptor::initialise(const GRAMR *a_gr_amr_ptr,
                                 std::string a_python_script_path)
{
    if (!a_gr_amr_ptr)
    {
        std::cerr << "CatalystAdaptor::initialise: failed to initalise due to"
                     " invalid GRAMR pointer"
                  << std::endl;
        m_initialised = false;
        return;
    }
    m_gr_amr_ptr = a_gr_amr_ptr;
}
void CatalystAdaptor::finalise();

void CatalystAdaptor::update_vtk_grid();

void CatalystAdaptor::add_var();

void CatalystAdaptor::coprocess(double time);

#endif /* USE_CATALYST */
