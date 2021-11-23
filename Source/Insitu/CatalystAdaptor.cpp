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

CatalystAdaptor::~CatalystAdaptor()
{
    if (m_initialised)
    {
        finalise();
    }
}

void CatalystAdaptor::initialise(const GRAMR *a_gr_amr_ptr,
                                 std::string a_python_script_path)
{
    // don't initalise twice
    if (m_initialised)
        return;

    if (!a_gr_amr_ptr)
    {
        std::cerr << "CatalystAdaptor::initialise: failed to initalise due to"
                     " invalid GRAMR pointer"
                  << std::endl;
        m_initialised = false;
        return;
    }
    m_gr_amr_ptr = a_gr_amr_ptr;

    // Initialise VTK CP Processor
    if (!m_proc_ptr)
    {
        m_proc_ptr = vtkCPProcessor::New();
        m_proc_ptr->Initialize();
    }
    else
    {
        m_proc_ptr->RemoveAllPipelines();
    }

    // Create Python script pipeline and add it to the VTK CP Processor
    if (auto pipeline = vtkCPPythonScriptPipeline::CreateAndInitializePipeline(
            a_python_script_path.c_str()))
    {
        m_proc_ptr->AddPipeline(pipeline);
    }
    else
    {
        std::string pipeline_fail_warning =
            "Failed to set up pipeline for script:";
        pipeline_fail_warning += a_python_script_path;
        MayDay::Warning(pipeline_fail_warning.c_str());
    }

    m_initialised = true;
}

void CatalystAdaptor::finalise()
{
    if (m_proc_ptr)
    {
        m_proc_ptr->Delete();
        m_proc_ptr = nullptr;
    }
    if (m_vtk_grid_ptr)
    {
        m_vtk_grid_ptr->Delete();
        m_vtk_grid_ptr = nullptr;
    }
    m_initialised = false;
}

void CatalystAdaptor::update_vtk_grid() {}

void CatalystAdaptor::add_var() {}

void CatalystAdaptor::coprocess(double time) {}

#endif /* USE_CATALYST */
