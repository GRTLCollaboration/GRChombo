/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CatalystAdaptor.hpp"
#include "GRAMRLevel.hpp"

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

void CatalystAdaptor::build_vtk_grid()
{
    if (m_vtk_grid_ptr != nullptr)
        m_vtk_grid_ptr->Delete();

    m_vtk_grid_ptr = vtkOverlappingAMR::New();

    const auto gramrlevels = m_gr_amr_ptr->get_gramrlevels();

    // need the number of levels and the number of boxes on each level to
    // initialise the vtkOverlappingAMR object
    int num_levels = gramrlevels.size();
    int *boxes_per_level = new int[num_levels];
    for (int ilevel = 0; ilevel < num_levels; ++ilevel)
    {
        const GRAMRLevel *level = gramrlevels[ilevel];
        const DisjointBoxLayout &level_box_layout =
            level->getLevelData().disjointBoxLayout();
        boxes_per_level[ilevel] = level_box_layout.boxArray().size();
    }
    m_vtk_grid_ptr->Initialize(num_levels, boxes_per_level);

    // The origin is always at (0, 0, 0) in Chombo
    double origin_global[3] = {0., 0., 0.};
    m_vtk_grid_ptr->SetOrigin(origin_global);

    // now add all the boxes
    for (int ilevel = 0; ilevel < num_levels; ++ilevel)
    {
        const GRAMRLevel *level = gramrlevels[ilevel];
        const double dx = level->get_dx();
        const double dx_arr[3] = {dx, dx, dx};
        m_vtk_grid_ptr->SetSpacing(ilevel, dx_arr);
        m_vtk_grid_ptr->SetRefinementRatio(ilevel, level->refRatio());
        const GRLevelData &level_data = level->getLevelData();
        const DisjointBoxLayout &level_box_layout =
            level_data.disjointBoxLayout();
        LayoutIterator lit = level_box_layout.layoutIterator();
        int ibox;
        for (ibox = 0, lit.begin(); lit.ok(); ++lit, ++ibox)
        {
            // first get the box without ghosts
            Box box = level_box_layout[lit];

            // now grow the box to make the ghosted box
            Box ghosted_box = box;
            ghosted_box.grow(level_data.ghostVect());
            const IntVect &small_end = ghosted_box.smallEnd();
            const IntVect &big_end = ghosted_box.bigEnd();

            vtkAMRBox vtk_amr_box(small_end.dataPtr(), big_end.dataPtr());
            m_vtk_grid_ptr->SetAMRBox(ilevel, ibox, vtk_amr_box);

            bool local_box = (procID() == level_box_layout.procID(lit()));
            // only need to do the following for local boxes
            if (local_box)
            {
                vtkNew<vtkUniformGrid> vtk_uniform_grid_ptr;
                vtk_uniform_grid_ptr->SetOrigin(origin_global);
                vtk_uniform_grid_ptr->SetSpacing(dx_arr);
                vtk_uniform_grid_ptr->SetExtent(small_end[0], big_end[0] + 1,
                                                small_end[1], big_end[1] + 1,
                                                small_end[2], big_end[2] + 1);
                m_vtk_grid_ptr->SetDataSet(ilevel, ibox, vtk_uniform_grid_ptr);
            }
        }
    }

    // not sure if this is necessary but it was on the OverlappingAMR example
    m_vtk_grid_ptr->GenerateParentChildInformation();
}

void CatalystAdaptor::add_var() {}

void CatalystAdaptor::coprocess(double time) {}

#endif /* USE_CATALYST */
