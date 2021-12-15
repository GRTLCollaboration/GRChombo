/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CatalystAdaptor.hpp"
#include "GRAMRLevel.hpp"

#ifdef USE_CATALYST

CatalystAdaptor::CatalystAdaptor() {}

CatalystAdaptor::CatalystAdaptor(GRAMR *a_gr_amr_ptr,
                                 std::string a_python_script_path,
                                 int a_verbosity)
{
    initialise(a_gr_amr_ptr, a_python_script_path, a_verbosity);
}

CatalystAdaptor::~CatalystAdaptor()
{
    if (m_initialised)
    {
        finalise();
    }
}

void CatalystAdaptor::initialise(GRAMR *a_gr_amr_ptr,
                                 std::string a_python_script_path,
                                 int a_verbosity)
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

    m_verbosity = a_verbosity;

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
    if (m_verbosity)
    {
        pout() << "CatalystAdaptor::build_vtk_grid" << std::endl;
    }
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
            const IntVect &small_end = box.smallEnd();
            const IntVect &big_end = box.bigEnd();

            // now grow the box to make the ghosted box
            Box ghosted_box = box;
            ghosted_box.grow(level_data.ghostVect());
            const IntVect &small_ghosted_end = ghosted_box.smallEnd();
            const IntVect &big_ghosted_end = ghosted_box.bigEnd();

            // vtkAMRBox vtk_amr_box(small_ghosted_end.dataPtr(),
            //                       big_ghosted_end.dataPtr());
            double origin[3] = {dx_arr[0] * small_ghosted_end[0],
                                dx_arr[1] * small_ghosted_end[1],
                                dx_arr[2] * small_ghosted_end[2]};

            vtkAMRBox vtk_amr_box(origin, ghosted_box.size().dataPtr(), dx_arr,
                                  origin_global);
            m_vtk_grid_ptr->SetAMRBox(ilevel, ibox, vtk_amr_box);

            bool local_box = (procID() == level_box_layout.procID(lit()));
            // only need to do the following for local boxes
            if (local_box)
            {
                vtkNew<vtkUniformGrid> vtk_uniform_grid_ptr;
                vtk_uniform_grid_ptr->SetOrigin(origin_global);
                vtk_uniform_grid_ptr->SetSpacing(dx_arr);
                vtk_uniform_grid_ptr->SetExtent(
                    small_ghosted_end[0], big_ghosted_end[0] + 1,
                    small_ghosted_end[1], big_ghosted_end[1] + 1,
                    small_ghosted_end[2], big_ghosted_end[2] + 1);
                // add the ghost cell information
                int no_ghost[6] = {small_end[0], big_end[0] + 1,
                                   small_end[1], big_end[1] + 1,
                                   small_end[2], big_end[2] + 1};
                bool cell_data = true;
                vtk_uniform_grid_ptr->GenerateGhostArray(no_ghost, cell_data);
                m_vtk_grid_ptr->SetDataSet(ilevel, ibox, vtk_uniform_grid_ptr);
            }
        }
    }

    m_vtk_grid_ptr->Audit();
    // not sure if this is necessary but it was on the OverlappingAMR example
    // m_vtk_grid_ptr->GenerateParentChildInformation();
}

void CatalystAdaptor::add_vars(vtkCPInputDataDescription *a_input_data_desc)
{
    if (m_verbosity)
    {
        pout() << "CatalystAdaptor::add_vars" << std::endl;
    }

    std::array<bool, NUM_VARS> requested_evolution_vars;
    std::array<bool, NUM_DIAGNOSTIC_VARS> requested_diagnostic_vars;

    if (m_verbosity > 1)
    {
        pout() << "CatalystAdaptor Requested variables:\n";
    }

    for (int ivar = 0; ivar < NUM_VARS; ++ivar)
    {
        requested_evolution_vars[ivar] = a_input_data_desc->IsFieldNeeded(
            UserVariables::variable_names[ivar].c_str(), vtkDataObject::CELL);
        if (m_verbosity > 1 && requested_evolution_vars[ivar])
            pout() << UserVariables::variable_names[ivar] << " ";
    }
    for (int ivar = 0; ivar < NUM_DIAGNOSTIC_VARS; ++ivar)
    {
        requested_diagnostic_vars[ivar] = a_input_data_desc->IsFieldNeeded(
            DiagnosticVariables::variable_names[ivar].c_str(),
            vtkDataObject::CELL);
        if (m_verbosity > 1 && requested_diagnostic_vars[ivar])
            pout() << DiagnosticVariables::variable_names[ivar] << " ";
    }
    if (m_verbosity > 1)
        pout() << std::endl;

    vtkAMRInformation *amr_info = m_vtk_grid_ptr->GetAMRInfo();
    auto gramrlevels = m_gr_amr_ptr->get_gramrlevels();

    for (int ilevel = 0; ilevel < gramrlevels.size(); ++ilevel)
    {
        GRAMRLevel *level = gramrlevels[ilevel];
        // Unfortunately it doesn't seem possible to pass const pointers
        // to Catalyst
        GRLevelData &evolution_level_data = const_cast<GRLevelData &>(
            level->getLevelData(VariableType::evolution));
        GRLevelData &diagnostic_level_data = const_cast<GRLevelData &>(
            level->getLevelData(VariableType::diagnostic));

        const DisjointBoxLayout &level_box_layout =
            evolution_level_data.disjointBoxLayout();
        LayoutIterator lit = level_box_layout.layoutIterator();
        int ibox;
        for (ibox = 0, lit.begin(); lit.ok(); ++lit, ++ibox)
        {
            // only add data that we have locally
            bool local_box = (procID() == level_box_layout.procID(lit()));
            if (local_box)
            {
                vtkUniformGrid *vtk_uniform_grid_ptr =
                    m_vtk_grid_ptr->GetDataSet(ilevel, ibox);
                // hopefully this promotion works
                DataIndex dind(lit());
                FArrayBox &evolution_fab = evolution_level_data[dind];
                FArrayBox &diagnostic_fab = diagnostic_level_data[dind];

                for (int ivar = 0; ivar < NUM_VARS; ++ivar)
                {
                    if (requested_evolution_vars[ivar])
                    {
                        vtkDoubleArray *vtk_double_arr = fab_to_vtk_array(
                            evolution_fab, ivar,
                            UserVariables::variable_names[ivar]);
                        vtk_uniform_grid_ptr->GetCellData()->AddArray(
                            vtk_double_arr);
                    }
                }
                for (int ivar = 0; ivar < NUM_DIAGNOSTIC_VARS; ++ivar)
                {
                    if (requested_diagnostic_vars[ivar])
                    {
                        vtkDoubleArray *vtk_double_arr = fab_to_vtk_array(
                            diagnostic_fab, ivar,
                            DiagnosticVariables::variable_names[ivar]);
                        vtk_uniform_grid_ptr->GetCellData()->AddArray(
                            vtk_double_arr);
                    }
                }
            }
        }
    }
}

void CatalystAdaptor::coprocess(double a_time, unsigned int a_timestep)
{
    pout() << "CatalystAdaptor::coprocess at time " << a_time << " and step "
           << a_timestep << std::endl;

    vtkNew<vtkCPDataDescription> data_description;
    data_description->AddInput("input");
    data_description->SetTimeData(a_time, a_timestep);

    if (m_proc_ptr->RequestDataDescription(data_description) != 0)
    {
        build_vtk_grid();
        auto input_data_description =
            data_description->GetInputDescriptionByName("input");
        add_vars(input_data_description);
        input_data_description->SetGrid(m_vtk_grid_ptr);
        m_proc_ptr->CoProcess(data_description);
    }
}

vtkDoubleArray *CatalystAdaptor::fab_to_vtk_array(FArrayBox &a_fab, int a_var,
                                                  const std::string &a_name)
{
    vtkDoubleArray *out = vtkDoubleArray::New();
    vtkIdType num_cells = a_fab.size().product();
    out->SetNumberOfTuples(num_cells);
    out->SetName(a_name.c_str());
    // this prevents Catalyst from deallocating the Chombo
    // data pointers
    int save_data = 1;
    out->SetArray(a_fab.dataPtr(a_var), num_cells, save_data);
    return out;
}

#endif /* USE_CATALYST */
