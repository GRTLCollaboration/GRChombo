/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CatalystAdaptor.hpp"
#include "GRAMRLevel.hpp"

#ifdef USE_CATALYST

CatalystAdaptor::CatalystAdaptor() {}

CatalystAdaptor::CatalystAdaptor(
    GRAMR *a_gr_amr_ptr, const std::vector<std::string> &a_python_scripts,
    const std::vector<std::pair<int, VariableType>> &a_vars,
    bool a_abort_on_catalyst_error, int a_verbosity)
{
    initialise(a_gr_amr_ptr, a_python_scripts, a_vars,
               a_abort_on_catalyst_error, a_verbosity);
}

CatalystAdaptor::~CatalystAdaptor()
{
    if (m_initialised)
    {
        finalise();
    }
}

void CatalystAdaptor::initialise(
    GRAMR *a_gr_amr_ptr, const std::vector<std::string> &a_python_scripts,
    const std::vector<std::pair<int, VariableType>> &a_vars,
    bool a_abort_on_catalyst_error, int a_verbosity)
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
    m_vars = a_vars;
    m_abort_on_catalyst_error = a_abort_on_catalyst_error;

    // Initialise VTK CP Processor
    if (!m_proc_ptr)
    {
        m_proc_ptr = vtkCPProcessor::New();
        int m_proc_ptr_intialization_success = m_proc_ptr->Initialize();
        std::string error_msg = "Failed to initialize vtkCPProcessor in "
                                "CatalystAdaptor::initialise";
        catalyst_error_or_warning(m_proc_ptr_intialization_success, error_msg);
    }
    else
    {
        m_proc_ptr->RemoveAllPipelines();
    }

    // Create Python script pipeline and add it to the VTK CP Processor
    for (const std::string &script : a_python_scripts)
    {
        if (auto pipeline =
                vtkCPPythonScriptPipeline::CreateAndInitializePipeline(
                    script.c_str()))
        {
            int add_pipeline_success = m_proc_ptr->AddPipeline(pipeline);
            std::string add_pipeline_fail_msg =
                "Failed to add pipeline for script: ";
            add_pipeline_fail_msg += script;
            catalyst_error_or_warning(add_pipeline_success,
                                      add_pipeline_fail_msg);
        }
        else
        {
            std::string pipeline_fail_msg =
                "Failed to set up pipeline for script: ";
            pipeline_fail_msg += script;
            catalyst_error_or_warning(false, pipeline_fail_msg);
        }
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
    const IntVect &coarsest_ghost_vect =
        gramrlevels[0]->getLevelData().ghostVect();
    const double coarsest_dx = gramrlevels[0]->get_dx();
    RealVect ghosted_origin_global_vect = coarsest_ghost_vect;
    ghosted_origin_global_vect *= -coarsest_dx;
    m_vtk_grid_ptr->SetOrigin(ghosted_origin_global_vect.dataPtr());

    // now add all the boxes
    for (int ilevel = 0; ilevel < num_levels; ++ilevel)
    {
        // pout() << "========================================\n";
        // pout() << "Level: " << ilevel << std::endl;
        const GRAMRLevel *level = gramrlevels[ilevel];
        const double dx = level->get_dx();
        double dx_arr[3] = {dx, dx, dx};
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
            // RealVect level_origin = -dx * level_data.ghostVect();

            vtkAMRBox vtk_amr_box(origin, ghosted_box.size().dataPtr(), dx_arr,
                                  ghosted_origin_global_vect.dataPtr());
            // vtk_amr_box.Print(pout());
            // pout() << "\n";
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
                // int no_ghost[6] = {small_end[0], big_end[0] + 1,
                //    small_end[1], big_end[1] + 1,
                //    small_end[2], big_end[2] + 1};
                // bool cell_data = true;
                // vtk_uniform_grid_ptr->GenerateGhostArray(no_ghost,
                // cell_data);

                // vtk_uniform_grid_ptr->Initialize(
                // &vtk_amr_box, origin_global, dx_arr,
                // level_data.ghostVect().dataPtr());

                m_vtk_grid_ptr->SetDataSet(ilevel, ibox, vtk_uniform_grid_ptr);
            }
        }
    }

    m_vtk_grid_ptr->Audit();
#if DEBUG
    const double *vtk_grid_bounds = m_vtk_grid_ptr->GetAMRInfo()->GetBounds();
    pout() << "VTK Grid Bounds:"
           << "(" << vtk_grid_bounds[0] << "," << vtk_grid_bounds[2] << ","
           << vtk_grid_bounds[4] << "), (" << vtk_grid_bounds[1] << ","
           << vtk_grid_bounds[3] << "," << vtk_grid_bounds[5] << ")"
           << std::endl;
#endif

    // not sure if this is necessary but it was on the OverlappingAMR example
    m_vtk_grid_ptr->GenerateParentChildInformation();
}

void CatalystAdaptor::add_vars(vtkCPInputDataDescription *a_input_data_desc)
{
    if (m_verbosity)
    {
        pout() << "CatalystAdaptor::add_vars" << std::endl;
    }

    std::array<bool, NUM_VARS> requested_evolution_vars;
    std::array<bool, NUM_DIAGNOSTIC_VARS> requested_diagnostic_vars;

    if (m_verbosity)
    {
        pout() << "CatalystAdaptor Requested variables: ";
    }

    for (int ivar = 0; ivar < NUM_VARS; ++ivar)
    {
        requested_evolution_vars[ivar] = a_input_data_desc->IsFieldNeeded(
            UserVariables::variable_names[ivar].c_str(), vtkDataObject::CELL);
        if (m_vars.size() > 0)
        {
            bool pass_var =
                !(std::find(m_vars.begin(), m_vars.end(),
                            std::make_pair(ivar, VariableType::evolution)) ==
                  m_vars.end());
            requested_evolution_vars[ivar] &= pass_var;
        }
        if (m_verbosity && requested_evolution_vars[ivar])
            pout() << UserVariables::variable_names[ivar] << " ";
    }
    for (int ivar = 0; ivar < NUM_DIAGNOSTIC_VARS; ++ivar)
    {
        requested_diagnostic_vars[ivar] = a_input_data_desc->IsFieldNeeded(
            DiagnosticVariables::variable_names[ivar].c_str(),
            vtkDataObject::CELL);
        if (m_vars.size() > 0)
        {
            bool pass_var =
                !(std::find(m_vars.begin(), m_vars.end(),
                            std::make_pair(ivar, VariableType::diagnostic)) ==
                  m_vars.end());
            requested_diagnostic_vars[ivar] &= pass_var;
        }
        if (m_verbosity && requested_diagnostic_vars[ivar])
            pout() << DiagnosticVariables::variable_names[ivar] << " ";
    }
    if (m_verbosity)
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
        GRLevelData &diagnostic_level_data =
            (NUM_DIAGNOSTIC_VARS > 0)
                ? const_cast<GRLevelData &>(
                      level->getLevelData(VariableType::diagnostic))
                : evolution_level_data;

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
                FArrayBox &diagnostic_fab = (NUM_DIAGNOSTIC_VARS > 0)
                                                ? diagnostic_level_data[dind]
                                                : evolution_level_data[dind];

#if DEBUG
                vtkAMRBox vtk_box = m_vtk_grid_ptr->GetAMRBox(ilevel, ibox);
                // shift to account for different indexing
                IntVect vtk_origin_offset =
                    -ipow(2, ilevel) * evolution_level_data.ghostVect();
                vtk_box.Shift(vtk_origin_offset.dataPtr());
                const int *vtk_box_lo = vtk_box.GetLoCorner();
                const int *vtk_box_hi = vtk_box.GetHiCorner();

                const Box &evolution_box = evolution_fab.box();
                const IntVect &evolution_box_lo = evolution_box.smallEnd();
                const IntVect &evolution_box_hi = evolution_box.bigEnd();

                const Box &diagnostic_box = diagnostic_fab.box();
                const IntVect &diagnostic_box_lo = diagnostic_box.smallEnd();
                const IntVect &diagnostic_box_hi = diagnostic_box.bigEnd();
                // add 1 to hi corner of VTK box as this is what is meant by hi
                // corner in VTK
                bool all_same_boxes =
                    (vtk_box_lo[0] == evolution_box_lo[0]) &&
                    (vtk_box_lo[1] == evolution_box_lo[1]) &&
                    (vtk_box_lo[2] == evolution_box_lo[2]) &&
                    (vtk_box_hi[0] + 1 == evolution_box_hi[0]) &&
                    (vtk_box_hi[1] + 1 == evolution_box_hi[1]) &&
                    (vtk_box_hi[2] + 1 == evolution_box_hi[2]) &&
                    (vtk_box_lo[0] == diagnostic_box_lo[0]) &&
                    (vtk_box_lo[1] == diagnostic_box_lo[1]) &&
                    (vtk_box_lo[2] == diagnostic_box_lo[2]) &&
                    (vtk_box_hi[0] + 1 == diagnostic_box_hi[0]) &&
                    (vtk_box_hi[1] + 1 == diagnostic_box_hi[1]) &&
                    (vtk_box_hi[2] + 1 == diagnostic_box_hi[2]);
                if (!all_same_boxes)
                {
                    pout() << "Boxes do not agree: \n";
                    pout() << "vtk_box: ";
                    vtk_box.Print(pout());
                    pout() << "\n";
                    pout() << "evolution_box: " << evolution_box << "\n";
                    pout() << "diagnostic_box: " << diagnostic_box << std::endl;
                }
#endif

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
        // vtkNew<vtkOverlappingAMR> stripped_vtk_grid;
        // vtkAMRUtilities::StripGhostLayers(m_vtk_grid_ptr, stripped_vtk_grid);
        input_data_description->SetGrid(m_vtk_grid_ptr);
        int coprocess_success = m_proc_ptr->CoProcess(data_description);
        catalyst_error_or_warning(coprocess_success,
                                  "Error in Catalyst CoProcess");
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

void CatalystAdaptor::catalyst_error_or_warning(bool a_success,
                                                std::string a_msg)
{
    if (a_success)
        return;

    if (m_abort_on_catalyst_error)
        MayDay::Error(a_msg.c_str());
    else
        MayDay::Warning(a_msg.c_str());
}

#endif /* USE_CATALYST */