/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CATALYSTADAPTOR_HPP_
#define CATALYSTADAPTOR_HPP_

#ifdef USE_CATALYST

// Standard library includes
#include <string>

// Chombo includes
#include "Box.H"
#include "FArrayBox.H"
#include "MayDay.H"

// GRChombo includes
#include "UserVariables.hpp"

// ParaView/VTK/Catalyst includes
#include <vtkAMRBox.h>
#include <vtkAMRInformation.h>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkNew.h>
#include <vtkOverlappingAMR.h>
#include <vtkUniformGrid.h>
#include <vtkXMLPUniformGridAMRWriter.h>

// Chombo namespace
#include "UsingNamespace.H"

// Forward declaration of GRAMR class
class GRAMR;

// Class that interfaces with ParaView Catalyst for insitu visualisation
class CatalystAdaptor
{
  public:
    struct params_t
    {
        int verbosity;
        bool abort_on_error = false;
        bool remove_ghosts = false;

        // Pipeline parameters
        std::vector<std::string> python_scripts;
        std::string output_path;

        // VTK file writing parameters
        bool write_vtk_files = false;
        std::string vtk_file_prefix = "Catalyst_VTK_grid_";

        // variables to pass to Catalyst set by GRChombo parameter
        std::vector<std::pair<int, VariableType>> vars;
    };

    // empty constructor (doesn't call initialise)
    CatalystAdaptor();

    // full constructor (calls initialise)
    CatalystAdaptor(GRAMR *a_gr_amr_ptr, const params_t &a_params);

    // destructor
    ~CatalystAdaptor();

    // Initialisation/Finalisation
    void initialise(GRAMR *m_gr_amr_ptr, const params_t &a_params);
    void finalise();

    // do Catalyst processing
    void coprocess(double a_time, unsigned int a_timestep);

    // returns the variables that were last requested/sent to Catalyst
    const std::array<bool, NUM_VARS> &get_requested_evolution_vars();
    const std::array<bool, NUM_DIAGNOSTIC_VARS> &
    get_requested_diagnostic_vars();

  private:
    // update the AMR grid (no grid data)
    void build_vtk_grid();

    // send variables to catalyst
    void add_vars(vtkCPInputDataDescription *a_input_data_desc);

    // write VTK grid to a file
    void write_vtk_grid(unsigned int a_timestep);

    // directly passes the FAB pointer to the VTK array
    vtkDoubleArray *fab_to_vtk_array(FArrayBox &a_fab, int a_var,
                                     const std::string &a_name);

    // copies the data in the FAB to the VTK array without the ghosts
    vtkDoubleArray *fab_to_vtk_array_without_ghosts(FArrayBox &a_fab,
                                                    const Box &a_unghosted_box,
                                                    int a_var,
                                                    const std::string &a_name);

    // if a_success = false, either aborts or prints a warning depending on
    // m_abort_on_catalyst_error
    void catalyst_error_or_warning(bool a_success, std::string a_msg);

    GRAMR *m_gr_amr_ptr = nullptr;
    bool m_initialised = false;

    params_t m_p;

    // variables actually passed to Catalyst on last CoProcess
    std::array<bool, NUM_VARS> m_requested_evolution_vars;
    std::array<bool, NUM_DIAGNOSTIC_VARS> m_requested_diagnostic_vars;

    vtkCPProcessor *m_proc_ptr = nullptr;
    vtkOverlappingAMR *m_vtk_grid_ptr = nullptr;
};

#endif /* USE_CATALYST */
#endif /* CATALYSTADAPTOR_HPP_ */
