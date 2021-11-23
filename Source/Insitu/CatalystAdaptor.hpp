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
#include "MayDay.H"

// GRChombo includes
#include "GRAMR.hpp"
#include "UserVariables.hpp"

// ParaView/VTK/Catalyst includes
#include <vtkAMRBox.h>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkNew.h>
#include <vtkOverlappingAMR.h>
#include <vtkUniformGrid.h>

// Chombo namespace
#include "UsingNamespace.H"

// Forward declaration of GRAMR class
class GRAMR;

// Class that interfaces with ParaView Catalyst for insitu visualisation
class CatalystAdaptor
{
  public:
    // empty constructor (doesn't call initialise)
    CatalystAdaptor();

    // full constructor (calls initialise)
    CatalystAdaptor(const GRAMR *a_gr_amr_ptr,
                    std::string a_python_script_path);

    // destructor
    ~CatalystAdaptor();

    // Initialisation/Finalisation
    void initialise(const GRAMR *m_gr_amr_ptr,
                    std::string a_python_script_path);
    void finalise();

    // update the AMR grid (no grid data)
    void update_vtk_grid();

    // send variables to catalyst
    void add_var();

    // do Catalyst processing
    void coprocess(double time);

  private:
    bool m_initialised = false;
    const GRAMR *m_gr_amr_ptr = nullptr;
    std::string m_python_script_path;

    vtkCPProcessor *m_proc_ptr = nullptr;
    vtkOverlappingAMR *m_vtk_grid_ptr = nullptr;
};

#endif /* USE_CATALYST */
#endif /* CATALYSTADAPTOR_HPP_ */
