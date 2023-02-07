/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to LICENSE, in Chombo's root directory.
 */
#endif

// Chombo includes
#include "parstream.H" //Gives us pout()

// General includes:
#include <algorithm>
#include <cmath>
#include <cstdio>

#include "GRAMR.hpp"

#include "GRParmParse.hpp"
#include "MultiLevelTask.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "DefaultLevelFactory.hpp"
#include "InsituTestLevel.hpp"
#include "UserVariables.hpp"

// VTK/Paraview includes
#ifdef USE_CATALYST
#include <vtkImageDifference.h>
#include <vtkNew.h>
#include <vtkPNGReader.h>
#endif

// Chombo namespace
#include "UsingNamespace.H"

#ifdef USE_CATALYST

// A simple function to compare two timespecs
bool operator<=(const timespec &lhs, const timespec &rhs)
{
    if (lhs.tv_sec == rhs.tv_sec)
        return lhs.tv_nsec <= rhs.tv_nsec;
    else
        return lhs.tv_sec <= rhs.tv_sec;
}

int runInsituTest(int argc, char *argv[])
{
    int status = 0;

    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    std::string in_string = argv[argc - 1];
    pout() << in_string << std::endl;
    char const *in_file = argv[argc - 1];
    GRParmParse pp(0, argv + argc, NULL, in_file);
    SimulationParameters sim_params(pp);

    // Delete the generated PNG file if it already exists
    std::remove(sim_params.generated_png_file.c_str());

    GRAMR gr_amr;
    DefaultLevelFactory<InsituTestLevel> insitu_test_level_factory(gr_amr,
                                                                   sim_params);
    setupAMRObject(gr_amr, insitu_test_level_factory);

    auto task = [](GRAMRLevel *a_level)
    {
        if (a_level->time() == 0.0)
        {
            a_level->catalystCoProcess();
        }
    };

    MultiLevelTaskPtr<> call_task(task);
    call_task.execute(gr_amr);
    gr_amr.run(sim_params.stop_time, sim_params.max_steps);
    gr_amr.conclude();

    // read the newly generated PNG and the expected PNG
    vtkNew<vtkPNGReader> generated_png_reader;
    generated_png_reader->SetFileName(sim_params.generated_png_file.c_str());
    vtkNew<vtkPNGReader> valid_png_reader;
    valid_png_reader->SetFileName(sim_params.valid_png_file.c_str());

    // Calculate the difference
    vtkNew<vtkImageDifference> image_difference;
    image_difference->SetInputConnection(generated_png_reader->GetOutputPort());
    image_difference->SetImageConnection(valid_png_reader->GetOutputPort());
    image_difference->Update();

    double error = image_difference->GetThresholdedError();

    pout() << "error           = " << error << "\n";
    pout() << "error_threshold = " << sim_params.error_threshold << std::endl;

    status = (error < sim_params.error_threshold) ? 0 : 1;

    return status;
}
#endif

int main(int argc, char *argv[])
{
#ifdef USE_CATALYST
    mainSetup(argc, argv);

    int status = runInsituTest(argc, argv);
#else
    int status = -1;
#endif

    if (status == 0)
        pout() << "Catalyst Insitu test passed." << std::endl;
    else if (status == -1)
    {
        pout() << "Catalyst Insitu test skipped (USE_CATALYST undefined)."
               << std::endl;
    }
    else if (status == -2)
    {
        pout() << "Catalyst Insitu test skipped (ParaView version < 5.9)."
               << std::endl;
    }
    else
        pout() << "Catalyst Insitu test failed with return code " << status
               << std::endl;
#ifdef USE_CATALYST
    mainFinalize();
#endif

    return std::max(status, 0);
}
