/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "CH_Timer.H"
#include "parstream.H" //Gives us pout()

// System includes
#include <chrono>
#include <iostream>

// Our includes
#include "DefaultLevelFactory.hpp"
#include "GRParmParse.hpp"
#include "MultiLevelTask.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"
// TPAMR.hpp includes BHAMR.hpp
#include "TPAMR.hpp" // TPAMR code conditional compiled on USE_TWOPUNCTURES

// Problem specific includes:
#include "BinaryBHLevel.hpp"

// Chombo namespace
#include "UsingNamespace.H"

int runGRChombo(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);

    if (sim_params.just_check_params)
        return 0;

#ifdef USE_TWOPUNCTURES
    TPAMR bh_amr;
    bh_amr.set_two_punctures_parameters(sim_params.tp_params);
    // Run TwoPunctures solver
    bh_amr.m_two_punctures.Run();
#else
    BHAMR bh_amr;
#endif

    // must be before 'setupAMRObject' to define punctures for tagging criteria
    if (sim_params.track_punctures)
    {
        // the tagging criterion used in this example means that the punctures
        // should be on the max level but let's fill ghosts on the level below
        // too just in case
        int puncture_tracker_min_level = sim_params.max_level - 1;
        bh_amr.m_puncture_tracker.initial_setup(
            {sim_params.bh1_params.center, sim_params.bh2_params.center},
            "punctures", sim_params.data_path, puncture_tracker_min_level);
    }

    // The line below selects the problem that is simulated
    // (To simulate a different problem, define a new child of AMRLevel
    // and an associated LevelFactory)
    DefaultLevelFactory<BinaryBHLevel> binary_bh_level_fact(bh_amr, sim_params);
    setupAMRObject(bh_amr, binary_bh_level_fact);

    // call this after amr object setup so grids known
    // and need it to stay in scope throughout run
    AMRInterpolator<Lagrange<4>> interpolator(
        bh_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    bh_amr.set_interpolator(
        &interpolator); // also sets puncture_tracker interpolator

    // must be after interpolator is set
    if (sim_params.track_punctures)
        bh_amr.m_puncture_tracker.restart_punctures();

    using Clock = std::chrono::steady_clock;
    using Minutes = std::chrono::duration<double, std::ratio<60, 1>>;

    std::chrono::time_point<Clock> start_time = Clock::now();

    // Add a scheduler to call specificPostTimeStep on every AMRLevel at t=0
    auto task = [](GRAMRLevel *level)
    {
        if (level->time() == 0.)
            level->specificPostTimeStep();
    };
    // call 'now' really now
    MultiLevelTaskPtr<> call_task(task);
    call_task.execute(bh_amr);
    // or call at post-plotLevel, at every 'some_interval'
    // int some_interval = 10;
    // bool reverse_levels = true;
    // MultiLevelTaskPtr<> call_task(task, reverse_levels, some_interval);
    // bh_amr.schedule(call_task);

    bh_amr.run(sim_params.stop_time, sim_params.max_steps);

    auto now = Clock::now();
    auto duration = std::chrono::duration_cast<Minutes>(now - start_time);
    pout() << "Total simulation time (mins): " << duration.count() << ".\n";

    bh_amr.conclude();

    CH_TIMER_REPORT(); // Report results when running with Chombo timers.

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRChombo(argc, argv);

    if (status == 0)
        pout() << "GRChombo finished." << std::endl;
    else
        pout() << "GRChombo failed with return code " << status << std::endl;

    mainFinalize();
    return status;
}
