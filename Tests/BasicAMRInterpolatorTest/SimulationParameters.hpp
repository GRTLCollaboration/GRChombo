/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ParmParse.H"

class SimulationParameters
{
  public:
    // For the Interpolator test we don't need any parameters
    SimulationParameters(ParmParse &pp) { readParams(pp); }

    void readParams(ParmParse &pp)
    {
        pp.get("verbosity", verbosity);
        // Grid setup
        pp.get("L", L);
        pp.getarr("isPeriodic", isPeriodic, 0, SpaceDim);
        pp.get("num_ghosts", num_ghosts);
    }
    int verbosity;
    Real L; // Physical sidelength of the grid
    int num_ghosts;
    std::vector<bool> isPeriodic;
    double regrid_threshold = 0;
    int tag_buffer_size = 0;
    bool ignore_checkpoint_name_mismatch = false;
    double dt_multiplier = 0.2; // Doesn't matter for this test
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
