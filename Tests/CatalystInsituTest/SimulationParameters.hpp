/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"

class SimulationParameters : public ChomboParameters
{
  public:
    SimulationParameters(GRParmParse &pp) : ChomboParameters(pp)
    {
        pp.load("error_threshold", error_threshold, 0.15);
        pp.load("valid_png_file", valid_png_file);
        pp.load("generated_png_file", generated_png_file);
        if (!output_path.empty())
        {
            generated_png_file = output_path + generated_png_file;
        }
    }

    double error_threshold;
    std::string valid_png_file;
    std::string generated_png_file;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
