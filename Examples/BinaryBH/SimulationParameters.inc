#ifndef SIMULATIONPARAMETERS_HPP_
#error  "This file must only be included through SimulationParameters.hpp"
#endif

void auto_read_params(GRParmParse& pp)
{
    pp.load("verbosity", verbosity);
    //Grid setup
    pp.load("L", L);
    pp.load("regrid_threshold", regrid_threshold);
    pp.load("isPeriodic", isPeriodic);
    pp.load("num_ghosts", num_ghosts);
    pp.load("tag_buffer_size", tag_buffer_size);
    pp.load("dt_multiplier", dt_multiplier);

    // Lapse evolution
    pp.load("lapse_advec_coeff",lapse_advec_coeff);
    pp.load("lapse_coeff", lapse_coeff);
    pp.load("lapse_power", lapse_power);

    // Shift Evolution
    pp.load("shift_advec_coeff",shift_advec_coeff);
    pp.load("shift_Gamma_coeff",shift_Gamma_coeff);
    pp.load("eta",eta);
    pp.load("SpatialBetaDriverRadius",SpatialBetaDriverRadius);

    //CCZ4 parameters
    pp.load("kappa1",kappa1);
    pp.load("kappa2",kappa2);
    pp.load("kappa3",kappa3);
    pp.load("covariantZ4",covariantZ4);

    //Dissipation
    pp.load("sigma",sigma);

    //Initial data
    pp.load("massA", massA);
    pp.load("centerA", centerA);
    pp.load("momentumA", momentumA);
    pp.load("massB", massB);
    pp.load("centerB", centerB);
    pp.load("momentumB", momentumB);

    //Fill in BinaryBHParameters
    bh1_params.mass = massA;
    bh1_params.center = centerA;
    bh1_params.momentum = momentumA;
    bh2_params.mass = massB;
    bh2_params.center = centerB;
    bh2_params.momentum = momentumB;

    //Misc
    pp.load("nan_check", nan_check);
    pp.load("ignore_checkpoint_name_mismatch", ignore_checkpoint_name_mismatch, false);
}

////// General parameters
int verbosity;
Real L; //Physical sidelength of the grid
int num_ghosts;
int tag_buffer_size; //Amount the tagged region is grown by
std::array<bool, CH_SPACEDIM> isPeriodic;
int nan_check;
bool ignore_checkpoint_name_mismatch; //ignore mismatch of variable names between restart file and program
Real dt_multiplier;

///// Problem specific parameters
Real regrid_threshold;
//Lapse evolution
Real lapse_advec_coeff, lapse_power, lapse_coeff;
//ShiftEvolution
Real shift_advec_coeff, shift_Gamma_coeff, eta, SpatialBetaDriverRadius;
//CCZ4 parameters
Real kappa1, kappa2, kappa3;
int covariantZ4;
Real sigma; //Dissipation
//Initial data
Real massA, massB;
std::array<double, CH_SPACEDIM> centerA, centerB;
std::array<double, CH_SPACEDIM> momentumA, momentumB;
