#ifndef SIMULATIONPARAMETERS_HPP_
#error  "This file must only be included through SimulationParameters.hpp"
#endif

void auto_read_params(GRParmParse& pp)
{
    pp.load("verbosity", verbosity);
    //Grid setup
    pp.load("L", L);
    pp.load("regrid_threshold", regrid_threshold);
    pp.load("isPeriodic", isPeriodic, {true, true, true});
    pp.load("num_ghosts", num_ghosts, 3);
    pp.load("tag_buffer_size", tag_buffer_size, 3);
    pp.load("dt_multiplier", dt_multiplier);

    // Lapse evolution
    pp.load("lapse_advec_coeff", lapse_advec_coeff);
    pp.load("lapse_coeff", lapse_coeff);
    pp.load("lapse_power", lapse_power);

    // Shift Evolution
    pp.load("shift_advec_coeff", shift_advec_coeff);
    pp.load("shift_Gamma_coeff",shift_Gamma_coeff);
    pp.load("eta",eta);

    //CCZ4 parameters
    pp.load("formulation", formulation);
    pp.load("kappa1",kappa1);
    pp.load("kappa2",kappa2);
    pp.load("kappa3",kappa3);
    pp.load("covariantZ4",covariantZ4);

    //Dissipation
    pp.load("sigma",sigma);

    //Initial Kerr data
    pp.load("kerr_mass", kerr_mass);
    pp.load("kerr_spin", kerr_spin);
    pp.load("kerr_center", kerr_center, {0.5*L, 0.5*L, 0.5*L});

    //Misc
    pp.load("nan_check", nan_check, 1);
    pp.load("ignore_checkpoint_name_mismatch", ignore_checkpoint_name_mismatch, false);
}

// General parameters
int verbosity;
Real L; //Physical sidelength of the grid
int num_ghosts; //Need min 3 for KO dissipation
int tag_buffer_size; //Amount the tagged region is grown by
std::array<bool, CH_SPACEDIM> isPeriodic;
int nan_check;
bool ignore_checkpoint_name_mismatch; //ignore mismatch of variable names between restart file and program
Real dt_multiplier;

// Problem specific parameters
Real regrid_threshold;
//Lapse evolution
Real lapse_power, lapse_coeff, lapse_advec_coeff;
//ShiftEvolution
Real shift_advec_coeff, shift_Gamma_coeff, eta;
//CCZ4 parameters
int formulation;
double kappa1, kappa2, kappa3;
int covariantZ4;
Real sigma; //Dissipation
//Initial data for Kerr
Real kerr_mass, kerr_spin;
std::array<double, CH_SPACEDIM> kerr_center;
