/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(PUNCTURETRACKER_HPP_)
#error "This file should only be included through PunctureTracker.hpp"
#endif

#ifndef PUNCTURETRACKER_IMPL_HPP_
#define PUNCTURETRACKER_IMPL_HPP_

//! Set punctures post restart
void PunctureTracker::restart_punctures(
    BHAMR &a_gramr,
    std::vector<std::array<double, CH_SPACEDIM>> initial_puncture_coords) const
{
    if (m_time == 0)
    {
        // if it is the first timestep, use the param values
        // rather than look for the output file, e.g. for when
        // restart from IC solver checkpoint
        set_initial_punctures(a_gramr, initial_puncture_coords);
    }
    else
    {
        // look for the current puncture location in the
        // puncture output file (it needs to exist!)
        read_in_punctures(a_gramr);
    }
}

//! set and write initial puncture locations
void PunctureTracker::set_initial_punctures(
    BHAMR &a_gramr,
    std::vector<std::array<double, CH_SPACEDIM>> initial_puncture_coords) const
{
    // first set the puncture data, initial shift is always zero
    std::vector<std::array<double, CH_SPACEDIM>> initial_shift;
    initial_shift.resize(m_num_punctures);
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        FOR1(i) { initial_shift[ipuncture][i] = 0.0; }
    }
    a_gramr.set_puncture_data(initial_puncture_coords, initial_shift);

    // now the write out to a new file
    bool first_step = true;
    SmallDataIO punctures_file(m_punctures_filename, m_dt, m_time,
                               m_restart_time, SmallDataIO::APPEND, first_step);
    std::vector<std::string> header1_strings(CH_SPACEDIM * m_num_punctures);
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        header1_strings[CH_SPACEDIM * ipuncture + 0] = "x";
        header1_strings[CH_SPACEDIM * ipuncture + 1] = "y";
        header1_strings[CH_SPACEDIM * ipuncture + 2] = "z";
    }
    punctures_file.write_header_line(header1_strings);

    // use a vector for the write out
    std::vector<double> puncture_vector =
        get_puncture_vector(initial_puncture_coords);
    punctures_file.write_time_data_line(puncture_vector);
}

//! Set punctures post restart
void PunctureTracker::read_in_punctures(BHAMR &a_gramr) const
{
    std::vector<std::array<double, CH_SPACEDIM>> puncture_coords;
    puncture_coords.resize(m_num_punctures);

    // read them in from the Punctures file at current time m_time
    // NB opening in APPEND mode allows reading where m_restart_time
    // is greater than zero and m_time < m_restart_time + m_dt
    bool first_step = false;
    SmallDataIO punctures_file(m_punctures_filename, m_dt, m_time,
                               m_restart_time, SmallDataIO::APPEND, first_step);
    // NB need to give the get function an empty vector to fill
    std::vector<double> puncture_vector;
    punctures_file.get_specific_data_line(puncture_vector, m_time);
    // check the data returned is the right size
    CH_assert(puncture_vector.size() == m_num_punctures * CH_SPACEDIM);
    // remove any duplicate data from the file
    const bool keep_m_time_data = true;
    punctures_file.remove_duplicate_time_data(keep_m_time_data);

    // convert vector to list of coords
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        puncture_coords[ipuncture] = {
            puncture_vector[ipuncture * CH_SPACEDIM + 0],
            puncture_vector[ipuncture * CH_SPACEDIM + 1],
            puncture_vector[ipuncture * CH_SPACEDIM + 2]};
    }

    // set the coordinates and get the current shift
    std::vector<std::array<double, CH_SPACEDIM>> current_shift;
    get_interp_shift(current_shift, a_gramr, puncture_coords);
    a_gramr.set_puncture_data(puncture_coords, current_shift);

    // print out values into pout files
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        pout() << "Puncture " << ipuncture
               << " restarted at : " << puncture_coords[ipuncture][0] << " "
               << puncture_coords[ipuncture][1] << " "
               << puncture_coords[ipuncture][2] << endl;
        pout() << " with shift vector : " << current_shift[ipuncture][0] << " "
               << current_shift[ipuncture][1] << " "
               << current_shift[ipuncture][2] << endl;
        pout() << "at time = " << m_time << endl;
    }
}

//! Execute the tracking and write out
void PunctureTracker::execute_tracking(BHAMR &a_gramr,
                                       const bool write_punctures) const
{
    CH_TIME("PunctureTracker::execute_tracking");
    // get puncture coordinates and old shift value
    std::vector<std::array<double, CH_SPACEDIM>> puncture_coords =
        a_gramr.get_puncture_coords();
    std::vector<std::array<double, CH_SPACEDIM>> old_shift =
        a_gramr.get_puncture_shift();
    CH_assert(puncture_coords.size() == m_num_punctures); // sanity check

    // new shift value
    std::vector<std::array<double, CH_SPACEDIM>> new_shift;
    get_interp_shift(new_shift, a_gramr, puncture_coords);

    // update puncture locations using second order update
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        FOR1(i)
        {
            puncture_coords[ipuncture][i] +=
                -0.5 * m_dt *
                (new_shift[ipuncture][i] + old_shift[ipuncture][i]);
        }
    }
    a_gramr.set_puncture_data(puncture_coords, new_shift);

    // print them out
    if (write_punctures)
    {
        bool first_step = false;
        SmallDataIO punctures_file(m_punctures_filename, m_dt, m_time,
                                   m_restart_time, SmallDataIO::APPEND,
                                   first_step);

        // use a vector for the write out
        std::vector<double> puncture_vector =
            get_puncture_vector(puncture_coords);
        punctures_file.write_time_data_line(puncture_vector);
    }
}

//! Use the interpolator to get the value of the shift at
//! given coords
void PunctureTracker::get_interp_shift(
    std::vector<std::array<double, CH_SPACEDIM>> &interp_shift, BHAMR &a_gramr,
    std::vector<std::array<double, CH_SPACEDIM>> puncture_coords) const
{
    CH_TIME("PunctureTracker::get_interp_shift");
    // resize the vector to the number of punctures
    interp_shift.resize(m_num_punctures);

    // refresh interpolator
    a_gramr.m_interpolator->refresh();

    // set up shift and coordinate holders
    std::vector<double> interp_shift1(m_num_punctures);
    std::vector<double> interp_shift2(m_num_punctures);
    std::vector<double> interp_shift3(m_num_punctures);
    std::vector<double> interp_x(m_num_punctures);
    std::vector<double> interp_y(m_num_punctures);
    std::vector<double> interp_z(m_num_punctures);

    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        interp_x[ipuncture] = puncture_coords[ipuncture][0];
        interp_y[ipuncture] = puncture_coords[ipuncture][1];
        interp_z[ipuncture] = puncture_coords[ipuncture][2];
    }

    // setup query
    InterpolationQuery query(m_num_punctures);
    query.setCoords(0, interp_x.data())
        .setCoords(1, interp_y.data())
        .setCoords(2, interp_z.data())
        .addComp(c_shift1, interp_shift1.data())
        .addComp(c_shift2, interp_shift2.data())
        .addComp(c_shift3, interp_shift3.data());

    // engage!
    a_gramr.m_interpolator->interp(query);

    // put the shift values into the output array
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        interp_shift[ipuncture] = {interp_shift1[ipuncture],
                                   interp_shift2[ipuncture],
                                   interp_shift3[ipuncture]};
    }
}

//! get a vector of the puncture coords - used for write out
std::vector<double> PunctureTracker::get_puncture_vector(
    std::vector<std::array<double, CH_SPACEDIM>> puncture_coords) const
{
    std::vector<double> puncture_vector;
    puncture_vector.resize(m_num_punctures * CH_SPACEDIM);
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        puncture_vector[ipuncture * CH_SPACEDIM + 0] =
            puncture_coords[ipuncture][0];
        puncture_vector[ipuncture * CH_SPACEDIM + 1] =
            puncture_coords[ipuncture][1];
        puncture_vector[ipuncture * CH_SPACEDIM + 2] =
            puncture_coords[ipuncture][2];
    }
    return puncture_vector;
}

#endif /* PUNCTURETRACKER_IMPL_HPP_ */
