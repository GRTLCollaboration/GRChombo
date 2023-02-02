/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "GRAMRLevel.hpp"

GRAMRLevel::GRAMRLevel(GRAMR &gr_amr, const SimulationParameters &a_p,
                       int a_verbosity)
    : m_gr_amr(gr_amr), m_p(a_p), m_verbosity(a_verbosity),
      m_num_ghosts(a_p.num_ghosts)
{
    if (m_verbosity)
        pout() << "GRAMRLevel constructor" << endl;
}

GRAMRLevel::~GRAMRLevel() {}

void GRAMRLevel::define(AMRLevel *a_coarser_level_ptr,
                        const Box &a_problem_domain, int a_level,
                        int a_ref_ratio)
{
    ProblemDomain physdomain(a_problem_domain);
    define(a_coarser_level_ptr, physdomain, a_level, a_ref_ratio);

    // Also define the boundaries
    m_boundaries.define(m_dx, m_p.center, m_p.boundary_params, physdomain,
                        m_num_ghosts);
}

void GRAMRLevel::define(AMRLevel *a_coarser_level_ptr,
                        const ProblemDomain &a_problem_domain, int a_level,
                        int a_ref_ratio)
{
    if (m_verbosity)
        pout() << "GRAMRLevel::define " << a_level << endl;

    AMRLevel::define(a_coarser_level_ptr, a_problem_domain, a_level,
                     a_ref_ratio);

    if (a_coarser_level_ptr)
    {
        GRAMRLevel *coarser_level_ptr = gr_cast(a_coarser_level_ptr);
        m_dx = coarser_level_ptr->m_dx / Real(coarser_level_ptr->m_ref_ratio);
    }
    else
    {
        m_dx = m_p.L / (a_problem_domain.domainBox().longside());
    }

    // Also define the boundaries
    m_boundaries.define(m_dx, m_p.center, m_p.boundary_params, a_problem_domain,
                        m_num_ghosts);
}

void GRAMRLevel::defineNewGrids(const Vector<Box> &a_new_grids)
{
    const DisjointBoxLayout level_domain = m_grids = loadBalance(a_new_grids);

    // reshape state with new grids
    IntVect iv_ghosts = m_num_ghosts * IntVect::Unit;
    m_state_new.define(level_domain, NUM_VARS, iv_ghosts);

    // maintain interlevel stuff
    defineExchangeCopier(level_domain);
    m_coarse_average.define(level_domain, NUM_VARS, m_ref_ratio);
    m_fine_interp.define(level_domain, NUM_VARS, m_ref_ratio, m_problem_domain);

    if (NUM_DIAGNOSTIC_VARS > 0)
    {
        m_state_diagnostics.define(level_domain, NUM_DIAGNOSTIC_VARS,
                                   iv_ghosts);
    }
    if (m_p.use_truncation_error_tagging)
    {
        m_fine_interp_truncation_error.define(level_domain,
                                              m_p.num_truncation_error_vars,
                                              m_ref_ratio, m_problem_domain);
        m_state_truncation_error.define(
            level_domain, 2 * m_p.num_truncation_error_vars, IntVect::Zero);
        m_state_truncation_error_old.define(
            level_domain, m_p.num_truncation_error_vars, IntVect::Zero);
        aliasLevelData(m_state_truncation_error_coarse_alias,
                       &m_state_truncation_error,
                       Interval(m_p.num_truncation_error_vars,
                                2 * m_p.num_truncation_error_vars - 1));
    }

    if (m_coarser_level_ptr != nullptr)
    {
        GRAMRLevel *coarser_gr_amr_level_ptr = gr_cast(m_coarser_level_ptr);
        m_patcher.define(level_domain, coarser_gr_amr_level_ptr->m_grids,
                         NUM_VARS, coarser_gr_amr_level_ptr->problemDomain(),
                         m_ref_ratio, m_num_ghosts);
        if (NUM_DIAGNOSTIC_VARS > 0)
        {
            m_patcher_diagnostics.define(
                level_domain, coarser_gr_amr_level_ptr->m_grids,
                NUM_DIAGNOSTIC_VARS, coarser_gr_amr_level_ptr->problemDomain(),
                m_ref_ratio, m_num_ghosts);
        }
    }
}

/// Do casting from AMRLevel to GRAMRLevel and stop if this isn't possible
const GRAMRLevel *GRAMRLevel::gr_cast(const AMRLevel *const amr_level_ptr)
{
    const GRAMRLevel *gr_amr_level_ptr =
        dynamic_cast<const GRAMRLevel *>(amr_level_ptr);
    if (gr_amr_level_ptr == nullptr)
    {
        MayDay::Error("in GRAMRLevel::gr_cast: amr_level_ptr is not castable "
                      "to GRAMRLevel*");
    }
    return gr_amr_level_ptr;
}

/// Do casting from AMRLevel to GRAMRLevel and stop if this isn't possible
GRAMRLevel *GRAMRLevel::gr_cast(AMRLevel *const amr_level_ptr)
{
    return const_cast<GRAMRLevel *>(
        gr_cast(static_cast<const AMRLevel *const>(amr_level_ptr)));
}

const GRLevelData &GRAMRLevel::getLevelData(VariableType var_type) const
{
    if (var_type == VariableType::evolution)
        return m_state_new;
    else
        return m_state_diagnostics;
}

bool GRAMRLevel::contains(const std::array<double, CH_SPACEDIM> &point) const
{
    const Box &domainBox = problemDomain().domainBox();
    for (int i = 0; i < CH_SPACEDIM; ++i)
    {
        if (point[i] < domainBox.smallEnd(i) - m_num_ghosts ||
            point[i] > domainBox.bigEnd(i) + m_num_ghosts)
        {
            return false;
        }
    }
    return true;
}

// advance by one timestep
Real GRAMRLevel::advance()
{
    CH_TIME("GRAMRLevel::advance");

    // Work out roughly how fast the evolution is going since restart
    double speed = (m_time - m_restart_time) / m_gr_amr.get_walltime();

    // Get information on number of boxes on this level (helps with better load
    // balancing)
    const DisjointBoxLayout &level_domain = m_state_new.disjointBoxLayout();
    int nbox = level_domain.dataIterator().size();
    int total_nbox = level_domain.size();
    int rank;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
    rank = 0;
#endif
    if (rank == 0)
    {
        pout() << "GRAMRLevel::advance level " << m_level << " at time "
               << m_time << " (" << speed << " M/hr)"
               << ". Boxes on this rank: " << nbox << " / " << total_nbox
               << endl;
    }

    // copy soln to old state to save it
    m_state_new.copyTo(m_state_new.interval(), m_state_old,
                       m_state_old.interval());

    // Specifically copy boundary cells
    copyBdyGhosts(m_state_new, m_state_old);

    // The level classes take flux-register parameters, use dummy ones here
    LevelFluxRegister *coarser_fr = nullptr;
    LevelFluxRegister *finer_fr = nullptr;
    // Default coarser level pointers to an empty LevelData
    // (Chombo usually checks isDefined() rather than != nullptr so we shouldn't
    // use the nullptr)
    const GRLevelData null_gr_level_data;
    const GRLevelData *coarser_data_old = &null_gr_level_data;
    const GRLevelData *coarser_data_new = &null_gr_level_data;

    Real t_coarser_old = 0.0;
    Real t_coarser_new = 0.0;

    // A coarser level exists
    if (m_coarser_level_ptr != nullptr)
    {
        GRAMRLevel *coarser_gr_amr_level_ptr = gr_cast(m_coarser_level_ptr);
        coarser_data_old = &coarser_gr_amr_level_ptr->m_state_old;
        coarser_data_new = &coarser_gr_amr_level_ptr->m_state_new;

        t_coarser_new = coarser_gr_amr_level_ptr->m_time;
        t_coarser_old = t_coarser_new - coarser_gr_amr_level_ptr->m_dt;
    }

    if (m_finer_level_ptr != nullptr)
    {
        GRAMRLevel *fine_gr_amr_level_ptr = gr_cast(m_finer_level_ptr);
        RK4LevelAdvance(m_state_new, m_state_old,
                        fine_gr_amr_level_ptr->m_patcher.getTimeInterpolator(),
                        *coarser_data_old, t_coarser_old, *coarser_data_new,
                        t_coarser_new, *coarser_fr, *finer_fr, m_time, m_dt,
                        *this);
    }
    else
    {
        RK4LevelAdvance(m_state_new, m_state_old, *coarser_data_old,
                        t_coarser_old, *coarser_data_new, t_coarser_new,
                        *coarser_fr, *finer_fr, m_time, m_dt, *this);
    }

    specificAdvance();
    // enforce solution BCs - in case of updates in specificAdvance
    fillBdyGhosts(m_state_new);

    m_time += m_dt;
    return m_dt;
}

// things to do after a timestep
void GRAMRLevel::postTimeStep()
{
    if (m_verbosity)
        pout() << "GRAMRLevel::postTimeStep " << m_level << endl;

    if (m_p.use_truncation_error_tagging && m_gr_amr.need_to_regrid(m_level))
    {
        for (int ivar = 0; ivar < m_p.num_truncation_error_vars; ++ivar)
        {
            const int var = m_p.truncation_error_vars[ivar].first;
            // copy state vars now before coarse fine averaging for
            // truncation error estimation
            m_state_new.copyTo(Interval(var, var), m_state_truncation_error_old,
                               Interval(ivar, ivar));
        }
    }

    if (m_finer_level_ptr != nullptr)
    {
        GRAMRLevel *finer_gr_amr_level_ptr = gr_cast(m_finer_level_ptr);
        finer_gr_amr_level_ptr->m_coarse_average.averageToCoarse(
            m_state_new, finer_gr_amr_level_ptr->m_state_new);
        // Synchronise times to avoid floating point errors for finer levels
        finer_gr_amr_level_ptr->time(m_time);
    }

    specificPostTimeStep();

    // enforce solution BCs - this is required after the averaging
    // and postentially after specificPostTimeStep actions
    fillBdyGhosts(m_state_new);

    if (m_verbosity)
        pout() << "GRAMRLevel::postTimeStep " << m_level << " finished" << endl;
}

// things to do before tagging cells using truncation tagging
void GRAMRLevel::preTagCellsTruncationTagging()
{
    CH_TIME("GRAMRLevel::preTagCellsTruncationTagging");
    if (m_coarser_level_ptr != nullptr)
    {
        GRAMRLevel *coarser_gr_amr_level_ptr = gr_cast(m_coarser_level_ptr);
        // interpolate from coarser level
        m_fine_interp_truncation_error.interpToFine(
            m_state_truncation_error_coarse_alias,
            coarser_gr_amr_level_ptr->m_state_truncation_error_old);
    }
    for (int ivar = 0; ivar < m_p.num_truncation_error_vars; ++ivar)
    {
        const int var = m_p.truncation_error_vars[ivar].first;
        // copy these grid variables to m_state_truncation_error so that we
        // can have all the quantities needed to compute truncation error
        // in one object (m_state_truncation_error_coarse_alias is an alias
        // to the second half)
        m_state_new.copyTo(Interval(var, var), m_state_truncation_error,
                           Interval(ivar, ivar));
    }
}

// things to do before tagging cells
void GRAMRLevel::preTagCells()
{
    CH_TIME("GRAMRLevel::preTagCells");
    if (m_p.use_truncation_error_tagging)
    {
        preTagCellsTruncationTagging();
    }
    else
        fillAllEvolutionGhosts(); // We need filled ghost cells to calculate
                                  // gradients etc
}

// tag cells that need to be refined
void GRAMRLevel::tagCells(IntVectSet &a_tags)
{
    tagCellsImplem(a_tags, m_p.use_truncation_error_tagging);
}

// create tags at initialization
void GRAMRLevel::tagCellsInit(IntVectSet &a_tags)
{
    // can't use truncation error tagging for the initial grids
    tagCellsImplem(a_tags, false);
}

// create tags
void GRAMRLevel::tagCellsImplem(IntVectSet &a_tags,
                                bool a_use_truncation_error_tagging)
{
    CH_TIME("GRAMRLevel::tagCellsImplem");
    if (m_verbosity)
        pout() << "GRAMRLevel::tagCellsImplem " << m_level << endl;

    preTagCells();

    IntVectSet local_tags;

    const DisjointBoxLayout &level_domain = m_state_new.disjointBoxLayout();
    DataIterator dit0 = level_domain.dataIterator();
    int nbox = dit0.size();
    for (int ibox = 0; ibox < nbox; ++ibox)
    {
        DataIndex di = dit0[ibox];
        const Box &b = level_domain[di];

        // mod gradient
        FArrayBox tagging_criterion(b, 1);
        if (a_use_truncation_error_tagging)
        {
            const FArrayBox &state_truncation_error_fab =
                m_state_truncation_error[di];
            computeTruncationError(tagging_criterion,
                                   state_truncation_error_fab);
        }
        else
        {
            const FArrayBox &state_fab_diagnostics = m_state_diagnostics[di];
            const FArrayBox &state_fab = m_state_new[di];
            computeTaggingCriterion(tagging_criterion, state_fab);
            computeDiagnosticsTaggingCriterion(tagging_criterion,
                                               state_fab_diagnostics);
        }

        const IntVect &smallEnd = b.smallEnd();
        const IntVect &bigEnd = b.bigEnd();

        const int xmin = smallEnd[0];
        const int ymin = smallEnd[1];
        const int zmin = smallEnd[2];

        const int xmax = bigEnd[0];
        const int ymax = bigEnd[1];
        const int zmax = bigEnd[2];

        const double threshold =
            (a_use_truncation_error_tagging)
                ? m_p.truncation_error_regrid_thresholds[m_level]
                : m_p.regrid_thresholds[m_level];

#pragma omp parallel for collapse(3) schedule(static) default(shared)
        for (int iz = zmin; iz <= zmax; ++iz)
            for (int iy = ymin; iy <= ymax; ++iy)
                for (int ix = xmin; ix <= xmax; ++ix)
                {
                    IntVect iv(ix, iy, iz);
                    if (tagging_criterion(iv, 0) >= threshold)
                    {
// local_tags |= is not thread safe.
#pragma omp critical
                        {
                            local_tags |= iv;
                        }
                    }
                }
    }

    local_tags.grow(m_p.tag_buffer_size);

    // Need to do this in two steps unless a IntVectSet::operator &=
    // (ProblemDomain) operator is defined
    Box local_tags_box = local_tags.minBox();
    local_tags_box &= m_problem_domain;
    local_tags &= local_tags_box;

    a_tags = local_tags;
}

void GRAMRLevel::computeTruncationError(
    FArrayBox &tagging_criterion, const FArrayBox &a_state_truncation_error)
{
    BoxLoops::loop(
        TruncationErrorTagging(m_p.num_truncation_error_vars, m_level),
        a_state_truncation_error, tagging_criterion);
}

// regrid
void GRAMRLevel::regrid(const Vector<Box> &a_new_grids)
{
    CH_TIME("GRAMRLevel::regrid");

    if (m_verbosity)
        pout() << "GRAMRLevel::regrid " << m_level << endl;

    m_level_grids = a_new_grids;

    mortonOrdering(m_level_grids);

    // save data for later copy, including boundary cells
    m_state_new.copyTo(m_state_new.interval(), m_state_old,
                       m_state_old.interval());

    // Specifically copy boundary cells
    copyBdyGhosts(m_state_new, m_state_old);

    // define new grids, LevelDatas, patchers and interpolators
    defineNewGrids(a_new_grids);

    if (m_coarser_level_ptr != nullptr)
    {
        GRAMRLevel *coarser_gr_amr_level_ptr = gr_cast(m_coarser_level_ptr);
        // interpolate from coarser level
        m_fine_interp.interpToFine(m_state_new,
                                   coarser_gr_amr_level_ptr->m_state_new);

        // also interpolate fine boundary cells
        if (m_p.boundary_params.nonperiodic_boundaries_exist)
        {
            m_boundaries.interp_boundaries(
                m_state_new, coarser_gr_amr_level_ptr->m_state_new, Side::Hi);
            m_boundaries.interp_boundaries(
                m_state_new, coarser_gr_amr_level_ptr->m_state_new, Side::Lo);
        }
    }

    // copy from old state
    m_state_old.copyTo(m_state_old.interval(), m_state_new,
                       m_state_new.interval());

    // Specifically copy boundary cells (only if same layout, otherwise
    // interpolated)
    copyBdyGhosts(m_state_old, m_state_new);

    // enforce solution BCs (overwriting any interpolation)
    fillBdyGhosts(m_state_new);

    IntVect iv_ghosts = m_num_ghosts * IntVect::Unit;
    m_state_old.define(m_grids, NUM_VARS, iv_ghosts);
}

/// things to do after regridding
void GRAMRLevel::postRegrid(int a_base_level)
{
    // set m_restart_time to same as the coarser level
    if (m_level > a_base_level && m_coarser_level_ptr != nullptr)
    {
        GRAMRLevel *coarser_gr_amr_level_ptr = gr_cast(m_coarser_level_ptr);
        m_restart_time = coarser_gr_amr_level_ptr->m_restart_time;
    }
}

// initialize grid
void GRAMRLevel::initialGrid(const Vector<Box> &a_new_grids)
{
    CH_TIME("GRAMRLevel::initialGrid");

    if (m_verbosity)
        pout() << "GRAMRLevel::initialGrid " << m_level << endl;

    m_level_grids = a_new_grids;

    // define new grids, LevelDatas, patchers and interpolators
    defineNewGrids(a_new_grids);
    IntVect iv_ghosts = m_num_ghosts * IntVect::Unit;
    m_state_old.define(m_grids, NUM_VARS, iv_ghosts);
}

// things to do after initialization
void GRAMRLevel::postInitialize() { m_restart_time = 0.; }

// compute dt
Real GRAMRLevel::computeDt()
{
    if (m_verbosity)
        pout() << "GRAMRLevel::computeDt " << m_level << endl;
    return m_dt;
}

// compute dt with initial data
Real GRAMRLevel::computeInitialDt()
{
    if (m_verbosity)
        pout() << "GRAMRLevel::computeInitialDt " << m_level << endl;

    m_dt = m_initial_dt_multiplier * m_dx;
    return m_dt;
}

DisjointBoxLayout GRAMRLevel::loadBalance(const Vector<Box> &a_grids)
{
    CH_TIME("GRAMRLevel::loadBalance");

    // load balance and create boxlayout
    Vector<int> procMap;

    // appears to be faster for all procs to do the loadbalance (ndk)
    LoadBalance(procMap, a_grids);

    if (m_verbosity == 1)
    {
        pout() << "GRAMRLevel::::loadBalance" << endl;
    }
    else if (m_verbosity > 1)
    {
        pout() << "GRAMRLevel::::loadBalance: procesor map: " << endl;
        for (int igrid = 0; igrid < a_grids.size(); ++igrid)
        {
            pout() << igrid << ": " << procMap[igrid] << "  " << endl;
        }
        pout() << endl;
    }

    DisjointBoxLayout dbl(a_grids, procMap, m_problem_domain);
    dbl.close();

    return dbl;
}

// write checkpoint header
#ifdef CH_USE_HDF5
void GRAMRLevel::writeCheckpointHeader(HDF5Handle &a_handle) const
{
    if (m_verbosity)
        pout() << "GRAMRLevel::writeCheckpointHeader" << endl;

    HDF5HeaderData header;
    header.m_int["num_components"] = NUM_VARS;
    char comp_str[30];
    for (int comp = 0; comp < NUM_VARS; ++comp)
    {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = UserVariables::variable_names[comp];
    }
    header.writeToFile(a_handle);

    if (m_verbosity)
        pout() << header << endl;
}

void GRAMRLevel::writeCheckpointLevel(HDF5Handle &a_handle) const
{
    CH_TIME("GRAMRLevel::writeCheckpointLevel");

    if (m_verbosity)
        pout() << "GRAMRLevel::writeCheckpointLevel" << endl;

    char level_str[20];
    sprintf(level_str, "%d", m_level);
    const std::string label = std::string("level_") + level_str;

    a_handle.setGroup(label);

    HDF5HeaderData header;

    header.m_int["ref_ratio"] = m_ref_ratio;
    header.m_int["tag_buffer_size"] = m_p.tag_buffer_size;
    header.m_real["dx"] = m_dx;
    header.m_real["dt"] = m_dt;
    header.m_real["time"] = m_time;
    header.m_box["prob_domain"] = m_problem_domain.domainBox();

    // Setup the periodicity info
    for (int dir = 0; dir < SpaceDim; ++dir)
    {
        char dir_str[20];
        sprintf(dir_str, "%d", dir);
        const std::string periodic_label =
            std::string("is_periodic_") + dir_str;
        header.m_int[periodic_label] = m_problem_domain.isPeriodic(dir);
    }

    header.writeToFile(a_handle);

    if (m_verbosity)
        pout() << header << endl;

    write(a_handle, m_state_new.boxLayout());

    // only need to write ghosts when non periodic BCs exist
    IntVect ghost_vector = IntVect::Zero;
    if (m_p.boundary_params.nonperiodic_boundaries_exist)
    {
        ghost_vector = m_num_ghosts * IntVect::Unit;
    }
    write(a_handle, m_state_new, "data", ghost_vector);
}

void GRAMRLevel::readCheckpointHeader(HDF5Handle &a_handle)
{
    CH_TIME("GRAMRLevel::readCheckpointHeader");

    if (m_verbosity)
        pout() << "GRAMRLevel::readCheckpointHeader" << endl;

    HDF5HeaderData header;
    header.readFromFile(a_handle);

    if (m_verbosity)
        pout() << "hdf5 header data:" << endl;
    if (m_verbosity)
        pout() << header << endl;

    // read number of components
    if (header.m_int.find("num_components") == header.m_int.end())
    {
        MayDay::Error("GRAMRLevel::readCheckpointHeader: checkpoint file does "
                      "not have num_components");
    }
    int num_comps = header.m_int["num_components"];
    if (num_comps != NUM_VARS)
    {
        MayDay::Error("GRAMRLevel::readCheckpointHeader: num_components in "
                      "checkpoint file does not match solver");
    }

    // read component names
    std::string state_name;
    char comp_str[60];
    for (int comp = 0; comp < NUM_VARS; ++comp)
    {
        sprintf(comp_str, "component_%d", comp);
        if (header.m_string.find(comp_str) == header.m_string.end())
        {
            MayDay::Error("GRAMRLevel::readCheckpointHeader: checkpoint file "
                          "does not have enough component names");
        }
        state_name = header.m_string[comp_str];
        if (state_name != UserVariables::variable_names[comp])
        {
            if (m_p.ignore_checkpoint_name_mismatch)
            {
                MayDay::Warning("GRAMRLevel::readCheckpointHeader: state_name "
                                "mismatch error silenced by user.");
            }
            else
                MayDay::Error("GRAMRLevel::readCheckpointHeader: state_name in "
                              "checkpoint does not match solver");
        }
    }

    // can't regrid before the first timestep on a restart if using
    // truncation error tagging so defer it.
    if (m_p.use_truncation_error_tagging)
    {
        m_gr_amr.defer_regridding();
    }
}

void GRAMRLevel::readCheckpointLevel(HDF5Handle &a_handle)
{
    CH_TIME("GRAMRLevel::readCheckpointLevel");
    if (m_verbosity)
        pout() << "GRAMRLevel::readCheckpointLevel" << endl;

    char level_str[20];
    sprintf(level_str, "%d", m_level);
    const std::string label = std::string("level_") + level_str;

    a_handle.setGroup(label);

    HDF5HeaderData header;
    header.readFromFile(a_handle);

    if (m_verbosity)
        pout() << "hdf5 header data:" << endl;
    if (m_verbosity)
        pout() << header << endl;

    // read refinement ratio
    if (header.m_int.find("ref_ratio") == header.m_int.end())
    {
        MayDay::Error(
            "GRAMRLevel::readCheckpointLevel: file does not contain ref_ratio");
    }
    m_ref_ratio = header.m_int["ref_ratio"];

    if (m_verbosity)
        pout() << "read ref_ratio = " << m_ref_ratio << endl;

    // read dx
    if (header.m_real.find("dx") == header.m_real.end())
    {
        MayDay::Error(
            "GRAMRLevel::readCheckpointLevel: file does not contain dx");
    }
    m_dx = header.m_real["dx"];

    if (m_verbosity)
        pout() << "read dx = " << m_dx << endl;

    // Since we have fixed time steping it is better to take dt from the
    // parameter file
    computeInitialDt();
    if (m_verbosity)
        pout() << "dt = " << m_dt << endl;

    // read time
    if (header.m_real.find("time") == header.m_real.end())
    {
        MayDay::Error(
            "GRAMRLevel::readCheckpointLevel: file does not contain time");
    }
    m_time = header.m_real["time"];
    m_restart_time = m_time;
    if (m_verbosity)
        pout() << "read time = " << m_time << endl;

    // read problem domain
    if (header.m_box.find("prob_domain") == header.m_box.end())
    {
        MayDay::Error("GRAMRLevel::readCheckpointLevel: file does not contain "
                      "prob_domain");
    }
    Box domainBox = header.m_box["prob_domain"];

    // Get the periodicity info
    bool isPeriodic[SpaceDim] = {
        false}; // default to false unless other info is available
    for (int dir = 0; dir < SpaceDim; ++dir)
    {
        char dir_str[20];
        sprintf(dir_str, "%d", dir);
        const std::string periodic_label =
            std::string("is_periodic_") + dir_str;
        if (!(header.m_int.find(periodic_label) == header.m_int.end()))
        {
            isPeriodic[dir] = (header.m_int[periodic_label] == true);
        }
    }
    m_problem_domain = ProblemDomain(domainBox, isPeriodic);

    // read grids
    Vector<Box> grids;
    const int grid_status = read(a_handle, grids);
    if (grid_status != 0)
    {
        MayDay::Error("GRAMRLevel::readCheckpointLevel: file does not contain "
                      "a Vector<Box>");
    }

    // define new grids, LevelDatas, patchers and interpolators
    defineNewGrids(grids);

    if (m_verbosity)
        pout() << "read level domain: " << endl;
    LayoutIterator lit = m_grids.layoutIterator();
    for (lit.begin(); lit.ok(); ++lit)
    {
        const Box &b = m_grids[lit()];
        if (m_verbosity)
            pout() << lit().intCode() << ": " << b << endl;
        m_level_grids.push_back(b);
    }
    if (m_verbosity)
        pout() << endl;

    bool redefine_data = false;
    Interval comps(0, NUM_VARS - 1);
    const int data_status = read<FArrayBox>(a_handle, m_state_new, "data",
                                            m_grids, comps, redefine_data);
    if (data_status != 0)
    {
        MayDay::Error("GRAMRLevel::readCheckpointLevel: file does not contain "
                      "state data");
    }
    IntVect iv_ghosts = m_num_ghosts * IntVect::Unit;
    m_state_old.define(m_grids, NUM_VARS, iv_ghosts);
}

void GRAMRLevel::writePlotLevel(HDF5Handle &a_handle) const
{
    if (m_verbosity)
        pout() << "GRAMRLevel::writePlotLevel" << endl;

    // number and index of states to print
    const std::vector<std::pair<int, VariableType>> &plot_states =
        m_p.plot_vars;
    int num_states = plot_states.size();

    if (num_states > 0)
    {
        // Setup the level string
        char levelStr[20];
        sprintf(levelStr, "%d", m_level);
        const std::string label = std::string("level_") + levelStr;

        a_handle.setGroup(label);

        // Setup the level header information
        HDF5HeaderData header;

        header.m_int["ref_ratio"] = m_ref_ratio;
        header.m_int["tag_buffer_size"] = m_p.tag_buffer_size;
        header.m_real["dx"] = m_dx;
        header.m_real["dt"] = m_dt;
        header.m_real["time"] = m_time;
        header.m_box["prob_domain"] = m_problem_domain.domainBox();

        // Setup the periodicity info
        for (int dir = 0; dir < SpaceDim; ++dir)
        {
            char dir_str[20];
            sprintf(dir_str, "%d", dir);
            const std::string periodic_label =
                std::string("is_periodic_") + dir_str;
            header.m_int[periodic_label] = m_problem_domain.isPeriodic(dir);
        }

        // Write the header for this level
        header.writeToFile(a_handle);

        if (m_verbosity)
            pout() << header << endl;

        const DisjointBoxLayout &levelGrids = m_state_new.getBoxes();
        IntVect iv_ghosts = m_num_ghosts * IntVect::Unit;
        LevelData<FArrayBox> plot_data(levelGrids, num_states, iv_ghosts);

        // only need to write ghosts when non periodic BCs exist
        IntVect ghost_vector = IntVect::Zero;
        if (m_p.write_plot_ghosts)
        {
            ghost_vector = m_num_ghosts * IntVect::Unit;
            Box grown_domain_box = m_problem_domain.domainBox();
            grown_domain_box.grow(ghost_vector);
            Copier boundary_copier;
            boundary_copier.ghostDefine(
                m_state_new.disjointBoxLayout(), plot_data.disjointBoxLayout(),
                grown_domain_box, ghost_vector, ghost_vector);
            for (int comp = 0; comp < num_states; comp++)
            {
                Interval currentComp(comp, comp);
                if (plot_states[comp].second == VariableType::evolution)
                {
                    Interval plotComps(plot_states[comp].first,
                                       plot_states[comp].first);
                    m_state_new.copyTo(plotComps, plot_data, currentComp,
                                       boundary_copier);
                }
                else
                {
                    Interval plotComps(plot_states[comp].first,
                                       plot_states[comp].first);
                    if (NUM_DIAGNOSTIC_VARS > 0)
                    {
                        m_state_diagnostics.copyTo(plotComps, plot_data,
                                                   currentComp);
                    }
                }
            }
        }

        else
        {
            for (int comp = 0; comp < num_states; comp++)
            {
                Interval currentComp(comp, comp);
                if (plot_states[comp].second == VariableType::evolution)
                {
                    Interval plotComps(plot_states[comp].first,
                                       plot_states[comp].first);
                    m_state_new.copyTo(plotComps, plot_data, currentComp);
                }
                else
                {
                    Interval plotComps(plot_states[comp].first,
                                       plot_states[comp].first);
                    if (NUM_DIAGNOSTIC_VARS > 0)
                    {
                        m_state_diagnostics.copyTo(plotComps, plot_data,
                                                   currentComp);
                    }
                }
            }
        }

        plot_data.exchange(plot_data.interval());

        // Write the data for this level
        write(a_handle, levelGrids);
        write(a_handle, plot_data, "data", ghost_vector);
    }
}

void GRAMRLevel::writePlotHeader(HDF5Handle &a_handle) const
{
    if (m_verbosity)
        pout() << "GRAMRLevel::writePlotHeader" << endl;

    // number and index of states to print.
    const std::vector<std::pair<int, VariableType>> &plot_states =
        m_p.plot_vars;
    int num_states = plot_states.size();

    if (num_states > 0)
    {
        // Setup the number of components
        HDF5HeaderData header;
        header.m_int["num_components"] = num_states;

        // Setup the component names
        char compStr[30];
        for (int comp = 0; comp < num_states; ++comp)
        {
            sprintf(compStr, "component_%d", comp);
            if (plot_states[comp].second == VariableType::evolution)
            {
                header.m_string[compStr] =
                    UserVariables::variable_names[plot_states[comp].first];
            }
            else
            {
                header.m_string[compStr] =
                    DiagnosticVariables::variable_names[plot_states[comp]
                                                            .first];
            }
        }

        // Write the header
        header.writeToFile(a_handle);

        if (m_verbosity)
            pout() << header << endl;
    }
    else
    {
        MayDay::Warning("GRAMRLevel::writePlotLevel: A plot interval is "
                        "provided but no components are selected for plotting. "
                        "Plot files will be empty.");
    }
}
#endif /*ifdef CH_USE_HDF5*/

void GRAMRLevel::evalRHS(GRLevelData &rhs, GRLevelData &soln,
                         LevelFluxRegister &fineFR, LevelFluxRegister &crseFR,
                         const GRLevelData &oldCrseSoln, Real oldCrseTime,
                         const GRLevelData &newCrseSoln, Real newCrseTime,
                         Real time, Real fluxWeight)
{
    CH_TIME("GRAMRLevel::evalRHS");
    if (m_verbosity)
        pout() << "GRAMRLevel::evalRHS" << endl;

    soln.exchange(m_exchange_copier);

    if (oldCrseSoln.isDefined())
    {
        // "time" falls between the old and the new coarse times
        Real alpha = (time - oldCrseTime) / (newCrseTime - oldCrseTime);

        // Assuming RK4, we know that there can only be 5 different alpha so fix
        // them with tolerance to prevent floating point problems
        Real eps = 0.01;
        if (abs(alpha) < eps)
            alpha = 0.0;
        else if (abs(alpha - 0.25) < eps)
            alpha = 0.25;
        else if (abs(alpha - 0.5) < eps)
            alpha = 0.5;
        else if (abs(alpha - 0.75) < eps)
            alpha = 0.75;
        else if (abs(alpha - 1.) < eps)
            alpha = 1.0;
        else
        {
            pout() << "alpha: " << alpha << endl;
            MayDay::Error(
                "Time interpolation coefficient is incompatible with RK4.");
        }

        // Interpolate ghost cells from next coarser level in space and time
        m_patcher.fillInterp(soln, alpha, 0, 0, NUM_VARS);
    }

    fillBdyGhosts(soln);

    specificEvalRHS(soln, rhs, time); // Call the problem specific rhs

    // evolution of the boundaries according to conditions
    if (m_p.boundary_params.nonperiodic_boundaries_exist)
    {
        m_boundaries.fill_rhs_boundaries(Side::Lo, soln, rhs);
        m_boundaries.fill_rhs_boundaries(Side::Hi, soln, rhs);
    }
}

// implements soln += dt*rhs
void GRAMRLevel::updateODE(GRLevelData &soln, const GRLevelData &rhs, Real dt)
{
    CH_TIME("GRAMRLevel::updateODE");
    // m_grown_grids will include outer boundary ghosts in the case of
    // nonperiodic BCs but will just be the problem domain otherwise.
    soln.plus(rhs, dt, m_grown_grids);

    specificUpdateODE(soln, rhs, dt);
    fillBdyGhosts(soln);
}

// define data holder newSoln based on existingSoln,
// including ghost cell specification
void GRAMRLevel::defineSolnData(GRLevelData &newSoln,
                                const GRLevelData &existingSoln)
{
    newSoln.define(existingSoln.disjointBoxLayout(), existingSoln.nComp(),
                   existingSoln.ghostVect());
}

// define data holder for RHS based on existingSoln including ghost cell
// specification (which in most cases is no ghost cells)
void GRAMRLevel::defineRHSData(GRLevelData &newRHS,
                               const GRLevelData &existingSoln)
{
    // only need ghosts for non periodic boundary case
    IntVect ghost_vector = IntVect::Zero;
    if (m_p.boundary_params.nonperiodic_boundaries_exist)
    {
        ghost_vector = m_num_ghosts * IntVect::Unit;
    }
    newRHS.define(existingSoln.disjointBoxLayout(), existingSoln.nComp(),
                  ghost_vector);
}

/// copy data in src into dest
void GRAMRLevel::copySolnData(GRLevelData &dest, const GRLevelData &src)
{
    src.copyTo(src.interval(), dest, dest.interval());

    // Specifically copy boundary cells if non periodic as
    // cells outside the domain are not copied by default
    copyBdyGhosts(src, dest);
}

double GRAMRLevel::get_dx() const { return m_dx; }

bool GRAMRLevel::at_level_timestep_multiple(int a_level) const
{
    const double coarsest_dt = m_p.coarsest_dx * m_p.dt_multiplier;
    double target_dt = coarsest_dt;
    for (int ilevel = 0; ilevel < a_level; ++ilevel)
    {
        target_dt /= m_p.ref_ratios[ilevel];
    }
    // get difference to nearest multiple of target_dt
    const double time_remainder = remainder(m_time, target_dt);
    return (abs(time_remainder) < m_gr_amr.timeEps() * coarsest_dt);
}

void GRAMRLevel::fillAllGhosts(const VariableType var_type,
                               const Interval &a_comps)
{
    if (var_type == VariableType::evolution)
    {
        Interval comps(a_comps.begin(),
                       std::min<int>(NUM_VARS - 1, a_comps.end()));
        fillAllEvolutionGhosts(comps);
    }
    else if (var_type == VariableType::diagnostic)
    {
        Interval comps(a_comps.begin(),
                       std::min<int>(NUM_DIAGNOSTIC_VARS - 1, a_comps.end()));
        fillAllDiagnosticsGhosts(comps);
    }
}

void GRAMRLevel::fillAllEvolutionGhosts(const Interval &a_comps)
{
    CH_TIME("GRAMRLevel::fillAllEvolutionGhosts()");
    if (m_verbosity)
        pout() << "GRAMRLevel::fillAllEvolutionGhosts" << endl;

    // If there is a coarser level then interpolate undefined ghost cells
    if (m_coarser_level_ptr != nullptr)
    {
        GRAMRLevel *coarser_gr_amr_level_ptr = gr_cast(m_coarser_level_ptr);
        m_patcher.fillInterp(m_state_new, coarser_gr_amr_level_ptr->m_state_new,
                             a_comps.begin(), a_comps.begin(), a_comps.size());
    }
    fillIntralevelGhosts(a_comps);
}

void GRAMRLevel::fillAllDiagnosticsGhosts(const Interval &a_comps)
{
    CH_TIME("GRAMRLevel::fillAllDiagnosticsGhosts");
    if (m_verbosity)
        pout() << "GRAMRLevel::fillAllDiagnosticsGhosts" << endl;

    // If there is a coarser level then interpolate undefined ghost cells
    if (m_coarser_level_ptr != nullptr)
    {
        GRAMRLevel *coarser_gr_amr_level_ptr = gr_cast(m_coarser_level_ptr);
        m_patcher_diagnostics.fillInterp(
            m_state_diagnostics, coarser_gr_amr_level_ptr->m_state_diagnostics,
            a_comps.begin(), a_comps.begin(), a_comps.size());
    }
    m_state_diagnostics.exchange(a_comps, m_exchange_copier);

    // We should always fill the boundary ghosts to avoid nans
    // if we have non periodic directions
    if (m_p.boundary_params.nonperiodic_boundaries_exist)
    {
        m_boundaries.fill_diagnostic_boundaries(Side::Hi, m_state_diagnostics,
                                                a_comps);
        m_boundaries.fill_diagnostic_boundaries(Side::Lo, m_state_diagnostics,
                                                a_comps);
    }
}

void GRAMRLevel::fillIntralevelGhosts(const Interval &a_comps)
{
    m_state_new.exchange(a_comps, m_exchange_copier);
    fillBdyGhosts(m_state_new);
}

void GRAMRLevel::fillBdyGhosts(GRLevelData &a_state, const Interval &a_comps)
{
    // enforce solution BCs after filling ghosts
    if (m_p.boundary_params.boundary_solution_enforced)
    {
        m_boundaries.fill_solution_boundaries(Side::Hi, a_state, a_comps);
        m_boundaries.fill_solution_boundaries(Side::Lo, a_state, a_comps);
    }
}

void GRAMRLevel::copyBdyGhosts(const GRLevelData &a_src, GRLevelData &a_dest)
{
    // Specifically copy boundary cells if non periodic as
    // cells outside the domain are not copied by default
    if (m_p.boundary_params.nonperiodic_boundaries_exist)
    {
        m_boundaries.copy_boundary_cells(Side::Hi, a_src, a_dest);
        m_boundaries.copy_boundary_cells(Side::Lo, a_src, a_dest);
    }
}

void GRAMRLevel::defineExchangeCopier(const DisjointBoxLayout &a_level_grids)
{
    // if there are Sommerfeld BCs, expand boxes along those sides
    if (m_p.boundary_params.nonperiodic_boundaries_exist)
    {
        m_boundaries.expand_grids_to_boundaries(m_grown_grids, a_level_grids);
    }
    else
    { // nothing to do if periodic BCs
        m_grown_grids = a_level_grids;
    }

    IntVect iv_ghosts = m_num_ghosts * IntVect::Unit;
    m_exchange_copier.exchangeDefine(m_grown_grids, iv_ghosts);
}
