/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef USE_AHFINDER

#include "PETScCommunicator.hpp"

void PETScCommunicator::set_num_ranks(int a_num_ranks)
{
    if (m_initialized)
        return;

#ifdef CH_MPI
    if (m_mpi_group != MPI_GROUP_NULL)
        MPI_Group_free(&m_mpi_group);

    if (m_mpi_comm != MPI_COMM_NULL)
        MPI_Comm_free(&m_mpi_comm);

    if (a_num_ranks > 0) // otherwise use the whole Chombo communicator
    {
        // don't make more ranks than there are processes
        int size = numProc(); // declared in Chombo's SPMD.H
        a_num_ranks = std::min(a_num_ranks, size);

        int rank;
        if (procID() == 0) // declared in Chombo's SPMD.H
            std::cout << "Using PETSc with " << a_num_ranks << " ranks"
                      << std::endl;

        MPI_Group MPI_GROUP_WORLD;
        MPI_Comm_group(Chombo_MPI::comm, &MPI_GROUP_WORLD);

        // (TF): for things like 'SmallDataIO', rank 0 of Chombo_MPI::comm
        // is the one that prints (not rank 0 from PETSC_COMM_WORLD) so pay
        // attention before taking rank 0 of Chombo out of PETSc. As of
        // 26/05/2019, I tested with 'int range[3] = { 1, a_num_ranks , 1 }'
        // and all the AH stuff was compatible with doing so
        // Not true anymore (19/08/2020) due to how
        // ApparentHorizon::write_coords_file uses rank 0 of PETSc to transfer
        // data and to write, but SmallDataIO uses rank 0 of Chombo to write.
        int range[3] = {0, a_num_ranks - 1, 1};
        MPI_Group_range_incl(MPI_GROUP_WORLD, 1,
                             reinterpret_cast<int(*)[3]>(&range), &m_mpi_group);
        MPI_Group_free(&MPI_GROUP_WORLD);

        MPI_Comm_create(Chombo_MPI::comm, m_mpi_group, &m_mpi_comm);
    }
#endif
}

//! initialize PETSc and its MPI sub-communicator
PetscErrorCode PETScCommunicator::initialize(int a_num_ranks)
{
    set_num_ranks(a_num_ranks);

    PetscErrorCode err = 0;
    if (!m_initialized)
    {
#ifdef CH_MPI
        if (m_mpi_comm != MPI_COMM_NULL)
            PETSC_COMM_WORLD = m_mpi_comm;
        else // use Chombo communicator if no 'set_num_ranks' was called (or
             // if called with 'a_num_ranks'<=0)
            PETSC_COMM_WORLD = Chombo_MPI::comm;
#endif

        if (PETScCommunicator::is_rank_active())
            err = PetscInitializeNoArguments();

        if (!err)
            m_initialized = true;
    }
    return err;
}

//! finalize PETSc
PetscErrorCode PETScCommunicator::finalize()
{
    PetscErrorCode err = 0;
    if (m_initialized)
    {
        if (PETScCommunicator::is_rank_active())
            err = PetscFinalize();

        if (!err)
            m_initialized = false;
    }
    return err;
}

//! true if part of PETSc MPI sub-communicator
bool PETScCommunicator::is_rank_active()
{
#ifdef CH_MPI
    if (m_mpi_group == MPI_GROUP_NULL)
        return true;
    else
    {
        int rank;
        MPI_Group_rank(m_mpi_group, &rank);
        return rank != MPI_UNDEFINED;
    }
#else
    return true;
#endif
}

/////////////////////////////////////////////////////////
// initialize some "static" variables of PETScCommunicator class
/////////////////////////////////////////////////////////

bool PETScCommunicator::m_initialized = false;

#ifdef CH_MPI
MPI_Group PETScCommunicator::m_mpi_group = MPI_GROUP_NULL;
MPI_Comm PETScCommunicator::m_mpi_comm = MPI_COMM_NULL;
#endif

#endif
