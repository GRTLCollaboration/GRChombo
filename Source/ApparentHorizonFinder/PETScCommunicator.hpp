/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _PETSCCOMMUNICATOR_HPP_
#define _PETSCCOMMUNICATOR_HPP_

// Chombo includes
#include "SPMD.H" //Chombo_MPI::comm

#ifdef CH_MPI
#include <mpi.h>
#endif
#include <petsc.h>

// Chombo namespace
#include "UsingNamespace.H"

//! Class to control PETSc MPI sub-communicator
class PETScCommunicator
{
  public:
    //! true if part of PETSc MPI sub-communicator
    static bool is_rank_active();

    //! initialize PETSc and its MPI sub-communicator
    static PetscErrorCode initialize(int a_num_ranks);

    //! finalize PETSc
    static PetscErrorCode finalize();

  private:
    static bool m_initialized; //!< is initialized?

#ifdef CH_MPI
    static MPI_Group m_mpi_group; //!< set of MPI ranks for PETSc
    static MPI_Comm m_mpi_comm;   //!< MPI sub-communicator
#endif

    //! define number of ranks of PETSc sub-communicator
    static void set_num_ranks(int a_num_ranks);
}; // namespace AHFinder

#endif /* _PETSCCOMMUNICATOR_HPP_ */
