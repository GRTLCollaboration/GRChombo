/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MPICONTEXT_IMPL_HPP_
#define MPICONTEXT_IMPL_HPP_

inline MPIContext::MPIContext()
    : m_num_process(comm_size()), m_rank(comm_rank()), m_query(m_num_process),
      m_answer(m_num_process), m_async_active(false)
{
}

inline int MPIContext::queryCount(int rank) { return m_query.count(rank); }

inline int MPIContext::totalQueryCount() { return m_query.totalCount(); }

inline int MPIContext::answerCount(int rank) { return m_answer.count(rank); }

inline int MPIContext::totalAnswerCount() { return m_answer.totalCount(); }

inline int MPIContext::queryDispl(int rank) { return m_query.displ(rank); }

inline int MPIContext::answerDispl(int rank) { return m_answer.displ(rank); }

inline void MPIContext::setQueryCount(int rank, int count)
{
    CH_assert(!m_async_active);
    m_query.setCount(rank, count);
}

inline void MPIContext::incrementQueryCount(int rank)
{
    CH_assert(!m_async_active);
    m_query.incrementCount(rank);
}

inline void MPIContext::clearQueryCounts()
{
    CH_assert(!m_async_active);
    m_query.clearCounts();
}

inline void MPIContext::exchangeLayout()
{
    CH_assert(!m_async_active);
#ifdef CH_MPI
    MPI_Alltoall(m_query.countsPtr(), 1, MPI_INT, m_answer.countsPtr(), 1,
                 MPI_INT, Chombo_MPI::comm);
#else
    *m_answer.countsPtr() = *m_query.countsPtr();
#endif
    m_answer.updateDirty();
}

#ifdef CH_MPI
inline void MPIContext::asyncBegin()
{
    CH_assert(!m_async_active);
    m_async_active = true;
}

inline void MPIContext::asyncExchangeQuery(void *sendbuf, void *recvbuf,
                                           MPI_Datatype type)
{
    CH_assert(m_async_active);
    MPI_Request req;
    m_mpi_requests.push_back(req);

#if MPI_VERSION >= 3 && !defined(OPEN_MPI)
    MPI_Ialltoallv(sendbuf, m_query.countsPtr(), m_query.displsPtr(), type,
                   recvbuf, m_answer.countsPtr(), m_answer.displsPtr(), type,
                   Chombo_MPI::comm, &m_mpi_requests.back());
#else
    MPI_Alltoallv(sendbuf, m_query.countsPtr(), m_query.displsPtr(), type,
                  recvbuf, m_answer.countsPtr(), m_answer.displsPtr(), type,
                  Chombo_MPI::comm);
#endif
}

inline void MPIContext::asyncExchangeAnswer(void *sendbuf, void *recvbuf,
                                            MPI_Datatype type)
{
    CH_assert(m_async_active);
    MPI_Request req;
    m_mpi_requests.push_back(req);

#if MPI_VERSION >= 3 && !defined(OPEN_MPI)
    MPI_Ialltoallv(sendbuf, m_answer.countsPtr(), m_answer.displsPtr(), type,
                   recvbuf, m_query.countsPtr(), m_query.displsPtr(), type,
                   Chombo_MPI::comm, &m_mpi_requests.back());
#else
    MPI_Alltoallv(sendbuf, m_answer.countsPtr(), m_answer.displsPtr(), type,
                  recvbuf, m_query.countsPtr(), m_query.displsPtr(), type,
                  Chombo_MPI::comm);
#endif
}

inline void MPIContext::asyncEnd()
{
    CH_assert(m_async_active);
    m_async_active = false;

#if MPI_VERSION >= 3 && !defined(OPEN_MPI)
    MPI_Waitall(m_mpi_requests.size(), &m_mpi_requests[0], MPI_STATUSES_IGNORE);
#endif

    m_mpi_requests.clear();
}
#endif /* ifdef CH_MPI */

#endif /* MPICONTEXT_IMPL_HPP_ */
