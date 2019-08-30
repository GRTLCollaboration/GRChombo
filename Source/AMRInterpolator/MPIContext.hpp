/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MPICONTEXT_HPP_
#define MPICONTEXT_HPP_

#include "MPILayout.hpp"
#include <vector>

class MPIContext
{
  public:
    const int m_num_process;
    const int m_rank;

    MPIContext();

    // Getters
    inline int queryCount(int rank);
    inline int totalQueryCount();
    inline int answerCount(int rank);
    inline int totalAnswerCount();
    inline int queryDispl(int rank);
    inline int answerDispl(int rank);

    // Setters
    inline void setQueryCount(int rank, int count);
    inline void incrementQueryCount(int rank);
    inline void clearQueryCounts();

    void exchangeLayout();

#ifdef CH_MPI
    // MPI asynchronous comms
    inline void asyncBegin();
    inline void asyncExchangeQuery(void *sendbuf, void *recvbuf,
                                   MPI_Datatype type);
    inline void asyncExchangeAnswer(void *sendbuf, void *recvbuf,
                                    MPI_Datatype type);
    inline void asyncEnd();
#endif

    // MPI utils
    static int comm_size();
    static int comm_rank();

  private:
    MPILayout m_query;
    MPILayout m_answer;

    bool m_query_dirty;

    bool m_async_active;
#ifdef CH_MPI
    std::vector<MPI_Request> m_mpi_requests;
#endif
};

inline int MPIContext::comm_size()
{
    int out = 1;
#ifdef CH_MPI
    MPI_Comm_size(Chombo_MPI::comm, &out);
#endif
    return out;
}

inline int MPIContext::comm_rank()
{
    int out = 0;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &out);
#endif
    return out;
}

#include "MPIContext.impl.hpp"

#endif /* MPICONTEXT_HPP_ */
