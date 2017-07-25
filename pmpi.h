
#ifndef __PMPI_H__
#define __PMPI_H__

#include "mpi.h"
#include <ctf.hpp>
namespace CTF {
  void set_context(MPI_Comm);
}

class PMPI_Timer {
  public: 
    /* \brief number of bytes max(sent,recv'ed) for each call made locally */
    double my_bytes;
    /* \brief duration of communication time for each call made locally */
    double my_comm_time;
    /* \brief duration of idle time for each call made locally */
    double my_bar_time;
    /* \brief number of bytes max(sent,recv'ed) for each call made along the critical path */
    double crit_bytes;
    /* \brief duration of communication time for each call made along the critical path */
    double crit_comm_time;
    /* \brief number of bytes max(sent,recv'ed) for each call made along the critical path */
    double crit_bar_time;

    /* \brief name of collective */
    char const * name;

    /* \brief time when start() was last called, set to -1.0 initially and after stop() */
    double last_start_time;
    /* \brief cm with which start() was last called */
    double last_cm;
    /* \brief nbr_pe with which start() was last called */
    double last_nbr_pe;

    /**
     * \brief timer constructor, initializes vars
     * \param[in] name symbol name of MPI routine
     */
    PMPI_Timer(char const * name);

    /**
     * \brief starts timer for MPI call with nbytes bytes over communicator cm
     * \param[in] name symbol name of MPI routine
     * \param[in] cm MPI_Communicator on which MPI routine is called
     * \param[in] nbe_pe neighbor processor (only used for p2p routines)
     */
    void start(int64_t nbytes=1, MPI_Comm cm=MPI_COMM_WORLD, int nbr_pe=-1);

    /**
     * \brief stop timer, record time (use last_*, ensure last_start_time != -1., set last_start_time to -1
     */
    void stop();

    /**
     * \brief computes max critical path costs over given communicator (used internally and can be used at end of execution 
     * \param[in] cm communicator over which we want to get the maximum cost
     */
    void compute_max_crit(MPI_Comm cm=MPI_COMM_WORLD);
};

PMPI_Timer bcast_timer, reduce_timer, allred_timer, sendrecv_timer;


//int64_t my_bcast_bytes;
//int64_t my_reduce_bytes;
//int64_t my_allreduce_bytes;
//int64_t my_allgather_bytes;
//int64_t my_allgatherv_bytes;
//int64_t my_alltoall_bytes;
//int64_t my_alltoallv_bytes;
//int64_t my_send_bytes;
//int64_t my_recv_bytes;
//int64_t my_sendrecv_bytes;

#ifdef PMPI
#ifdef PMPI
#define MPI_Bcast(...)                                            \
  { CTF::Timer __t("MPI_Bcast");                                   \
              __t.start();                                        \
    PMPI_Bcast(__VA_ARGS__);                                      \
              __t.stop(); }
#define MPI_Reduce(...)                                           \
  { CTF::Timer __t("MPI_Reduce");                                  \
              __t.start();                                        \
    PMPI_Reduce(__VA_ARGS__);                                     \
              __t.stop(); }
#define MPI_Wait(...)                                             \
  { CTF::Timer __t("MPI_Wait");                                    \
              __t.start();                                        \
    PMPI_Wait(__VA_ARGS__);                                       \
              __t.stop(); }
#define MPI_Send(...)                                             \
  { CTF::Timer __t("MPI_Send");                                    \
              __t.start();                                        \
    PMPI_Send(__VA_ARGS__);                                       \
              __t.stop(); }
#define MPI_Allreduce(...)                                        \
  { CTF::Timer __t("MPI_Allreduce");                               \
              __t.start();                                        \
    PMPI_Allreduce(__VA_ARGS__);                                  \
              __t.stop(); }
#define MPI_Allgather(...)                                        \
  { CTF::Timer __t("MPI_Allgather");                               \
              __t.start();                                        \
    PMPI_Allgather(__VA_ARGS__);                                  \
              __t.stop(); }
#define MPI_Scatter(...)                                          \
  { CTF::Timer __t("MPI_Scatter");                                 \
              __t.start();                                        \
    PMPI_Scatter(__VA_ARGS__);                                    \
              __t.stop(); }
#define MPI_Alltoall(...)                                         \
  { CTF::Timer __t("MPI_Alltoall");                                \
              __t.start();                                        \
    PMPI_Alltoall(__VA_ARGS__);                                   \
              __t.stop(); }
#define MPI_Alltoallv(...)                                        \
  { CTF::Timer __t("MPI_Alltoallv");                               \
              __t.start();                                        \
    PMPI_Alltoallv(__VA_ARGS__);                                  \
              __t.stop(); }
#define MPI_Gatherv(...)                                          \
  { CTF::Timer __t("MPI_Gatherv");                                 \
              __t.start();                                        \
    PMPI_Gatherv(__VA_ARGS__);                                    \
              __t.stop(); }
#define MPI_Scatterv(...)                                         \
  { CTF::Timer __t("MPI_Scatterv");                                \
              __t.start();                                        \
   PMPI_Scatterv(__VA_ARGS__);                                    \
              __t.stop(); }
#define MPI_Waitall(...)                                          \
  { CTF::Timer __t("MPI_Waitall");                                 \
              __t.start();                                        \
    PMPI_Waitall(__VA_ARGS__);                                    \
              __t.stop(); }
#define MPI_Barrier(...)                                          \
  { CTF::Timer __t("MPI_Barrier");                                 \
              __t.start();                                        \
    PMPI_Barrier(__VA_ARGS__);                                    \
              __t.stop(); }
#endif

#endif
