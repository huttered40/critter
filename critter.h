
#ifndef __CRITTER_H__
#define __CRITTER_H__

#include "mpi.h"
#include <stdint.h>
#include <functional>
#include <cmath>

class Critter {
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

    /* \brief comm cost #messages term*/
    double crit_msg;
    /* \brief comm cost #words term */
    double crit_wrd;
    /* \brief name of collective */
    char * name;
    
    /* \brief function for cost model of collective, takes (msg_size_in_bytes, number_processors) and returns (latency_cost, bandwidth_cost) */
    std::function< std::pair<double,double>(int64_t,int) > cost_func;

    /* \brief time when start() was last called, set to -1.0 initially and after stop() */
    double last_start_time;
    /* \brief cm with which start() was last called */
    MPI_Comm last_cm;
    /* \brief nbr_pe with which start() was last called */
    int last_nbr_pe;
    /* \brief nbr_pe2 with which start() was last called */
    int nbr_pe2;

    /**
     * \brief timer constructor, initializes vars
     * \param[in] name symbol name of MPI routine
     * \param[in] function for cost model of collective, takes (msg_size_in_bytes, number_processors) and returns (latency_cost, bandwidth_cost) 
     */
    Critter(char const * name,
            std::function< std::pair<double,double>(int64_t,int) > 
              cost_func = [](int64_t n, int p){ 
                return std::pair<double,double>(1.,n); 
              });
 
    /**
     * \brief timer copy constructor, copies name
     * \param[in] t other timer
     */
    Critter(Critter const & t);
    
   
    /**
     * \brief timer destructor, frees name
     */ 
    ~Critter();

    /**
     * \brief starts timer for MPI call with nelem elements of type t over communicator cm, performs barrier over cm
     * \param[in] name symbol name of MPI routine
     * \param[in] cm MPI_Communicator on which MPI routine is called
     * \param[in] nbe_pe neighbor processor (only used for p2p routines)
     */
    void start(int64_t nelem=1, MPI_Datatype t=MPI_CHAR, MPI_Comm cm=MPI_COMM_WORLD, int nbr_pe=-1, int nbr_pe2=-1);

    /**
     * \brief stop timer, record time (use last_*, ensure last_start_time != -1., set last_start_time to -1), performs barrier over last cm
     */
    void stop();

    /**
     * \brief computes max critical path costs over given communicator (used internally and can be used at end of execution
     * \param[in] cm communicator over which we want to get the maximum cost
     */
    void compute_max_crit(MPI_Comm cm=MPI_COMM_WORLD);
    
    /**
     * \brief prints timer data for critical path measurements
     */
    void print_crit();

    /**
     * \brief prints timer data for local measurements
     */
    void print_local();

    /**
     * \brief evaluates communication cost model as specifed by cost_func
     * \return pair (latency cost, bandwidth cost)
     */
    std::pair<double,double> get_crit_cost();

  private:
    /**
     * \brief common initialization of variables for construtors
     */
    void init();
};

static
Critter MPI_Bcast_critter("MPI_Bcast", 
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }), 
        MPI_Reduce_critter("MPI_Reduce", 
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }), 
        MPI_Allreduce_critter("MPI_Allreduce",
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }), 
        MPI_Gather_critter("MPI_Gather",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 

                          }), 
        MPI_Allgather_critter("MPI_Allgather",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Scatter_critter("MPI_Scatter",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 

        MPI_Reduce_scatter_critter("MPI_Reduce_scatter",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Alltoall_critter("MPI_Alltoall",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n); 
                          }), 
        MPI_Alltoallv_critter("MPI_Alltoallv",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n); 
                          }), 
        MPI_Sendrecv_critter("MPI_Sendrecv",
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }); 

#define NUM_CRITTERS 10

static
Critter * critter_list[NUM_CRITTERS] = {
        &MPI_Bcast_critter,
        &MPI_Reduce_critter,
        &MPI_Allreduce_critter,
        &MPI_Scatter_critter,
        &MPI_Gather_critter,
        &MPI_Allgather_critter,
        &MPI_Reduce_scatter_critter,
        &MPI_Alltoall_critter,
        &MPI_Alltoallv_critter,
        &MPI_Sendrecv_critter };

#define MPI_Finalize() \
   do { \
    int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank); \
    for (int i=0; i<NUM_CRITTERS; i++){ \
      critter_list[i]->compute_max_crit(MPI_COMM_WORLD); \
      if (myrank == 0) { \
        critter_list[i]->print_crit(); \
        critter_list[i]->print_local(); \
      } \
    } PMPI_Finalize(); } while (0)

#define MPI_Bcast(buf, nelem, t, root, cm)                                            \
  { MPI_Bcast_critter.start(nelem, t, cm);                                        \
    PMPI_Bcast(buf, nelem, t, root, cm);                                      \
    MPI_Bcast_critter.stop(); }

#define MPI_Allreduce(sbuf, rbuf, nelem, t, op, cm)                                            \
  { MPI_Allreduce_critter.start(nelem, t, cm);                                        \
    PMPI_Allreduce(sbuf, rbuf, nelem, t, op, cm);                                      \
    MPI_Allreduce_critter.stop(); }

#define MPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm)                                            \
  { MPI_Reduce_critter.start(nelem, t, cm);                                        \
    PMPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm);                                      \
    MPI_Reduce_critter.stop(); }

#define MPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm)                                            \
  { assert(rt==st); MPI_Scatter_critter.start(std::max(scount,rcount), st, cm);                                        \
    PMPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm);                                      \
    MPI_Scatter_critter.stop(); }

#define MPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm)                                            \
  { assert(rt==st); MPI_Gather_critter.start(std::max(scount,rcount), st, cm);                                        \
    PMPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm);                                      \
    MPI_Gather_critter.stop(); }

#define MPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm)                                            \
  { assert(rt==st); MPI_Allgather_critter.start(std::max(scount,rcount), st, cm);                                        \
    PMPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm);                                      \
    MPI_Allgather_critter.stop(); }

#define MPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm)                                            \
  { int64_t tot_recv=0; \
    int p; MPI_Comm_size(cm, &p); \
    for (int i=0; i<p; i++){ tot_recv += rcounts[i]; } \\
    MPI_Reduce_scatter_critter.start(tot_recv, t, cm);                                        \
    PMPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm);                                      \
    MPI_Reduce_scatter_critter.stop(); }

#define MPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm)                                            \
  { assert(rt==st); MPI_Alltoall_critter.start(std::max(scount,rcount), st, cm);                                        \
    PMPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm);                                      \
    MPI_Alltoall_critter.stop(); }

#define MPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm)                                            \
  { assert(rt==st); \
    int64_t tot_send=0, tot_recv=0; \
    int p; MPI_Comm_size(cm, &p); \
    for (int i=0; i<p; i++){ tot_send += scounts[i]; tot_recv += rcounts[i]; } \\
    MPI_Alltoallv_critter.start(std::max(tot_send,tot_recv), st, cm);                                        \
    PMPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm);                                      \
    MPI_Alltoallv_critter.stop(); }


//#define MPI_Reduce(...)                                           \
//  { Critter __t("MPI_Reduce");                                  \
//              __t.start();                                        \
//    PMPI_Reduce(__VA_ARGS__);                                     \
//              __t.stop(); }
//#define MPI_Wait(...)                                             \
//  { Critter __t("MPI_Wait");                                    \
//              __t.start();                                        \
//    PMPI_Wait(__VA_ARGS__);                                       \
//              __t.stop(); }
//#define MPI_Send(...)                                             \
//  { Critter __t("MPI_Send");                                    \
//              __t.start();                                        \
//    PMPI_Send(__VA_ARGS__);                                       \
//              __t.stop(); }
//#define MPI_Allreduce(...)                                        \
//  { Critter __t("MPI_Allreduce");                               \
//              __t.start();                                        \
//    PMPI_Allreduce(__VA_ARGS__);                                  \
//              __t.stop(); }
//#define MPI_Allgather(...)                                        \
//  { Critter __t("MPI_Allgather");                               \
//              __t.start();                                        \
//    PMPI_Allgather(__VA_ARGS__);                                  \
//              __t.stop(); }
//#define MPI_Scatter(...)                                          \
//  { Critter __t("MPI_Scatter");                                 \
//              __t.start();                                        \
//    PMPI_Scatter(__VA_ARGS__);                                    \
//              __t.stop(); }
//#define MPI_Alltoall(...)                                         \
//  { Critter __t("MPI_Alltoall");                                \
//              __t.start();                                        \
//    PMPI_Alltoall(__VA_ARGS__);                                   \
//              __t.stop(); }
//#define MPI_Alltoallv(...)                                        \
//  { Critter __t("MPI_Alltoallv");                               \
//              __t.start();                                        \
//    PMPI_Alltoallv(__VA_ARGS__);                                  \
//              __t.stop(); }
//#define MPI_Gatherv(...)                                          \
//  { Critter __t("MPI_Gatherv");                                 \
//              __t.start();                                        \
//    PMPI_Gatherv(__VA_ARGS__);                                    \
//              __t.stop(); }
//#define MPI_Scatterv(...)                                         \
//  { Critter __t("MPI_Scatterv");                                \
//              __t.start();                                        \
//   PMPI_Scatterv(__VA_ARGS__);                                    \
//              __t.stop(); }
//#define MPI_Waitall(...)                                          \
//  { Critter __t("MPI_Waitall");                                 \
//              __t.start();                                        \
//    PMPI_Waitall(__VA_ARGS__);                                    \
//              __t.stop(); }
//#define MPI_Barrier(...)                                          \
//  { Critter __t("MPI_Barrier");                                 \
//              __t.start();                                        \
//    PMPI_Barrier(__VA_ARGS__);                                    \
//              __t.stop(); }

#endif
