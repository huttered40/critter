
#ifndef __CRITTER_H__
#define __CRITTER_H__

#include "mpi.h"
#include <fstream>
#include <stdint.h>
#include <functional>
#include <map>
#include <cmath>
#include <assert.h>

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

    /* Running sums in order to calculate averages */
    double my_bytesSum, my_comm_timeSum, my_bar_timeSum, crit_bytesSum, crit_comm_timeSum, crit_bar_timeSum, crit_msgSum, crit_wrdSum;

    /* \brief function for cost model of collective, takes (msg_size_in_bytes, number_processors) and returns (latency_cost, bandwidth_cost) */
    std::function< std::pair<double,double>(int64_t,int) > cost_func;

    /* \brief time when start() was last called, set to -1.0 initially and after stop() */
    double last_start_time;
    /* \brief cm with which start() was last called */
    MPI_Comm last_cm;
    /* \brief nbr_pe with which start() was last called */
    int last_nbr_pe;
    /* \brief nbr_pe2 with which start() was last called */
    int last_nbr_pe2;

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
     * \param[in] nbe_pe2 second neighbor processor (only used for p2p routines)
     * \param[in] is_async whether the call is asynchronous (used only for p2p routines)
     */
    void start(int64_t nelem=1, MPI_Datatype t=MPI_CHAR, MPI_Comm cm=MPI_COMM_WORLD, int nbr_pe=-1, int nbr_pe2=-1, bool is_async=false);

    /**
     * \brief stop timer, record time (use last_*, ensure last_start_time != -1., set last_start_time to -1), performs barrier over last cm
     */
    void stop();

    /**
     * \brief computes max critical path costs over given communicator (used internally and can be used at end of execution
     * \param[in] cm communicator over which we want to get the maximum cost
     * \param[in] nbe_pe neighbor processor (only used for p2p routines)
     * \param[in] nbe_pe2 second neighbor processor (only used for p2p routines)
     */
    void compute_max_crit(MPI_Comm cm=MPI_COMM_WORLD, int nbr_pe=-1, int nbr_pe2=-1);

    void compute_avg_crit_update();
    
    /**
     * \brief prints timer data for critical path measurements
     */
    void print_crit(std::ofstream& ptr);

    /**
     * \brief prints averaged timer data over all 'numIter' iterations for critical path measurements
     */
    void print_crit_avg(std::ofstream& fptr, int numIter);

    /**
     * \brief prints timer data for local measurements
     */
    void print_local();

    /**
     * \brief evaluates communication cost model as specifed by cost_func
     * \return pair (latency cost, bandwidth cost)
     */
    std::pair<double,double> get_crit_cost();

    /**
     * \brief common initialization of variables for construtors
     */
    void init();

    /**
     * \brief common initialization of variables ultimately used to find the average of each critter routine
     */
    void initSums();
/*
  private:
    void init();
*/
};
#define NUM_CRITTERS 16

extern
Critter * critter_list[NUM_CRITTERS];

/* \brief request/Critter dictionary for asynchronous messages */
extern
std::map<MPI_Request, Critter*> critter_req;

extern
Critter MPI_Barrier_critter, 
        MPI_Bcast_critter, 
        MPI_Reduce_critter, 
        MPI_Allreduce_critter, 
        MPI_Gather_critter, 
        MPI_Gatherv_critter, 
        MPI_Allgather_critter, 
        MPI_Allgatherv_critter, 
        MPI_Scatter_critter, 
        MPI_Scatterv_critter, 
        MPI_Reduce_scatter_critter, 
        MPI_Alltoall_critter, 
        MPI_Alltoallv_critter, 
        MPI_Send_critter, 
        MPI_Recv_critter, 
        MPI_Isend_critter, 
        MPI_Irecv_critter, 
        MPI_Sendrecv_critter; 

void compute_all_max_crit(MPI_Comm cm, int nbr_pe, int nbr_pe2);
void compute_all_avg_crit_updates();

// Instead of printing out each Critter for each iteration individually, I will save them for each iteration, print out the iteration, and then clear before next iteration
extern std::map<std::string,std::tuple<double,double,double,double,double> > saveCritterInfo;

#define Critter_Clear()                                   \
   do {                                                  \
    assert(critter_req.size() == 0);                     \
    for (int i=0; i<NUM_CRITTERS; i++){                  \
      critter_list[i]->init();                   \
    }                                   \
  } while (0)

#define Critter_Print(ARG1, ARG2)            \
   do {                                                  \
    assert(critter_req.size() == 0);                     \
    int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);  \
    if (myrank == 0)					 \
    { printf("\nCRITTER\n");}				 \
    compute_all_max_crit(MPI_COMM_WORLD,-1,-1);          \
    compute_all_avg_crit_updates();			 \
    if (myrank == 0) {                                   \
      printf("\t\t comm_bytes\t comm_time\t bar_time "); \
      printf("\t msg_cost \t wrd_cost\n");               \
    }                                                    \
    for (int i=0; i<NUM_CRITTERS; i++){                  \
      if (myrank == 0) {                                 \
        critter_list[i]->print_crit(ARG1);      \
      }                                                  \
    }							 \
    							\
    if (ARG2 == 0)					\
    {							\
      ARG1 << "Input";				\
      for (auto& it : saveCritterInfo)			 \
      {							\
        ARG1 << "\t" << it.first;	\
      }							\
      ARG1 << "\n";					\
    }							\
    ARG1 << ARG2;				\
    for (auto& it : saveCritterInfo)			 \
    {							\
      ARG1 << "\t" << std::get<0>(it.second);	\
    }							\
    ARG1 << "\n";							\
    ARG1 << ARG2;				\
    for (auto& it : saveCritterInfo)			 \
    {							\
      ARG1 << "\t" << std::get<1>(it.second);	\
    }							\
    ARG1 << "\n";							\
    ARG1 << ARG2;				\
    for (auto& it : saveCritterInfo)			 \
    {							\
      ARG1 << "\t" << std::get<2>(it.second);	\
    }							\
    ARG1 << "\n";							\
    ARG1 << ARG2;				\
    for (auto& it : saveCritterInfo)			 \
    {							\
      ARG1 << "\t" << std::get<3>(it.second);	\
    }							\
    ARG1 << "\n";							\
    ARG1 << ARG2;				\
    for (auto& it : saveCritterInfo)			 \
    {							\
      ARG1 << "\t" << std::get<4>(it.second);	\
    }							\
    ARG1 << "\n";							\
    saveCritterInfo.clear();				\
  } while (0)

/*
#define MPI_Finalize()                                   \
   do {                                                  \
    assert(critter_req.size() == 0);                     \
    int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);  \
    compute_all_max_crit(MPI_COMM_WORLD,-1,-1);          \
    if (myrank == 0) {                                   \
      printf("\t\t comm_bytes\t comm_time\t bar_time "); \
      printf("\t msg_cost \t wrd_cost\n");               \
    }                                                    \
    for (int i=0; i<NUM_CRITTERS; i++){                  \
      if (myrank == 0) {                                 \
        critter_list[i]->print_crit();                   \
      }                                                  \
    } PMPI_Finalize();                                   \
  } while (0)
*/

#define MPI_Barrier(cm)                            \
  do { MPI_Barrier_critter.start(0, MPI_CHAR, cm); \
    PMPI_Barrier(cm);                              \
    MPI_Barrier_critter.stop();                    \
  } while (0)

#define MPI_Bcast(buf, nelem, t, root, cm)    \
  do { MPI_Bcast_critter.start(nelem, t, cm); \
    PMPI_Bcast(buf, nelem, t, root, cm);      \
    MPI_Bcast_critter.stop();                 \
  } while (0)

#define MPI_Allreduce(sbuf, rbuf, nelem, t, op, cm) \
  do { MPI_Allreduce_critter.start(nelem, t, cm);   \
    PMPI_Allreduce(sbuf, rbuf, nelem, t, op, cm);   \
    MPI_Allreduce_critter.stop();                   \
  } while (0)

#define MPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm) \
  do { MPI_Reduce_critter.start(nelem, t, cm);         \
    PMPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm);   \
    MPI_Reduce_critter.stop();                         \
  } while (0)

#define MPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm)                                    \
  do { assert(rt==st); MPI_Scatter_critter.start(std::max((int64_t)scount,(int64_t)rcount), st, cm); \
    PMPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm);                                      \
    MPI_Scatter_critter.stop();                                                                      \
  } while (0)

#define MPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm)                                    \
  do { assert(rt==st); MPI_Gather_critter.start(std::max((int64_t)scount,(int64_t)rcount), st, cm); \
    PMPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm);                                      \
    MPI_Gather_critter.stop();                                                                      \
  } while (0)

#define MPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm)                                          \
  do { assert(rt==st); MPI_Allgather_critter.start(std::max((int64_t)scount,(int64_t)rcount), st, cm); \
    PMPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm);                                            \
    MPI_Allgather_critter.stop();                                                                      \
  } while (0)

#define MPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm) \
  do { int64_t tot_recv=0;                                 \
    int p; MPI_Comm_size(cm, &p);                          \
    for (int i=0; i<p; i++){ tot_recv += rcounts[i]; }     \
    MPI_Reduce_scatter_critter.start(tot_recv, t, cm);     \
    PMPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm);   \
    MPI_Reduce_scatter_critter.stop();                     \
  } while (0)

#define MPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm)                                          \
  do { assert(rt==st); MPI_Alltoall_critter.start(std::max((int64_t)scount,(int64_t)rcount), st, cm); \
    PMPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm);                                            \
    MPI_Alltoall_critter.stop();                                                                      \
  } while (0)

#define MPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm)     \
  do { assert(rt==st);                                                        \
    int64_t tot_recv=0;                                                       \
    int p; MPI_Comm_size(cm, &p);                                             \
    for (int i=0; i<p; i++){ tot_recv += rcounts[i]; }                        \
    MPI_Allgatherv_critter.start(std::max((int64_t)scount,tot_recv), st, cm); \
    PMPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm);       \
    MPI_Allgatherv_critter.stop();                                            \
  } while (0)


#define MPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm)  \
  do { assert(rt==st);                                                        \
    int64_t tot_recv=0;                                                       \
    int r, p; MPI_Comm_rank(cm, &r); MPI_Comm_size(cm, &p);                   \
    if (r == root) for (int i=0; i<p; i++){ tot_recv += ((int*)rcounts)[i]; } \
    MPI_Gatherv_critter.start(std::max((int64_t)scount,tot_recv), st, cm);    \
    PMPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm);    \
    MPI_Gatherv_critter.stop();                                               \
  } while (0)

#define MPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm)  \
  do { assert(rt==st);                                                        \
    int64_t tot_send=0;                                                       \
    int r, p; MPI_Comm_rank(cm, &r); MPI_Comm_size(cm, &p);                   \
    if (r == root) for (int i=0; i<p; i++){ tot_send += ((int*)scounts)[i]; } \
    MPI_Scatterv_critter.start(std::max(tot_send,(int64_t)rcount), st, cm);   \
    PMPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm);    \
    MPI_Scatterv_critter.stop();                                              \
  } while (0)

#define MPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm) \
  do { assert(rt==st);                                                             \
    int64_t tot_send=0, tot_recv=0;                                                \
    int p; MPI_Comm_size(cm, &p);                                                  \
    for (int i=0; i<p; i++){ tot_send += scounts[i]; tot_recv += rcounts[i]; }     \
    MPI_Alltoallv_critter.start(std::max(tot_send,tot_recv), st, cm);              \
    PMPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm);   \
    MPI_Alltoallv_critter.stop();                                                  \
  } while (0)


#define MPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status) \
  do { assert(st == rt);                                                                \
    MPI_Sendrecv_critter.start(std::max(scnt,rcnt), st, cm, dest, src);                 \
    PMPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status);   \
    MPI_Sendrecv_critter.stop();                                                        \
  } while (0)

#if 0
#define MPI_Send(buf, nelem, t, dest, tag, cm)                                            \
  do { MPI_Send_critter.start(nelem, t, cm, dest);                                        \
    PMPI_Send(buf, nelem, t, dest, tag, cm);                                      \
    MPI_Send_critter.stop(); \
  } while (0)

#define MPI_Recv(buf, nelem, t, src, tag, cm, status)                                            \
  do { MPI_Recv_critter.start(nelem, t, cm, src);                                        \
    PMPI_Recv(buf, nelem, t, src, tag, cm, status);                                      \
    MPI_Recv_critter.stop(); \
  } while (0)

#define MPI_Irecv(buf, nelem, t, src, tag, cm, req)                                            \
  do { MPI_Irecv_critter.start(nelem, t, cm, src, -1, 1);                                        \
    PMPI_Irecv(buf, nelem, t, src, tag, cm, req);                                      \
    critter_req[*req] = &MPI_Irecv_critter; \
  } while (0)

#define MPI_Isend(buf, nelem, t, dest, tag, cm, req)                                            \
  do { MPI_Isend_critter.start(nelem, t, cm, dest, -1, 1);                                        \
    PMPI_Isend(buf, nelem, t, dest, tag, cm, req);                                      \
    critter_req[*req] = &MPI_Isend_critter; \
  } while (0)

#define MPI_Wait(req, stat) \
  do { std::map<MPI_Request, Critter*>::iterator it = critter_req.find(*req); \
    if (it == critter_req.end()) *(int*)NULL = 1; \
    assert(it != critter_req.end()); \
    PMPI_Wait(req, stat); \
    it->second->stop(); critter_req.erase(it);  \
  } while (0)

#define MPI_Waitany(cnt, reqs, indx, stat) \
  do { PMPI_Waitany(cnt, reqs, indx, stat); \
    std::map<MPI_Request, Critter*>::iterator it = critter_req.find((reqs)[*(indx)]); \
    if (it != critter_req.end()) { it->second->stop(); critter_req.erase(it); } \
  } while (0)

#define MPI_Waitall(cnt, reqs, stats) \
  do { int __indx; MPI_Status __stat; for (int i=0; i<cnt; i++){ MPI_Waitany(cnt, reqs, &__indx, &__stat); if ((MPI_Status*)stats != (MPI_Status*)MPI_STATUSES_IGNORE) ((MPI_Status*)stats)[__indx] = __stat; } \
  } while (0)
#endif
#endif
