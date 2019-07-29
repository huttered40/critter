
#ifndef __CRITTER_H__
#define __CRITTER_H__

#include "mpi.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <stdint.h>
#include <functional>
#include <map>
#include <cmath>
#include <assert.h>

namespace Critter{

class _Critter {
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
    std::string name;

     /* \brief duration of computation time for each call made locally,
 *        Used to save the local computation time between _Critter methods, so that we can synchronize it after communication */
    double my_comp_time;

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
    _Critter(std::string name,
            std::function< std::pair<double,double>(int64_t,int) > 
              cost_func = [](int64_t n, int p){ 
                return std::pair<double,double>(1.,n); 
              });
 
    /**
     * \brief timer copy constructor, copies name
     * \param[in] t other timer
     */
    _Critter(_Critter const & t);
    
   
    /**
     * \brief timer destructor, frees name
     */ 
    ~_Critter();

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
    void print_crit(std::ofstream& ptr, std::string name);

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
#define NUM_CRITTERS 17

extern _Critter * critter_list[NUM_CRITTERS];

/* \brief request/Critter dictionary for asynchronous messages */
extern std::map<MPI_Request,_Critter*> critter_req;

extern double totalCritComputationTime;
extern double curComputationTimer;
extern double totalOverlapTime;			// Updated at each BSP step
extern double totalCommunicationTime;
extern double totalIdleTime;
// Instead of printing out each Critter for each iteration individually, I will save them for each iteration, print out the iteration, and then clear before next iteration
extern std::map<std::string,std::tuple<double,double,double,double,double,double,double,double>> saveCritterInfo;
extern std::map<std::string,std::vector<std::string>> AlgCritters;

extern _Critter MPI_Barrier_critter, 
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
         MPI_Sendrecv_critter, 
         MPI_Sendrecv_replace_critter; 


void FillAlgCritterList();
bool InAlgCritterList(std::string AlgName, std::string CritterName);
void compute_all_max_crit(MPI_Comm cm, int nbr_pe, int nbr_pe2);
void compute_all_avg_crit_updates();
void reset();
void print(std::ofstream& Stream, std::string AlgName, int NumPEs, size_t NumInputs, size_t* Inputs, const char** InputNames);
}

/*
#define MPI_Finalize()\
   do {\
    assert(critter_req.size() == 0);\
    int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);\
    compute_all_max_crit(MPI_COMM_WORLD,-1,-1);\
    if (myrank == 0) {\
    }\
    for (int i=0; i<NUM_CRITTERS; i++){\
      if (myrank == 0) {\
        critter_list[i]->print_crit();\
      }\
    } PMPI_Finalize();\
  } while (0)
*/

#define MPI_Barrier(cm)\
  do { Critter::MPI_Barrier_critter.start(0, MPI_CHAR, cm);\
    PMPI_Barrier(cm);\
    Critter::MPI_Barrier_critter.stop();\
  } while (0)

#define MPI_Bcast(buf, nelem, t, root, cm)\
  do { Critter::MPI_Bcast_critter.start(nelem, t, cm);\
    PMPI_Bcast(buf, nelem, t, root, cm);\
    Critter::MPI_Bcast_critter.stop();\
  } while (0)

#define MPI_Allreduce(sbuf, rbuf, nelem, t, op, cm)\
  do { Critter::MPI_Allreduce_critter.start(nelem, t, cm);\
    PMPI_Allreduce(sbuf, rbuf, nelem, t, op, cm);\
    Critter::MPI_Allreduce_critter.stop();\
  } while (0)

#define MPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm)\
  do { Critter::MPI_Reduce_critter.start(nelem, t, cm);\
    PMPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm);\
    Critter::MPI_Reduce_critter.stop();\
  } while (0)

#define MPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
  do { assert(rt==st); Critter::MPI_Scatter_critter.start(std::max((int64_t)scount,(int64_t)rcount), st, cm);\
    PMPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
    Critter::MPI_Scatter_critter.stop();\
  } while (0)

#define MPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
  do { assert(rt==st);\
    int pSize; MPI_Comm_size(cm, &pSize);\
    int64_t recvBufferSize = std::max((int64_t)scount,(int64_t)rcount) * pSize;\
    Critter::MPI_Gather_critter.start(recvBufferSize, st, cm);\
    PMPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
    Critter::MPI_Gather_critter.stop();\
  } while (0)

#define MPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm)\
  do { assert(rt==st);\
    int pSize; MPI_Comm_size(cm, &pSize);\
    int64_t recvBufferSize = std::max((int64_t)scount,(int64_t)rcount) * pSize;\
    Critter::MPI_Allgather_critter.start(recvBufferSize, st, cm);\
    PMPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm);\
    Critter::MPI_Allgather_critter.stop();\
  } while (0)

#define MPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm)\
  do { int64_t tot_recv=0;\
    int p; MPI_Comm_size(cm, &p);\
    for (int i=0; i<p; i++){ tot_recv += rcounts[i]; }\
    Critter::MPI_Reduce_scatter_critter.start(tot_recv, t, cm);\
    PMPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm);\
    Critter::MPI_Reduce_scatter_critter.stop();\
  } while (0)

#define MPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm)\
  do { assert(rt==st); Critter::MPI_Alltoall_critter.start(std::max((int64_t)scount,(int64_t)rcount), st, cm);\
    PMPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm);\
    Critter::MPI_Alltoall_critter.stop();\
  } while (0)

#define MPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm)\
  do { assert(rt==st);\
    int64_t tot_recv=0;\
    int p; MPI_Comm_size(cm, &p);\
    for (int i=0; i<p; i++){ tot_recv += rcounts[i]; }\
    Critter::MPI_Allgatherv_critter.start(std::max((int64_t)scount,tot_recv), st, cm);\
    PMPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm);\
    Critter::MPI_Allgatherv_critter.stop();\
  } while (0)


#define MPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm)\
  do { assert(rt==st);\
    int64_t tot_recv=0;\
    int r, p; MPI_Comm_rank(cm, &r); MPI_Comm_size(cm, &p);\
    if (r == root) for (int i=0; i<p; i++){ tot_recv += ((int*)rcounts)[i]; }\
    Critter::MPI_Gatherv_critter.start(std::max((int64_t)scount,tot_recv), st, cm);\
    PMPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm);\
    Critter::MPI_Gatherv_critter.stop();\
  } while (0)

#define MPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm)  \
  do { assert(rt==st);                                                        \
    int64_t tot_send=0;                                                       \
    int r, p; MPI_Comm_rank(cm, &r); MPI_Comm_size(cm, &p);                   \
    if (r == root) for (int i=0; i<p; i++){ tot_send += ((int*)scounts)[i]; } \
    Critter::MPI_Scatterv_critter.start(std::max(tot_send,(int64_t)rcount), st, cm);   \
    PMPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm);    \
    Critter::MPI_Scatterv_critter.stop();\
  } while (0)

#define MPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm)\
  do { assert(rt==st);\
    int64_t tot_send=0, tot_recv=0;\
    int p; MPI_Comm_size(cm, &p);\
    for (int i=0; i<p; i++){ tot_send += scounts[i]; tot_recv += rcounts[i]; }\
    Critter::MPI_Alltoallv_critter.start(std::max(tot_send,tot_recv), st, cm);\
    PMPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm);\
    Critter::MPI_Alltoallv_critter.stop();\
  } while (0)


#define MPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status)\
  do { assert(st == rt);\
    Critter::MPI_Sendrecv_critter.start(std::max(scnt,rcnt), st, cm, dest, src);\
    PMPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status);\
    Critter::MPI_Sendrecv_critter.stop();\
  } while (0)

#define MPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status)\
    do {\
    Critter::MPI_Sendrecv_replace_critter.start(scnt, st, cm, dest, src);\
    PMPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status);\
    Critter::MPI_Sendrecv_replace_critter.stop();\
  } while (0)

#define MPI_Comm_split(comm1, arg2, arg3, comm2)\
    do {\
    volatile double curTime = MPI_Wtime();\
    double localCompTime = curTime - Critter::curComputationTimer;\
    curTime = MPI_Wtime();\
    PMPI_Comm_split(comm1, arg2, arg3, comm2);\
    double localCommTime = MPI_Wtime() - curTime;\
    double localTotalTime = localCompTime + localCommTime;\
    std::vector<double> critterVec(3);\
    std::vector<double> localVec(3);\
    localVec[0] = localCompTime; localVec[1] = localCommTime; localVec[2] = localTotalTime;\
    PMPI_Allreduce(&localVec[0], &critterVec[0], 3, MPI_DOUBLE, MPI_MAX, comm1);\
    Critter::totalCritComputationTime += critterVec[0];\
    Critter::totalCommunicationTime += critterVec[1];\
    Critter::totalOverlapTime += (critterVec[0] + critterVec[1] - critterVec[2]);\
    Critter::curComputationTimer = MPI_Wtime();\
  } while (0)

#define MPI_Send(buf, nelem, t, dest, tag, cm)\
  do { Critter::MPI_Send_critter.start(nelem, t, cm, dest);\
    PMPI_Send(buf, nelem, t, dest, tag, cm);\
    Critter::MPI_Send_critter.stop();\
  } while (0)

#define MPI_Recv(buf, nelem, t, src, tag, cm, status)\
  do { Critter::MPI_Recv_critter.start(nelem, t, cm, src);\
    PMPI_Recv(buf, nelem, t, src, tag, cm, status);\
    Critter::MPI_Recv_critter.stop();\
  } while (0)

#define MPI_Irecv(buf, nelem, t, src, tag, cm, req)\
  do { Critter::MPI_Irecv_critter.start(nelem, t, cm, src, -1, 1);\
    PMPI_Irecv(buf, nelem, t, src, tag, cm, req);\
    Critter::critter_req[*req] = &Critter::MPI_Irecv_critter;\
  } while (0)

#define MPI_Isend(buf, nelem, t, dest, tag, cm, req)\
  do { Critter::MPI_Isend_critter.start(nelem, t, cm, dest, -1, 1);\
    PMPI_Isend(buf, nelem, t, dest, tag, cm, req);\
    Critter::critter_req[*req] = &Critter::MPI_Isend_critter;\
  } while (0)

#define MPI_Wait(req, stat)\
  do { std::map<MPI_Request, Critter*>::iterator it = Critter::critter_req.find(*req);\
    if (it == Critter::critter_req.end()) *(int*)NULL = 1;\
    assert(it != Critter::critter_req.end());\
    PMPI_Wait(req, stat);\
    it->second->stop(); Critter::critter_req.erase(it);\
  } while (0)

#define MPI_Waitany(cnt, reqs, indx, stat)\
  do { PMPI_Waitany(cnt, reqs, indx, stat);\
    std::map<MPI_Request, Critter*>::iterator it = Critter::critter_req.find((reqs)[*(indx)]);\
    if (it != Critter::critter_req.end()) { it->second->stop(); Critter::critter_req.erase(it); }\
  } while (0)

#define MPI_Waitall(cnt, reqs, stats)\
  do { int __indx; MPI_Status __stat; for (int i=0; i<cnt; i++){ MPI_Waitany(cnt, reqs, &__indx, &__stat); if ((MPI_Status*)stats != (MPI_Status*)MPI_STATUSES_IGNORE) ((MPI_Status*)stats)[__indx] = __stat; }\
  } while (0)
#endif