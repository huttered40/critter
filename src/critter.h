
#ifndef __CRITTER_H__
#define __CRITTER_H__

#include "mpi.h"
#include <fstream>
#include <iostream>
#include "iomanip"
#include <vector>
#include <stdint.h>
#include <functional>
#include <map>
#include <cmath>
#include <assert.h>

namespace critter{
namespace internal{
class _critter {
  public: 
    /* \brief name of collective */
    std::string name;
    /* \brief integer tag of collective */
    int tag;

    /* \brief number of bytes max(sent,recv'ed) for each call made locally */
    double my_bytes;
    /* \brief duration of communication time for each call made locally */
    double my_comm_time;
    /* \brief duration of idle time for each call made locally */
    double my_bar_time;
    /* \brief duration of communication time for each call made locally */
    double my_msg;
    /* \brief duration of idle time for each call made locally */
    double my_wrd;

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

    /* \brief duration of computation time for each call made locally,
     *   Used to save the local computation time between _critter methods, so that we can synchronize it after communication */
    double save_comp_time;

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
    _critter(std::string name, int tag,
            std::function< std::pair<double,double>(int64_t,int) > 
              cost_func = [](int64_t n, int p){ 
                return std::pair<double,double>(1.,n); 
              });
 
    /**
     * \brief timer copy constructor, copies name
     * \param[in] t other timer
     */
    _critter(_critter const & t);
    
   
    /**
     * \brief timer destructor, frees name
     */ 
    ~_critter();

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
    void get_crit_data(double* Container);
    void set_crit_data(double* Container);
    void get_avg_data(double* Container);
    void set_avg_data(double* Container, int CommSize);

    void compute_avg_update();
    
    /**
     * \brief appends timer data for critical path and average measurements to a global array
     */
    void save_crit();

    /**
     * \brief evaluates communication cost model as specifed by cost_func
     * \return pair (latency cost, bandwidth cost)
     */
    std::pair<double,double> get_crit_cost();

    /**
     * \brief common initialization of variables for construtors
     */
    void init();

};
//#define NUM_CRITTERS 17
constexpr auto NumCritters=18;

extern _critter * critter_list[NumCritters];

struct double_int{
  double_int(){first=0; second=0;}
  double_int(double one, int two){first=one; second=two;}
  double first; int second;
};
struct int_int_double{
  int_int_double(){first=0; second=0; third=0;}
  int_int_double(int one, int two, double three){first=one; second=two; third=three;}
  int first; int second; double third;
};

/* \brief request/critter dictionary for asynchronous messages */
extern std::map<MPI_Request,_critter*> critter_req;
extern std::string StreamName,StreamTrackName,FileName;
extern bool track,flag,IsFirstIter,IsWorldRoot,NeedNewLine;
extern std::ofstream Stream,StreamTrack;
extern double ComputationTimer;
extern std::vector<std::vector<int_int_double>> CritterPaths;
extern std::array<double,14> CritterCostMetrics;	// NumBytes,CommTime,IdleTime,EstCommCost,EstSynchCost,CompTime,OverlapTime
// Instead of printing out each Critter for each iteration individually, I will save them for each iteration, print out the iteration, and then clear before next iteration
extern std::map<std::string,std::tuple<double,double,double,double,double,double,double,double,double,double>> saveCritterInfo;
extern std::map<std::string,std::vector<std::string>> AlgCritters;
extern _critter MPI_Barrier_critter, 
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
         MPI_Comm_split_critter, 
         MPI_Sendrecv_replace_critter; 
void compute_all_crit(MPI_Comm cm, int nbr_pe, int nbr_pe2);
void compute_all_avg(MPI_Comm cm);
}
void print(size_t NumData, double* Data);
void start();
void stop();
}

#define MPI_Init(argc, argv)\
  do {\
     PMPI_Init(argc,argv);\
     critter::internal::flag = 0;\
     critter::internal::FileName="";\
     critter::internal::StreamName="";\
     if (std::getenv("CRITTER_VIZ") != NULL){\
       critter::internal::flag = 1;\
       critter::internal::FileName = std::move(std::string(*argv[*argc-1]));\
       critter::internal::StreamName = critter::internal::FileName + ".txt";\
       critter::internal::StreamTrackName = critter::internal::FileName + "track.txt";\
     }\
     critter::internal::IsFirstIter = true;\
     critter::internal::NeedNewLine = false;\
     int rank;\
     MPI_Comm_rank(MPI_COMM_WORLD,&rank);\
     if (rank == 0){\
       critter::internal::IsWorldRoot = true;\
     } else {critter::internal::IsWorldRoot=false;}\
     if (critter::internal::flag == 1){\
       if (rank==0){\
         critter::internal::Stream.open(critter::internal::StreamName.c_str());\
         critter::internal::StreamTrack.open(critter::internal::StreamTrackName.c_str());\
       }\
     } else{\
     }\
  } while (0)

#define MPI_Init_thread(argc, argv, required, provided)\
  do{\
     PMPI_Init_thread(argc,argv,required,provided);\
     critter::internal::flag = 0;\
     critter::internal::FileName="";\
     critter::internal::StreamName="";\
     if (std::getenv("CRITTER_VIZ") != NULL){\
       critter::internal::flag = 1;\
       critter::internal::FileName = std::move(std::string(*argv[*argc-1]));\
       critter::internal::StreamName = critter::internal::FileName + ".txt";\
       critter::internal::StreamTrackName = critter::internal::FileName + "track.txt";\
     }\
     critter::internal::IsFirstIter = true;\
     critter::internal::NeedNewLine = false;\
     int rank;\
     MPI_Comm_rank(MPI_COMM_WORLD,&rank);\
     if (rank == 0){\
       critter::internal::IsWorldRoot = true;\
     } else {critter::internal::IsWorldRoot=false;}\
     if (critter::internal::flag == 1){\
       if (rank==0){\
         critter::internal::Stream.open(critter::internal::StreamName.c_str());\
         critter::internal::StreamTrack.open(critter::internal::StreamTrackName.c_str());\
       }\
     } else{\
     }\
   } while (0)

#define MPI_Finalize()\
  do {\
    if (critter::internal::IsWorldRoot){\
      if (critter::internal::flag == 1){\
        critter::internal::Stream.close();\
      }\
    }\
    PMPI_Finalize();\
    } while (0)

#define MPI_Barrier(cm)\
  do {\
    if (critter::internal::track){\
      critter::internal::MPI_Barrier_critter.start(0, MPI_CHAR, cm);\
      PMPI_Barrier(cm);\
      critter::internal::MPI_Barrier_critter.stop();}\
    else{\
      PMPI_Barrier(cm);\
    }\
  } while (0)

#define MPI_Bcast(buf, nelem, t, root, cm)\
  do {\
    if (critter::internal::track){\
      critter::internal::MPI_Bcast_critter.start(nelem, t, cm);\
      PMPI_Bcast(buf, nelem, t, root, cm);\
      critter::internal::MPI_Bcast_critter.stop();}\
    else{\
      PMPI_Bcast(buf, nelem, t, root, cm);\
    }\
  } while (0)

#define MPI_Allreduce(sbuf, rbuf, nelem, t, op, cm)\
  do {\
    if (critter::internal::track){\
      critter::internal::MPI_Allreduce_critter.start(nelem, t, cm);\
      PMPI_Allreduce(sbuf, rbuf, nelem, t, op, cm);\
      critter::internal::MPI_Allreduce_critter.stop();}\
    else{\
      PMPI_Allreduce(sbuf, rbuf, nelem, t, op, cm);\
    }\
  } while (0)

#define MPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm)\
  do {\
    if (critter::internal::track){\
      critter::internal::MPI_Reduce_critter.start(nelem, t, cm);\
      PMPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm);\
      critter::internal::MPI_Reduce_critter.stop();}\
    else{\
      PMPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm);\
    }\
  } while (0)

#define MPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
  do {\
    if (critter::internal::track){\
      assert(rt==st); critter::internal::MPI_Scatter_critter.start(std::max((int64_t)scount,(int64_t)rcount), st, cm);\
      PMPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
      critter::internal::MPI_Scatter_critter.stop();}\
    else{\
      PMPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
    }\
  } while (0)

#define MPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      int pSize; MPI_Comm_size(cm, &pSize);\
      int64_t recvBufferSize = std::max((int64_t)scount,(int64_t)rcount) * pSize;\
      critter::internal::MPI_Gather_critter.start(recvBufferSize, st, cm);\
      PMPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
      critter::internal::MPI_Gather_critter.stop();}\
    else{\
      PMPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
    }\
  } while (0)

#define MPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      int pSize; MPI_Comm_size(cm, &pSize);\
      int64_t recvBufferSize = std::max((int64_t)scount,(int64_t)rcount) * pSize;\
      critter::internal::MPI_Allgather_critter.start(recvBufferSize, st, cm);\
      PMPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm);\
      critter::internal::MPI_Allgather_critter.stop();}\
    else{\
      PMPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm);\
    }\
  } while (0)

#define MPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm)\
  do {\
    if (critter::internal::track){\
      int64_t tot_recv=0;\
      int p; MPI_Comm_size(cm, &p);\
      for (int i=0; i<p; i++){ tot_recv += rcounts[i]; }\
      critter::internal::MPI_Reduce_scatter_critter.start(tot_recv, t, cm);\
      PMPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm);\
      critter::internal::MPI_Reduce_scatter_critter.stop();}\
    else{\
      PMPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm);\
    }\
  } while (0)

#define MPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm)\
  do {\
    if (critter::internal::track){\
      assert(rt==st); critter::internal::MPI_Alltoall_critter.start(std::max((int64_t)scount,(int64_t)rcount), st, cm);\
      PMPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm);\
      critter::internal::MPI_Alltoall_critter.stop();}\
    else{\
      PMPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm);\
    }\
  } while (0)

#define MPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      int64_t tot_recv=0;\
      int p; MPI_Comm_size(cm, &p);\
      for (int i=0; i<p; i++){ tot_recv += rcounts[i]; }\
      critter::internal::MPI_Allgatherv_critter.start(std::max((int64_t)scount,tot_recv), st, cm);\
      PMPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm);\
      critter::internal::MPI_Allgatherv_critter.stop();}\
    else{\
      PMPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm);\
    }\
  } while (0)


#define MPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      int64_t tot_recv=0;\
      int r, p; MPI_Comm_rank(cm, &r); MPI_Comm_size(cm, &p);\
      if (r == root) for (int i=0; i<p; i++){ tot_recv += ((int*)rcounts)[i]; }\
      critter::internal::MPI_Gatherv_critter.start(std::max((int64_t)scount,tot_recv), st, cm);\
      PMPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm);\
      critter::internal::MPI_Gatherv_critter.stop();}\
    else{\
      PMPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm);\
    }\
  } while (0)

#define MPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      int64_t tot_send=0;\
      int r, p; MPI_Comm_rank(cm, &r); MPI_Comm_size(cm, &p);\
      if (r == root) for (int i=0; i<p; i++){ tot_send += ((int*)scounts)[i]; } \
      critter::internal::MPI_Scatterv_critter.start(std::max(tot_send,(int64_t)rcount), st, cm);\
      PMPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm);\
      critter::internal::MPI_Scatterv_critter.stop();}\
    else{\
      PMPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm);\
    }\
  } while (0)

#define MPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      int64_t tot_send=0, tot_recv=0;\
      int p; MPI_Comm_size(cm, &p);\
      for (int i=0; i<p; i++){ tot_send += scounts[i]; tot_recv += rcounts[i]; }\
      critter::internal::MPI_Alltoallv_critter.start(std::max(tot_send,tot_recv), st, cm);\
      PMPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm);\
      critter::internal::MPI_Alltoallv_critter.stop();}\
    else{\
      PMPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm);\
    }\
  } while (0)


#define MPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status)\
  do {\
    if (critter::internal::track){\
      assert(st == rt);\
      critter::internal::MPI_Sendrecv_critter.start(std::max(scnt,rcnt), st, cm, dest, src);\
      PMPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status);\
      critter::internal::MPI_Sendrecv_critter.stop();}\
    else{\
      PMPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status);\
    }\
  } while (0)

#define MPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status)\
    do {\
    if (critter::internal::track){\
      critter::internal::MPI_Sendrecv_replace_critter.start(scnt, st, cm, dest, src);\
      PMPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status);\
      critter::internal::MPI_Sendrecv_replace_critter.stop();}\
    else{\
      PMPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status);\
    }\
  } while (0)

/*
#define MPI_Comm_split(comm1, arg2, arg3, cm)\
  do {\
    if (critter::internal::track){\
      critter::internal::MPI_Comm_split_critter.start(0, MPI_CHAR, *cm);\
      PMPI_Comm_split(comm1, arg2, arg3, cm);\
      critter::internal::MPI_Comm_split_critter.stop();}\
    else{\
      PMPI_Comm_split(comm1, arg2, arg3, cm);\
    }\
  } while (0)
*/

#define MPI_Send(buf, nelem, t, dest, tag, cm)\
  do {\
    if (critter::internal::track){\
      critter::internal::MPI_Send_critter.start(nelem, t, cm, dest);\
      PMPI_Send(buf, nelem, t, dest, tag, cm);\
      critter::internal::MPI_Send_critter.stop();}\
    else{\
      PMPI_Send(buf, nelem, t, dest, tag, cm);\
    }\
  } while (0)

#define MPI_Recv(buf, nelem, t, src, tag, cm, status)\
  do {\
    if (critter::internal::track){\
      critter::internal::MPI_Recv_critter.start(nelem, t, cm, src);\
      PMPI_Recv(buf, nelem, t, src, tag, cm, status);\
      critter::internal::MPI_Recv_critter.stop();}\
    else{\
      PMPI_Recv(buf, nelem, t, src, tag, cm, status);\
    }\
  } while (0)

#define MPI_Irecv(buf, nelem, t, src, tag, cm, req)\
  do {\
    if (critter::internal::track){\
      critter::internal::MPI_Irecv_critter.start(nelem, t, cm, src, -1, 1);\
      PMPI_Irecv(buf, nelem, t, src, tag, cm, req);\
      critter::internal::critter_req[*req] = &critter::internal::MPI_Irecv_critter;}\
    else{\
      PMPI_Irecv(buf, nelem, t, src, tag, cm, req);\
    }\
  } while (0)

#define MPI_Isend(buf, nelem, t, dest, tag, cm, req)\
  do {\
    if (critter::internal::track){\
      critter::internal::MPI_Isend_critter.start(nelem, t, cm, dest, -1, 1);\
      PMPI_Isend(buf, nelem, t, dest, tag, cm, req);\
      critter::internal::critter_req[*req] = &critter::internal::MPI_Isend_critter;}\
    else{\
      PMPI_Isend(buf, nelem, t, dest, tag, cm, req);\
    }\
  } while (0)

#define MPI_Wait(req, stat)\
  do {\
    if (critter::internal::track){\
      std::map<MPI_Request,critter::internal::_critter*>::iterator it = critter::internal::critter_req.find(*req);\
      if (it == critter::internal::critter_req.end()) *(int*)NULL = 1;\
      assert(it != critter::internal::critter_req.end());\
      PMPI_Wait(req, stat);\
      it->second->stop(); critter::internal::critter_req.erase(it);}\
    else{\
      PMPI_Wait(req, stat);\
    }\
  } while (0)

#define MPI_Waitany(cnt, reqs, indx, stat)\
  do {\
    if (critter::internal::track){\
      PMPI_Waitany(cnt, reqs, indx, stat);\
      std::map<MPI_Request,critter::internal::_critter*>::iterator it = critter::internal::critter_req.find((reqs)[*(indx)]);\
      if (it != critter::internal::critter_req.end()) { it->second->stop(); critter::internal::critter_req.erase(it);}}\
    else{\
      PMPI_Waitany(cnt, reqs, indx, stat);\
    }\
  } while (0)

#define MPI_Waitall(cnt, reqs, stats)\
  do {\
    if (critter::internal::track){\
      int __indx; MPI_Status __stat; for (int i=0; i<cnt; i++){ MPI_Waitany(cnt, reqs, &__indx, &__stat); if ((MPI_Status*)stats != (MPI_Status*)MPI_STATUSES_IGNORE) ((MPI_Status*)stats)[__indx] = __stat;}}\
    else{\
      PMPI_Waitall(cnt, reqs, stats)\
    }\
  } while (0)
#endif
