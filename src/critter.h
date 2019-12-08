
#ifndef __CRITTER_H__
#define __CRITTER_H__

#include "mpi.h"
#include <fstream>
#include <iostream>
#include <utility>
#include "iomanip"
#include <vector>
#include <stdint.h>
#include <functional>
#include <tuple>
#include <map>
#include <cmath>
#include <assert.h>

namespace critter{

// User functions
void print(size_t num_data, double* data);
void start();
void stop();

namespace internal{

constexpr auto internal_tag 				= 1669220;	// arbitrary
constexpr auto list_size 				= 32;		// numbers of tracked MPI routines
constexpr auto num_critical_path_measures 		= 6;		// NumBytes,CommTime,EstCommCost,EstSynchCost,CompTime,RunTime
constexpr auto num_volume_measures 			= 7;		// NumBytes,CommTime,IdleTime,EstCommCost,EstSynchCost,CompTime,RunTime
constexpr auto num_tracker_critical_path_measures 	= 4;		// numbytes, commtime, estcomm, estsynch
constexpr auto num_tracker_volume_measures 		= 5;		// numbytes, commtime, barriertime, estcomm, estsynch

void update_critical_path(double* data);
void propagate_critical_path(MPI_Comm cm, int nbr_pe, int nbr_pe2);
void compute_volume(MPI_Comm cm);

class tracker{
  public: 
    /* \brief name of collective */
    std::string name;
    /* \brief integer tag of collective */
    int tag;

    /* \brief local number of bytes max(sent,recv'ed) */
    double* my_bytes;
    /* \brief local duration of communication time */
    double* my_comm_time;
    /* \brief local duration of idle time */
    double* my_bar_time;
    /* \brief local comm cost in #messages */
    double* my_msg;
    /* \brief local comm cost in #words */
    double* my_wrd;

    /* \brief number of bytes max(sent,recv'ed) along the critical_pathical path among all processes */
    double* critical_path_bytes;
    /* \brief duration of communication time along the critical_pathical path among all processes */
    double* critical_path_comm_time;
    /* \brief comm cost in #messages along the critical_pathical path among all processes */
    double* critical_path_msg;
    /* \brief comm cost in #words along the critical_pathical path among all processes */
    double* critical_path_wrd;

    /* \brief duration of computation time for each call made locally,
     *   used to save the local computation time between calls to ::start and ::stop variants */
    double save_comp_time;

    /* \brief function for cost model of collective, takes (msg_size_in_bytes, number_processors) and returns (latency_cost, bandwidth_cost) */
    std::function< std::pair<double,double>(int64_t,int) > cost_func;

    /* \brief time when start() was last called, set to -1.0 initially and after stop() */
    double last_start_time;
    /* \brief save barrier time across start_synch */
    double last_barrier_time;
    /* \brief cm with which start() was last called */
    MPI_Comm last_cm;
    /* \brief nbr_pe with which start() was last called */
    int last_nbr_pe;
    /* \brief nbr_pe2 with which start() was last called */
    int last_nbr_pe2;
    /* \brief helper for nonblocking comm */
    int64_t last_nbytes;
    /* \brief helper for nonblocking comm */
    int last_p;

    /**
     * \brief timer constructor, initializes vars
     * \param[in] name symbol name of MPI routine
     * \param[in] function for cost model of collective, takes (msg_size_in_bytes, number_processors) and returns (latency_cost, bandwidth_cost) 
     */
    tracker(std::string name, int tag,
            std::function< std::pair<double,double>(int64_t,int) > 
              cost_func = [](int64_t n, int p){ 
                return std::pair<double,double>(1.,n); 
              });
 
    /**
     * \brief timer copy constructor, copies name
     * \param[in] t other timer
     */
    tracker(tracker const& t);
   
    /**
     * \brief timer destructor, frees name
     */ 
    ~tracker();

    /**
     * \brief starts timer for MPI call with nelem elements of type t over communicator cm, performs barrier over cm
     * \param[in] name symbol name of MPI routine
     * \param[in] cm MPI_Communicator on which MPI routine is called
     * \param[in] nbe_pe neighbor processor (only used for p2p routines)
     * \param[in] nbe_pe2 second neighbor processor (only used for p2p routines)
     * \param[in] is_async whether the call is asynchronous (used only for p2p routines)
     */
    void start_synch(int64_t nelem=1, MPI_Datatype t=MPI_CHAR, MPI_Comm cm=MPI_COMM_WORLD, int nbr_pe=-1, int nbr_pe2=-1);

    /**
     * \brief starts timer for MPI call with nelem elements of type t over communicator cm, performs barrier over cm
     * \param[in] name symbol name of MPI routine
     * \param[in] cm MPI_Communicator on which MPI routine is called
     * \param[in] nbe_pe neighbor processor (only used for p2p routines)
     * \param[in] nbe_pe2 second neighbor processor (only used for p2p routines)
     * \param[in] is_async whether the call is asynchronous (used only for p2p routines)
     */
    void start_block(int64_t nelem=1, MPI_Datatype t=MPI_CHAR, MPI_Comm cm=MPI_COMM_WORLD, int nbr_pe=-1, int nbr_pe2=-1);

    /**
     * \brief starts timer for MPI call with nelem elements of type t over communicator cm, performs barrier over cm
     * \param[in] name symbol name of MPI routine
     * \param[in] cm MPI_Communicator on which MPI routine is called
     * \param[in] nbe_pe neighbor processor (only used for p2p routines)
     * \param[in] nbe_pe2 second neighbor processor (only used for p2p routines)
     * \param[in] is_async whether the call is asynchronous (used only for p2p routines)
     */
    void start_nonblock(MPI_Request* request, int64_t nelem=1, MPI_Datatype t=MPI_CHAR, MPI_Comm cm=MPI_COMM_WORLD, bool is_sender=false, int nbr_pe=-1, int nbr_pe2=-1);

    /**
     * \brief completes interception of synchronous communication protocol
     */
    void stop_synch();

    /**
     * \brief completes interception of blocking communication protocol
     * \param[in] is_sender whether or not the calling process is the sending process
     */
    void stop_block(bool is_sender);

    /**
     * \brief completes interception of nonblocking communication protocol
     */
    void stop_nonblock(MPI_Request* request, double comp_time, double comm_time);

    /**
     */
    void set_critical_path_costs(size_t idx);
    /**
     */
    void set_volume_costs();

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
     * \brief sets data members to point into global arrays
     */
    void set_cost_pointers();
};

extern tracker
         _MPI_Barrier,
         _MPI_Bcast,
         _MPI_Reduce,
         _MPI_Allreduce,
         _MPI_Gather,
         _MPI_Gatherv,
         _MPI_Allgather,
         _MPI_Allgatherv,
         _MPI_Scatter,
         _MPI_Scatterv,
         _MPI_Reduce_scatter,
         _MPI_Alltoall,
         _MPI_Alltoallv,
         _MPI_Send,
         _MPI_Recv,
         _MPI_Isend,
         _MPI_Irecv,
         _MPI_Sendrecv,
         _MPI_Sendrecv_replace,
         _MPI_Ibcast,
         _MPI_Iallreduce,
         _MPI_Ireduce,
         _MPI_Igather,
         _MPI_Igatherv,
         _MPI_Iallgather,
         _MPI_Iallgatherv,
         _MPI_Iscatter,
         _MPI_Iscatterv,
         _MPI_Ireduce_scatter,
         _MPI_Ialltoall,
         _MPI_Ialltoallv,
         _MPI_Ssend;

extern tracker* list[list_size];
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

extern std::string stream_name,file_name;
extern bool track,flag,is_first_iter,is_world_root,need_new_line;
extern std::ofstream stream;

extern double computation_timer;
extern std::map<MPI_Request,std::pair<MPI_Request,bool>> internal_comm_info;
extern std::map<MPI_Request,double*> internal_comm_message;
extern std::map<MPI_Request,std::pair<double,double>> internal_comm_data;
extern std::map<MPI_Request,tracker*> internal_comm_track;
constexpr auto critical_path_costs_size = num_critical_path_measures+num_tracker_critical_path_measures*num_critical_path_measures*list_size;
constexpr auto volume_costs_size = num_volume_measures+num_tracker_volume_measures*list_size;
extern std::array<double,critical_path_costs_size> critical_path_costs;
extern std::array<double,volume_costs_size> volume_costs;
extern std::map<std::string,std::vector<double>> save_info;
extern double new_cs[critical_path_costs_size];

}
}

#define MPI_Init(argc, argv)\
  do {\
     PMPI_Init(argc,argv);\
     critter::internal::track=0;\
     critter::internal::flag = 0;\
     critter::internal::file_name="";\
     critter::internal::stream_name="";\
     if (std::getenv("CRITTER_VIZ") != NULL){\
       critter::internal::flag = 1;\
       critter::internal::file_name = std::getenv("CRITTER_VIZ_FILE");\
       critter::internal::stream_name = critter::internal::file_name + ".txt";\
     }\
     critter::internal::is_first_iter = true;\
     critter::internal::need_new_line = false;\
     int rank;\
     MPI_Comm_rank(MPI_COMM_WORLD,&rank);\
     if (rank == 0){\
       critter::internal::is_world_root = true;\
     } else {critter::internal::is_world_root=false;}\
     if (critter::internal::flag == 1){\
       if (rank==0){\
         critter::internal::stream.open(critter::internal::stream_name.c_str());\
       }\
     } else{\
     }\
  } while (0)

#define MPI_Init_thread(argc, argv, required, provided)\
  do{\
     PMPI_Init_thread(argc,argv,required,provided);\
     critter::internal::track=0;\
     critter::internal::flag = 0;\
     critter::internal::file_name="";\
     critter::internal::stream_name="";\
     if (std::getenv("CRITTER_VIZ") != NULL){\
       critter::internal::flag = 1;\
       critter::internal::file_name = std::getenv("CRITTER_VIZ_FILE");\
       critter::internal::stream_name = critter::internal::file_name + ".txt";\
     }\
     critter::internal::is_first_iter = true;\
     critter::internal::need_new_line = false;\
     int rank;\
     MPI_Comm_rank(MPI_COMM_WORLD,&rank);\
     if (rank == 0){\
       critter::internal::is_world_root = true;\
     } else {critter::internal::is_world_root=false;}\
     if (critter::internal::flag == 1){\
       if (rank==0){\
         critter::internal::stream.open(critter::internal::stream_name.c_str());\
       }\
     } else{\
     }\
   } while (0)

#define MPI_Finalize()\
  do {\
    if (critter::internal::is_world_root){\
      if (critter::internal::flag == 1){\
        critter::internal::stream.close();\
      }\
    }\
    PMPI_Finalize();\
    } while (0)

#define MPI_Barrier(cm)\
  do {\
    if (critter::internal::track){\
      critter::internal::_MPI_Barrier.start_synch(0, MPI_CHAR, cm);\
      PMPI_Barrier(cm);\
      critter::internal::_MPI_Barrier.stop_synch();}\
    else{\
      PMPI_Barrier(cm);\
    }\
  } while (0)

#define MPI_Bcast(buf, nelem, t, root, cm)\
  do {\
    if (critter::internal::track){\
      critter::internal::_MPI_Bcast.start_synch(nelem, t, cm);\
      PMPI_Bcast(buf, nelem, t, root, cm);\
      critter::internal::_MPI_Bcast.stop_synch();\
    }\
    else{\
      PMPI_Bcast(buf, nelem, t, root, cm);\
    }\
  } while (0)

#define MPI_Allreduce(sbuf, rbuf, nelem, t, op, cm)\
  do {\
    if (critter::internal::track){\
      critter::internal::_MPI_Allreduce.start_synch(nelem, t, cm);\
      PMPI_Allreduce(sbuf, rbuf, nelem, t, op, cm);\
      critter::internal::_MPI_Allreduce.stop_synch();}\
    else{\
      PMPI_Allreduce(sbuf, rbuf, nelem, t, op, cm);\
    }\
  } while (0)

#define MPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm)\
  do {\
    if (critter::internal::track){\
      critter::internal::_MPI_Reduce.start_synch(nelem, t, cm);\
      PMPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm);\
      critter::internal::_MPI_Reduce.stop_synch();}\
    else{\
      PMPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm);\
    }\
  } while (0)

#define MPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      critter::internal::_MPI_Scatter.start_synch(std::max((int64_t)scount,(int64_t)rcount), st, cm);\
      PMPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
      critter::internal::_MPI_Scatter.stop_synch();}\
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
      critter::internal::_MPI_Gather.start_synch(recvBufferSize, st, cm);\
      PMPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
      critter::internal::_MPI_Gather.stop_synch();}\
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
      critter::internal::_MPI_Allgather.start_synch(recvBufferSize, st, cm);\
      PMPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm);\
      critter::internal::_MPI_Allgather.stop_synch();}\
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
      critter::internal::_MPI_Reduce_scatter.start_synch(tot_recv, t, cm);\
      PMPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm);\
      critter::internal::_MPI_Reduce_scatter.stop_synch();}\
    else{\
      PMPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm);\
    }\
  } while (0)

#define MPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      critter::internal::_MPI_Alltoall.start_synch(std::max((int64_t)scount,(int64_t)rcount), st, cm);\
      PMPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm);\
      critter::internal::_MPI_Alltoall.stop_synch();}\
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
      critter::internal::_MPI_Allgatherv.start_synch(std::max((int64_t)scount,tot_recv), st, cm);\
      PMPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm);\
      critter::internal::_MPI_Allgatherv.stop_synch();}\
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
      critter::internal::_MPI_Gatherv.start_synch(std::max((int64_t)scount,tot_recv), st, cm);\
      PMPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm);\
      critter::internal::_MPI_Gatherv.stop_synch();}\
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
      critter::internal::_MPI_Scatterv.start_synch(std::max(tot_send,(int64_t)rcount), st, cm);\
      PMPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm);\
      critter::internal::_MPI_Scatterv.stop_synch();}\
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
      critter::internal::_MPI_Alltoallv.start_synch(std::max(tot_send,tot_recv), st, cm);\
      PMPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm);\
      critter::internal::_MPI_Alltoallv.stop_synch();}\
    else{\
      PMPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm);\
    }\
  } while (0)


#define MPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status)\
  do {\
    if (critter::internal::track){\
      assert(st == rt); assert(stag != critter::internal::internal_tag); assert(rtag != critter::internal::internal_tag);\
      critter::internal::_MPI_Sendrecv.start_synch(std::max(scnt,rcnt), st, cm, dest, src);\
      PMPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status);\
      critter::internal::_MPI_Sendrecv.stop_synch();}\
    else{\
      PMPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status);\
    }\
  } while (0)

#define MPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status)\
    do {\
    if (critter::internal::track){\
      assert(stag != critter::internal::internal_tag); assert(rtag != critter::internal::internal_tag);\
      critter::internal::_MPI_Sendrecv_replace.start_synch(scnt, st, cm, dest, src);\
      PMPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status);\
      critter::internal::_MPI_Sendrecv_replace.stop_synch();}\
    else{\
      PMPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status);\
    }\
  } while (0)

#define MPI_Ssend(buf, nelem, t, dest, tag, cm)\
  do {\
    if (critter::internal::track){\
      assert(tag != critter::internal::internal_tag);\
      critter::internal::_MPI_Ssend.start_synch(nelem, t, cm, dest);\
      PMPI_Ssend(buf, nelem, t, dest, tag, cm);\
      critter::internal::_MPI_Ssend.stop_synch();\
    }\
    else{\
      PMPI_Ssend(buf, nelem, t, dest, tag, cm);\
    }\
  } while (0)

#define MPI_Send(buf, nelem, t, dest, tag, cm)\
  do {\
    if (critter::internal::track){\
      assert(tag != critter::internal::internal_tag);\
      critter::internal::_MPI_Send.start_block(nelem, t, cm, dest);\
      PMPI_Send(buf, nelem, t, dest, tag, cm);\
      critter::internal::_MPI_Send.stop_block(true);\
    }\
    else{\
      PMPI_Send(buf, nelem, t, dest, tag, cm);\
    }\
  } while (0)

#define MPI_Isend(buf, nelem, t, dest, tag, cm, req)\
  do {\
    if (critter::internal::track){\
      assert(tag != critter::internal::internal_tag);\
      PMPI_Isend(buf, nelem, t, dest, tag, cm, req);\
      critter::internal::_MPI_Isend.start_nonblock(req, nelem, t, cm, true, dest);\
      critter::internal::computation_timer = MPI_Wtime();\
    }\
    else{\
      PMPI_Isend(buf, nelem, t, dest, tag, cm, req);\
    }\
  } while (0)

#define MPI_Recv(buf, nelem, t, src, tag, cm, status)\
  do {\
    if (critter::internal::track){\
      assert(tag != critter::internal::internal_tag);\
      critter::internal::_MPI_Recv.start_block(nelem, t, cm, src);\
      PMPI_Recv(buf, nelem, t, src, tag, cm, status);\
      critter::internal::_MPI_Recv.stop_block(false);\
    }\
    else{\
      PMPI_Recv(buf, nelem, t, src, tag, cm, status);\
    }\
  } while (0)

#define MPI_Irecv(buf, nelem, t, src, tag, cm, req)\
  do {\
    if (critter::internal::track){\
      assert(tag != critter::internal::internal_tag);\
      PMPI_Irecv(buf, nelem, t, src, tag, cm, req);\
      critter::internal::_MPI_Irecv.start_nonblock(req, nelem, t, cm, false, src);\
      critter::internal::computation_timer = MPI_Wtime();\
    }\
    else{\
      PMPI_Irecv(buf, nelem, t, src, tag, cm, req);\
    }\
  } while (0)

#define MPI_Ibcast(buf, nelem, t, root, cm, req)\
  do {\
    if (critter::internal::track){\
      PMPI_Ibcast(buf, nelem, t, root, cm);\
      critter::internal::_MPI_Ibcast.start_nonblock(req, nelem, t, cm);\
      critter::internal::computation_timer = MPI_Wtime();\
    }\
    else{\
      PMPI_Ibcast(buf, nelem, t, root, cm, req);\
    }\
  } while (0)

#define MPI_Iallreduce(sbuf, rbuf, nelem, t, op, cm, req)\
  do {\
    if (critter::internal::track){\
      PMPI_Iallreduce(buf, nelem, t, root, cm, req);\
      critter::internal::_MPI_Iallreduce.start_nonblock(req, nelem, t, cm);\
      critter::internal::computation_timer = MPI_Wtime();\
    else{\
      PMPI_Iallreduce(sbuf, rbuf, nelem, t, op, cm, req);\
    }\
  } while (0)

#define MPI_Ireduce(sbuf, rbuf, nelem, t, op, root, cm, req)\
  do {\
    if (critter::internal::track){\
      PMPI_Ireduce(buf, nelem, t, root, cm, req);\
      critter::internal::_MPI_Iallreduce.start_nonblock(req, nelem, t, cm);\
      critter::internal::computation_timer = MPI_Wtime();\
    else{\
      PMPI_Ireduce(sbuf, rbuf, nelem, t, op, root, cm, req);\
    }\
  } while (0)

#define MPI_Igather(sbuf, scount, st, rbuf, rcount, rt, root, cm, req)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      int pSize; MPI_Comm_size(cm, &pSize);\
      int64_t recvBufferSize = std::max((int64_t)scount,(int64_t)rcount) * pSize;\
      PMPI_Igather(sbuf, scount, st, rbuf, rcount, rt, root, cm, req);\
      critter::internal::_MPI_Igather.start_nonblock(req, recvBufferSize, st, cm);\
    else{\
      PMPI_Igather(sbuf, scount, st, rbuf, rcount, rt, root, cm, req);\
    }\
  } while (0)

#define MPI_Igatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm, req)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      int64_t tot_recv=0;\
      int r, p; MPI_Comm_rank(cm, &r); MPI_Comm_size(cm, &p);\
      if (r == root) for (int i=0; i<p; i++){ tot_recv += ((int*)rcounts)[i]; }\
      PMPI_Igatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm, req);\
      critter::internal::_MPI_Igatherv.start_nonblock(req, std::max((int64_t)scount,tot_recv), st, cm);\
    else{\
      PMPI_Igatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm, req);\
    }\
  } while (0)

#define MPI_Iallgather(sbuf, scount, st, rbuf, rcount, rt, cm, req)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      int pSize; MPI_Comm_size(cm, &pSize);\
      int64_t recvBufferSize = std::max((int64_t)scount,(int64_t)rcount) * pSize;\
      PMPI_Iallgather(sbuf, scount, st, rbuf, rcount, rt, cm, req);\
      critter::internal::_MPI_Iallgather.start_nonblock(req, recvBufferSize, st, cm);\
    else{\
      PMPI_Iallgather(sbuf, scount, st, rbuf, rcount, rt, cm, req);\
    }\
  } while (0)

#define MPI_Iallgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm, req)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      int64_t tot_recv=0;\
      int p; MPI_Comm_size(cm, &p);\
      for (int i=0; i<p; i++){ tot_recv += rcounts[i]; }\
      PMPI_Iallgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm, req);\
      critter::internal::_MPI_Iallgatherv.start_nonblock(req, std::max((int64_t)scount,tot_recv), st, cm);\
    else{\
      PMPI_Iallgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm, req);\
    }\
  } while (0)

#define MPI_Iscatter(sbuf, scount, st, rbuf, rcount, rt, root, cm, req)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      PMPI_Iscatter(sbuf, scount, st, rbuf, rcount, rt, root, cm, req);\
      critter::internal::_MPI_Iscatter.start_nonblock(req, std::max((int64_t)scount,(int64_t)rcount), st, cm);\
    else{\
      PMPI_Iscatter(sbuf, scount, st, rbuf, rcount, rt, root, cm, req);\
    }\
  } while (0)

#define MPI_Iscatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm, req)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      int64_t tot_send=0;\
      int r, p; MPI_Comm_rank(cm, &r); MPI_Comm_size(cm, &p);\
      if (r == root) for (int i=0; i<p; i++){ tot_send += ((int*)scounts)[i]; } \
      PMPI_Iscatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm, req);\
      critter::internal::_MPI_Iscatterv.start_nonblock(req, std::max(tot_send,(int64_t)rcount), st, cm);\
    else{\
      PMPI_Iscatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm, req);\
    }\
  } while (0)

#define MPI_Ireduce_scatter(sbuf, rbuf, rcounts, t, op, cm, req)\
  do {\
    if (critter::internal::track){\
      int64_t tot_recv=0;\
      int p; MPI_Comm_size(cm, &p);\
      for (int i=0; i<p; i++){ tot_recv += rcounts[i]; }\
      PMPI_Ireduce_scatter(sbuf, rbuf, rcounts, t, op, cm, req);\
      critter::internal::_MPI_Ireduce_scatter.start_nonblock(req, tot_recv, t, cm);\
    else{\
      PMPI_Ireduce_scatter(sbuf, rbuf, rcounts, t, op, cm, req);\
    }\
  } while (0)

#define MPI_Ialltoall(sbuf, scount, st, rbuf, rcount, rt, cm, req)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      PMPI_Ialltoall(sbuf, scount, st, rbuf, rcount, rt, cm, req);\
      critter::internal::_MPI_Ialltoall.start_nonblock(req, std::max((int64_t)scount,(int64_t)rcount), st, cm);\
    else{\
      PMPI_Ialltoall(sbuf, scount, st, rbuf, rcount, rt, cm, req);\
    }\
  } while (0)

#define MPI_Ialltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm, req)\
  do {\
    if (critter::internal::track){\
      assert(rt==st);\
      int64_t tot_send=0, tot_recv=0;\
      int p; MPI_Comm_size(cm, &p);\
      for (int i=0; i<p; i++){ tot_send += scounts[i]; tot_recv += rcounts[i]; }\
      PMPI_Ialltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm, req);\
      critter::internal::_MPI_Ialltoallv.start_nonblock(req, std::max(tot_send,tot_recv), st, cm);\
    else{\
      PMPI_Ialltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm, req);\
    }\
  } while (0)

#define MPI_Wait(req, stat)\
  do {\
    if (critter::internal::track){\
      volatile double curTime = MPI_Wtime(); double save_comp_time = curTime - critter::internal::computation_timer;\
      auto comm_track_it = critter::internal::internal_comm_track.find(*req);\
      assert(comm_track_it != critter::internal::internal_comm_track.end());\
      auto comm_info_it = critter::internal::internal_comm_info.find(*req);\
      MPI_Request save_request = comm_info_it->first;\
      volatile double last_start_time = MPI_Wtime();\
      PMPI_Wait(req, stat);\
      curTime = MPI_Wtime(); double save_comm_time = curTime - last_start_time;\
      comm_track_it->second->stop_nonblock(&save_request, save_comp_time, save_comm_time);\
    }\
    else{\
      PMPI_Wait(req, stat);\
    }\
  } while (0)

#define MPI_Waitany(cnt, reqs, indx, stat)\
  do {\
    if (critter::internal::track){\
      volatile double curTime = MPI_Wtime(); double save_comp_time = curTime - critter::internal::computation_timer;\
      std::vector<MPI_Request> pt(cnt); for (int i=0;i<cnt;i++){pt[i]=(reqs)[i];}\
      volatile double last_start_time = MPI_Wtime();\
      PMPI_Waitany(cnt, reqs, indx, stat);\
      curTime = MPI_Wtime(); double save_comm_time = curTime - last_start_time;\
      MPI_Request request = pt[*indx];\
      auto comm_track_it = critter::internal::internal_comm_track.find(request);\
      assert(comm_track_it != critter::internal::internal_comm_track.end());\
      comm_track_it->second->stop_nonblock(&request, save_comp_time, save_comm_time);\
    }\
    else{\
      PMPI_Waitany(cnt, reqs, indx, stat);\
    }\
  } while (0)

#define MPI_Waitsome(incnt, reqs, outcnt, indices, stats)\
  do {\
    if (critter::internal::track){\
      volatile double curTime = MPI_Wtime(); double save_comp_time = curTime - critter::internal::computation_timer;\
      std::vector<MPI_Request> pt(cnt); for (int i=0;i<cnt;i++){pt[i]=(reqs)[i];}\
      volatile double last_start_time = MPI_Wtime();\
      PMPI_Waitsome(incnt, reqs, outcnt, indices, stats);\
      curTime = MPI_Wtime(); double save_comm_time = curTime - last_start_time;\
      for (int i=0; i<*outcnt; i++){\
        MPI_Request request = pt[indices[i]];\
        auto comm_track_it = critter::internal::internal_comm_track.find(request);\
        assert(comm_track_it != critter::internal::internal_comm_track.end());\
        comm_track_it->second->stop_nonblock(&request, save_comp_time, save_comm_time);\
      }\
    }\
    else{\
      PMPI_Waitany(cnt, reqs, indx, stat);\
    }\
  } while (0)

#define MPI_Waitall(cnt, reqs, stats)\
  do {\
    if (critter::internal::track){\
      int __indx; MPI_Status __stat; for (int i=0; i<cnt; i++){ MPI_Waitany(cnt, reqs, &__indx, &__stat); if ((MPI_Status*)stats != (MPI_Status*)MPI_STATUSES_IGNORE) ((MPI_Status*)stats)[__indx] = __stat;}\
    }\
    else{\
      PMPI_Waitall(cnt, reqs, stats);\
    }\
  } while (0)

#endif
