
#ifndef __CRITTER_H__
#define __CRITTER_H__

#include "mpi.h"
#include <fstream>
#include <iostream>
#include <utility>
#include "iomanip"
#include <vector>
#include <stack>
#include <array>
#include <stdint.h>
#include <functional>
#include <tuple>
#include <map>
#include <unordered_map>
#include <bitset>
#include <cmath>
#include <assert.h>

namespace critter{

// User functions
void start(size_t mode = 1);
void stop(size_t mode = 1, size_t factor = 1);

// User variables
// Note: `critical_path_breakdown_size` must equal `critical_path_breakdown.count()`. This will not be checked at compile time.
constexpr size_t critical_path_breakdown_size  = 1;
constexpr std::bitset<8> critical_path_breakdown(0b10000000); 	// RunTime,CompTime,DataMvtTime,SynchTime,CommTime,EstSynchCost,EstCommCost,NumBytes
constexpr int internal_tag                     = 1669220;	// arbitrary
constexpr bool p2p_blocking_comm_protocol      = true;		// 'false' for blocking p2p getting blockign protocol, 'true' for synchronous protocol
constexpr size_t max_timer_name_length = 50;			// max length of a symbol defining a timer
constexpr size_t max_num_symbols       = 50;			// max number of symbols to be tracked

// *****************************************************************************************************************************
namespace internal{

constexpr auto list_size 				= 32;		// numbers of tracked MPI routines
constexpr auto num_critical_path_measures 		= 8;		// NumBytes, EstCommCost, EstSynchCost, CommTime, SynchTime, DataMvtTime, CompTime,     RunTime
constexpr auto num_volume_measures 			= 9;		// NumBytes, EstCommCost, EstSynchCost, IdleTime, CommTime, SynchTime, DataMvtTime, CompTime,     RunTime
constexpr auto num_tracker_critical_path_measures 	= 6;		// Numbytes, EstCommCost, EstSynchCost, CommTime, SynchTime, DataMvtTime
constexpr auto num_tracker_volume_measures 		= 7;		// Numbytes, EstCommCost, EstSynchCost, IdleTime, CommTime, SynchTime, DataMvtTime,
constexpr auto num_ftimer_measures                      = 2;		// ExclusiveTime/Cost, InclusiveTime/Cost (NumCalls separate so as to avoid replication)

void update_critical_path(double* data);
void propagate_critical_path_synch(MPI_Comm cm, int partner);
void propagate_critical_path_blocking(MPI_Comm cm, int partner, bool is_sender);
void propagate_critical_path_nonblocking(double* data, MPI_Request internal_request, MPI_Comm cm, int partner, bool is_sender);
void compute_volume(MPI_Comm cm);

class tracker{
  public: 
    /* \brief name of collective */
    std::string name;
    /* \brief integer tag of collective */
    int tag;

    /* \brief local number of bytes max(sent,recv'ed) */
    double* my_bytes;
    /* \brief local duration of synchronization time */
    double* my_synch_time;
    /* \brief local duration of data mvt time */
    double* my_datamvt_time;
    /* \brief local duration of communication time */
    double* my_comm_time;
    /* \brief local duration of idle time */
    double* my_bar_time;
    /* \brief local comm cost in #messages */
    double* my_msg;
    /* \brief local comm cost in #words */
    double* my_wrd;

    /* \brief number of bytes max(sent,recv'ed) along the critical path among all processes */
    double* critical_path_bytes;
    /* \brief duration of synchronization time along the critical path among all processes */
    double* critical_path_synch_time;
    /* \brief number of time spent moving data about the network among all processes */
    double* critical_path_datamvt_time;
    /* \brief duration of communication time along the critical path among all processes */
    double* critical_path_comm_time;
    /* \brief comm cost in #messages along the critical path among all processes */
    double* critical_path_msg;
    /* \brief comm cost in #words along the critical path among all processes */
    double* critical_path_wrd;

    /* \brief duration of computation time for each call made locally,
     *   used to save the local computation time between calls to ::start and ::stop variants */
    double save_comp_time;
    volatile double save_time;

    /* \brief function for cost model of collective, takes (msg_size_in_bytes, number_processors) and returns (latency_cost, bandwidth_cost) */
    std::function< std::pair<double,double>(int64_t,int) > cost_func;

    /* \brief time when start() was last called, set to -1.0 initially and after stop() */
    volatile double last_start_time;
    /* \brief time when start() was last called, set to -1.0 initially and after stop() */
    volatile double last_synch_time;
    /* \brief save barrier time across start_synch */
    volatile double last_barrier_time;
    /* \brief cm with which start() was last called */
    MPI_Comm last_cm;
    /* \brief partner with which start() was last called */
    int last_partner;
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
    void start_synch(volatile double curTime, int64_t nelem, MPI_Datatype t=MPI_CHAR, MPI_Comm cm=MPI_COMM_WORLD, int partner=-1, int partner2=-1);
    void start_synch();

    /**
     * \brief starts timer for MPI call with nelem elements of type t over communicator cm, performs barrier over cm
     * \param[in] name symbol name of MPI routine
     * \param[in] cm MPI_Communicator on which MPI routine is called
     * \param[in] nbe_pe neighbor processor (only used for p2p routines)
     * \param[in] nbe_pe2 second neighbor processor (only used for p2p routines)
     * \param[in] is_async whether the call is asynchronous (used only for p2p routines)
     */
    void start_block(volatile double curTime, int64_t nelem, MPI_Datatype t=MPI_CHAR, MPI_Comm cm=MPI_COMM_WORLD, int partner=-1);
    void start_block();

    /**
     * \brief starts timer for MPI call with nelem elements of type t over communicator cm, performs barrier over cm
     * \param[in] name symbol name of MPI routine
     * \param[in] cm MPI_Communicator on which MPI routine is called
     * \param[in] nbe_pe neighbor processor (only used for p2p routines)
     * \param[in] nbe_pe2 second neighbor processor (only used for p2p routines)
     * \param[in] is_async whether the call is asynchronous (used only for p2p routines)
     */
    void start_nonblock(volatile double curTime, MPI_Request* request, int64_t nelem=1, MPI_Datatype t=MPI_CHAR, MPI_Comm cm=MPI_COMM_WORLD, bool is_sender=false, int partner=-1);
    void start_nonblock();

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
    void set_header();
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

// One instance for each unique symbol
class ftimer{
  public:
    ftimer() {}
    ftimer(std::string name_);
    void stop();
    void start();
    bool operator<(const ftimer& w) const ;

    std::string name;
    std::stack<double> start_timer;
    std::stack<std::unordered_map<std::string,std::array<double,num_critical_path_measures>>> exclusive_contributions;
    std::stack<std::array<double,num_critical_path_measures>> exclusive_measure;
    double* cp_numcalls; double* pp_numcalls; double* vol_numcalls;
    std::array<double*,num_critical_path_measures> cp_incl_measure;
    std::array<double*,num_critical_path_measures> cp_excl_measure;
    std::array<double*,num_critical_path_measures> pp_incl_measure;
    std::array<double*,num_critical_path_measures> pp_excl_measure;
    std::array<double*,num_critical_path_measures> vol_incl_measure;
    std::array<double*,num_critical_path_measures> vol_excl_measure;
    bool has_been_processed;
};


extern std::string stream_name,file_name;
extern bool flag,is_first_iter,is_world_root,need_new_line;
extern size_t mode;
extern std::ofstream stream;

extern double computation_timer;
extern std::map<MPI_Request,std::pair<MPI_Request,bool>> internal_comm_info;
extern std::map<MPI_Request,std::pair<MPI_Comm,int>> internal_comm_comm;
extern std::map<MPI_Request,double*> internal_comm_message;
extern std::map<MPI_Request,std::pair<double,double>> internal_comm_data;
extern std::map<MPI_Request,tracker*> internal_comm_track;
extern bool decisions[critical_path_breakdown_size];
constexpr auto critical_path_costs_size = num_critical_path_measures+num_tracker_critical_path_measures*critical_path_breakdown_size*list_size;
constexpr auto volume_costs_size = num_volume_measures+num_tracker_volume_measures*list_size;
extern std::array<double,critical_path_costs_size> critical_path_costs;
extern std::array<double,volume_costs_size> volume_costs;
extern std::array<double,num_critical_path_measures> max_per_process_costs;
extern std::map<std::string,std::vector<double>> save_info;
extern double new_cs[critical_path_costs_size];
extern double scratch_pad;
extern std::vector<char> synch_pad_send;
extern std::vector<char> synch_pad_recv;
extern std::array<char,max_timer_name_length*max_num_symbols> symbol_pad;
extern std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_local_cp;
extern std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_global_cp;
extern std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_pp;
extern std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_vol;
extern std::unordered_map<std::string,ftimer> symbol_timers;
extern std::stack<std::string> symbol_stack;
extern std::array<std::string,max_num_symbols> symbol_order;
extern std::array<std::string,num_critical_path_measures> critical_path_measure_names;
extern double_int timer_cp_info_sender[num_critical_path_measures];
extern double_int timer_cp_info_receiver[num_critical_path_measures];
}
}

// **********************************************************************************************************

#define TAU_START(ARG) do {\
  if (critter::internal::mode==2){\
  if (critter::internal::symbol_timers.find(#ARG) == critter::internal::symbol_timers.end()){\
    critter::internal::symbol_timers[#ARG] = critter::internal::ftimer(#ARG);\
    critter::internal::symbol_order[critter::internal::symbol_timers.size()-1] = #ARG;\
    critter::internal::symbol_timers[#ARG].start();\
  }\
  else{\
    critter::internal::symbol_timers[#ARG].start();\
  }}}while (0);

#define TAU_FSTART(ARG) do {\
  if (critter::internal::mode==2){\
  if (critter::internal::symbol_timers.find(#ARG) == critter::internal::symbol_timers.end()){\
    critter::internal::symbol_timers[#ARG] = critter::internal::ftimer(#ARG);\
    critter::internal::symbol_order[critter::internal::symbol_timers.size()-1] = #ARG;\
    critter::internal::symbol_timers[#ARG].start();\
  }\
  else{\
    critter::internal::symbol_timers[#ARG].start();\
  }}}while (0);

#define TAU_STOP(ARG) do {\
  if (critter::internal::mode==2){\
  if (critter::internal::symbol_timers.find(#ARG) == critter::internal::symbol_timers.end()){\
    assert(0);\
  }\
  else{\
    critter::internal::symbol_timers[#ARG].stop();\
  }}}while (0);

#define TAU_FSTOP(ARG) do {\
  if (critter::internal::mode==2){\
  if (critter::internal::symbol_timers.find(#ARG) == critter::internal::symbol_timers.end()){\
    assert(0);\
  }\
  else{\
    critter::internal::symbol_timers[#ARG].stop();\
  }}}while (0);

// **********************************************************************************************************

#define MPI_Init(argc, argv)\
  do {\
     assert(critter::critical_path_breakdown_size == critter::critical_path_breakdown.count());\
     PMPI_Init(argc,argv);\
     critter::internal::mode=0;\
     critter::internal::flag = 0;\
     critter::internal::file_name="";\
     critter::internal::stream_name="";\
     if (std::getenv("CRITTER_VIZ_FILE") != NULL){\
       critter::internal::flag = 1;\
       critter::internal::file_name = std::getenv("CRITTER_VIZ_FILE");\
       critter::internal::stream_name = critter::internal::file_name + ".txt";\
     }\
     critter::internal::is_first_iter = true;\
     critter::internal::need_new_line = false;\
     int _critter_rank,_critter_size;\
     MPI_Comm_rank(MPI_COMM_WORLD,&_critter_rank);\
     MPI_Comm_size(MPI_COMM_WORLD,&_critter_size);\
     critter::internal::synch_pad_send.resize(_critter_size);\
     critter::internal::synch_pad_recv.resize(_critter_size);\
     if (_critter_rank == 0){\
       critter::internal::is_world_root = true;\
     } else {critter::internal::is_world_root=false;}\
     if (critter::internal::flag == 1){\
       if (_critter_rank==0){\
         critter::internal::stream.open(critter::internal::stream_name.c_str());\
       }\
     } else{\
     }\
  } while (0)

#define MPI_Init_thread(argc, argv, required, provided)\
  do{\
     assert(critter::critical_path_breakdown_size == critter::critical_path_breakdown.count());\
     PMPI_Init_thread(argc,argv,required,provided);\
     critter::internal::mode=0;\
     critter::internal::flag = 0;\
     critter::internal::file_name="";\
     critter::internal::stream_name="";\
     if (std::getenv("CRITTER_VIZ_FILE") != NULL){\
       critter::internal::flag = 1;\
       critter::internal::file_name = std::getenv("CRITTER_VIZ_FILE");\
       critter::internal::stream_name = critter::internal::file_name + ".txt";\
     }\
     critter::internal::is_first_iter = true;\
     critter::internal::need_new_line = false;\
     int _critter_rank,_critter_size;\
     MPI_Comm_rank(MPI_COMM_WORLD,&_critter_rank);\
     MPI_Comm_size(MPI_COMM_WORLD,&_critter_size);\
     critter::internal::synch_pad_send.resize(_critter_size);\
     critter::internal::synch_pad_recv.resize(_critter_size);\
     if (_critter_rank == 0){\
       critter::internal::is_world_root = true;\
     } else {critter::internal::is_world_root=false;}\
     if (critter::internal::flag == 1){\
       if (_critter_rank==0){\
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
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      critter::internal::_MPI_Barrier.start_synch(_critter_curTime_, 0, MPI_CHAR, cm);\
      PMPI_Barrier(cm);\
      critter::internal::_MPI_Barrier.start_synch();\
      PMPI_Barrier(cm);\
      critter::internal::_MPI_Barrier.stop_synch();}\
    else{\
      PMPI_Barrier(cm);\
    }\
  } while (0)

#define MPI_Bcast(buf, nelem, t, root, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      critter::internal::_MPI_Bcast.start_synch(_critter_curTime_, nelem, t, cm);\
      PMPI_Bcast(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, root, cm);\
      critter::internal::_MPI_Bcast.start_synch();\
      PMPI_Bcast(buf, nelem, t, root, cm);\
      critter::internal::_MPI_Bcast.stop_synch();\
    }\
    else{\
      PMPI_Bcast(buf, nelem, t, root, cm);\
    }\
  } while (0)

#define MPI_Allreduce(sbuf, rbuf, nelem, t, op, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      critter::internal::_MPI_Allreduce.start_synch(_critter_curTime_, nelem, t, cm);\
      PMPI_Allreduce(MPI_IN_PLACE, &critter::internal::synch_pad_send[0], 1, MPI_CHAR, MPI_MAX, cm);\
      critter::internal::_MPI_Allreduce.start_synch();\
      PMPI_Allreduce(sbuf, rbuf, nelem, t, op, cm);\
      critter::internal::_MPI_Allreduce.stop_synch();}\
    else{\
      PMPI_Allreduce(sbuf, rbuf, nelem, t, op, cm);\
    }\
  } while (0)

#define MPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      critter::internal::_MPI_Reduce.start_synch(_critter_curTime_, nelem, t, cm);\
      PMPI_Reduce(&critter::internal::synch_pad_send[0], &critter::internal::synch_pad_recv[0], 1, MPI_CHAR, MPI_MAX, root, cm);\
      critter::internal::_MPI_Reduce.start_synch();\
      PMPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm);\
      critter::internal::_MPI_Reduce.stop_synch();} \
    else{\
      PMPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm);\
    }\
  } while (0)

#define MPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      int64_t _critter_sendBufferSize = std::max((int64_t)scount,(int64_t)rcount) * _critter_np;\
      critter::internal::_MPI_Scatter.start_synch(_critter_curTime_, _critter_sendBufferSize, st, cm);\
      PMPI_Scatter(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, &critter::internal::synch_pad_recv[0], 1, MPI_CHAR, root, cm);\
      critter::internal::_MPI_Scatter.start_synch();\
      PMPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
      critter::internal::_MPI_Scatter.stop_synch();}\
    else{\
      PMPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
    }\
  } while (0)

#define MPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      int64_t _critter_recvBufferSize = std::max((int64_t)scount,(int64_t)rcount) * _critter_np;\
      critter::internal::_MPI_Gather.start_synch(_critter_curTime_, _critter_recvBufferSize, st, cm);\
      PMPI_Gather(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, &critter::internal::synch_pad_recv[0], 1, MPI_CHAR, root, cm);\
      critter::internal::_MPI_Gather.start_synch();\
      PMPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
      critter::internal::_MPI_Gather.stop_synch();}\
    else{\
      PMPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
    }\
  } while (0)

#define MPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      int64_t _critter_recvBufferSize = std::max((int64_t)scount,(int64_t)rcount) * _critter_np;\
      critter::internal::_MPI_Allgather.start_synch(_critter_curTime_, _critter_recvBufferSize, st, cm);\
      PMPI_Allgather(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, &critter::internal::synch_pad_recv[0], 1, MPI_CHAR, cm);\
      critter::internal::_MPI_Allgather.start_synch();\
      PMPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm);\
      critter::internal::_MPI_Allgather.stop_synch();}\
    else{\
      PMPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm);\
    }\
  } while (0)

#define MPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      int64_t _critter_tot_recv=0;\
      int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      std::vector<int> _critter_rcounts(p,1);\
      for (int _critter_i=0; _critter_i<_critter_np; _critter_i++){ _critter_tot_recv += rcounts[_critter_i]; }\
      critter::internal::_MPI_Reduce_scatter.start_synch(_critter_curTime_, _critter_tot_recv, t, cm);\
      PMPI_Reduce_scatter(&critter::internal::synch_pad_send[0], &critter::internal::synch_pad_recv[0], &_critter_rcounts[0], MPI_CHAR, MPI_MAX, cm);\
      critter::internal::_MPI_Reduce_scatter.start_synch();\
      PMPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm);\
      critter::internal::_MPI_Reduce_scatter.stop_synch();}\
    else{\
      PMPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm);\
    }\
  } while (0)

#define MPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      int64_t _critter_recvBufferSize = std::max((int64_t)scount,(int64_t)rcount) * _critter_np;\
      critter::internal::_MPI_Alltoall.start_synch(_critter_curTime_,_critter_recvBufferSize, st, cm);\
      PMPI_Alltoall(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, &critter::internal::synch_pad_recv[0], 1, MPI_CHAR, cm);\
      critter::internal::_MPI_Alltoall.start_synch();\
      PMPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm);\
      critter::internal::_MPI_Alltoall.stop_synch();}\
    else{\
      PMPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm);\
    }\
  } while (0)

#define MPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int64_t _critter_tot_recv=0; int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      std::vector<int> _critter_rcounts(_critter_np,1); std::vector<int> _critter_rdisp(_critter_np,0);\
      for (int _critter_i=1; _critter_i<_critter_np; _critter_i++) _critter_rdisp[_critter_i]=_critter_rdisp[_critter_i-1]+1;\
      for (int _critter_i=0; _critter_i<_critter_np; _critter_i++){ _critter_tot_recv += rcounts[_critter_i]; }\
      critter::internal::_MPI_Allgatherv.start_synch(_critter_curTime_, std::max((int64_t)scount,_critter_tot_recv), st, cm);\
      PMPI_Allgatherv(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, &critter::internal::synch_pad_recv[0], &_critter_rcounts[0], &_critter_rdisp[0], MPI_CHAR, cm);\
      critter::internal::_MPI_Allgatherv.start_synch();\
      PMPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm);\
      critter::internal::_MPI_Allgatherv.stop_synch();}\
    else{\
      PMPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm);\
    }\
  } while (0)


#define MPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int64_t _critter_tot_recv=0;\
      int _critter_rank, _critter_np; MPI_Comm_rank(cm, &_critter_rank); MPI_Comm_size(cm, &_critter_np);\
      std::vector<int> _critter_rcounts(_critter_np,0);\
      if (_critter_rank == root) for (int _critter_i=0; _critter_i<_critter_np; _critter_i++){ _critter_tot_recv += ((int*)rcounts)[_critter_i]; }\
      critter::internal::_MPI_Gatherv.start_synch(_critter_curTime_, std::max((int64_t)scount,_critter_tot_recv), st, cm);\
      PMPI_Gatherv(nullptr, 0, st, nullptr, &_critter_rcounts[0], &_critter_rcounts[0], rt, root, cm);\
      critter::internal::_MPI_Gatherv.start_synch();\
      PMPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm);\
      critter::internal::_MPI_Gatherv.stop_synch();}\
    else{\
      PMPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm);\
    }\
  } while (0)

#define MPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int64_t _critter_tot_send=0;\
      int _critter_r, _critter_np; MPI_Comm_rank(cm, &_critter_rank); MPI_Comm_size(cm, &_critter_np);\
      std::vector<int> _critter_scounts(_critter_np,1); std::vector<int> _critter_sdisp(_critter_np,0); for (_critter_i=1; _critter_i<_critter_np; _critter_i++) _critter_sdisp[i]=_critter_sdisp[_critter_i-1]+1;\
      if (_critter_rank == root) for (int _critter_i=0; _critter_i<_critter_np; _critter_i++){ _critter_tot_send += ((int*)scounts)[_critter_i]; } \
      critter::internal::_MPI_Scatterv.start_synch(_critter_curTime_, std::max(_critter_tot_send,(int64_t)rcount), st, cm);\
      PMPI_Scatterv(&critter::internal::synch_pad_send[0], &_critter_scounts[0], &_critter_sdisp[0], MPI_CHAR, &critter::internal::synch_pad_recv[0], 1, MPI_CHAR, root, cm);\
      critter::internal::_MPI_Scatterv.start_synch();\
      PMPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm);\
      critter::internal::_MPI_Scatterv.stop_synch();}\
    else{\
      PMPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm);\
    }\
  } while (0)

#define MPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int64_t _critter_tot_send=0, _critter_tot_recv=0;\
      int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      for (int _critter_i=0; _critter_i<_critter_np; _critter_i++){ _critter_tot_send += scounts[_critter_i]; _critter_tot_recv += rcounts[_critter_i]; }\
      std::vector<int> _critter_counts(_critter_np,1); std::vector<int> _critter_disp(_critter_np,0);\
      for (int _critter_i=1; _critter_i<_critter_np; _critter_i++) _critter_disp[_critter_i]=_critter_disp[_critter_i-1]+1;\
      critter::internal::_MPI_Alltoallv.start_synch(_critter_curTime_, std::max(_critter_tot_send,_critter_tot_recv), st, cm);\
      PMPI_Alltoallv(&critter::internal::synch_pad_send[0], &_critter_counts[0], &_critter_disp[0], MPI_CHAR, &critter::internal::synch_pad_recv[0], &_critter_counts[0], &_critter_disp[0], MPI_CHAR, cm);\
      critter::internal::_MPI_Scatterv.start_synch();\
      PMPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm);\
      critter::internal::_MPI_Alltoallv.stop_synch();}\
    else{\
      PMPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm);\
    }\
  } while (0)


#define MPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(st == rt); assert(stag != critter::internal_tag); assert(rtag != critter::internal_tag);\
      critter::internal::_MPI_Sendrecv.start_synch(_critter_curTime_, std::max(scnt,rcnt), st, cm, dest, src);\
      PMPI_Sendrecv(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, dest, stag, &critter::internal::synch_pad_recv[0], 1, MPI_CHAR, src, rtag, cm, status);\
      critter::internal::_MPI_Sendrecv.start_synch();\
      PMPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status);\
      critter::internal::_MPI_Sendrecv.stop_synch();}\
    else{\
      PMPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status);\
    }\
  } while (0)

#define MPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(stag != critter::internal_tag); assert(rtag != critter::internal_tag);\
      critter::internal::_MPI_Sendrecv_replace.start_synch(_critter_curTime_, scnt, st, cm, dest, src);\
      PMPI_Sendrecv_replace(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, dest, stag, src, rtag, cm, status);\
      critter::internal::_MPI_Sendrecv_replace.start_synch();\
      PMPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status);\
      critter::internal::_MPI_Sendrecv_replace.stop_synch();}\
    else{\
      PMPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status);\
    }\
  } while (0)

#define MPI_Ssend(buf, nelem, t, dest, tag, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(tag != critter::internal_tag);\
      critter::internal::_MPI_Ssend.start_synch(_critter_curTime_, nelem, t, cm, dest);\
      PMPI_Ssend(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, dest, tag, cm);\
      critter::internal::_MPI_Ssend.start_synch();\
      PMPI_Ssend(buf, nelem, t, dest, tag, cm);\
      critter::internal::_MPI_Ssend.stop_synch();\
    }\
    else{\
      PMPI_Ssend(buf, nelem, t, dest, tag, cm);\
    }\
  } while (0)

#define MPI_Send(buf, nelem, t, dest, tag, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(tag != critter::internal_tag);\
      if (!critter::p2p_blocking_comm_protocol) {critter::internal::_MPI_Send.start_block(_critter_curTime_, nelem, t, cm, dest);}\
      else {critter::internal::_MPI_Send.start_synch(_critter_curTime_, nelem, t, cm, dest); }\
      PMPI_Send(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, dest, tag, cm);\
      if (!critter::p2p_blocking_comm_protocol) {critter::internal::_MPI_Send.start_block();}\
      else {critter::internal::_MPI_Send.start_synch(); }\
      PMPI_Send(buf, nelem, t, dest, tag, cm);\
      if (!critter::p2p_blocking_comm_protocol) {critter::internal::_MPI_Send.stop_block(true);}\
      else {critter::internal::_MPI_Send.stop_synch(); }\
    }\
    else{\
      PMPI_Send(buf, nelem, t, dest, tag, cm);\
    }\
  } while (0)

#define MPI_Recv(buf, nelem, t, src, tag, cm, status)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(tag != critter::internal_tag);\
      if (!critter::p2p_blocking_comm_protocol) {critter::internal::_MPI_Recv.start_block(_critter_curTime_, nelem, t, cm, src);}\
      else {critter::internal::_MPI_Recv.start_synch(_critter_curTime_, nelem, t, cm, src); }\
      PMPI_Recv(&critter::internal::synch_pad_recv[0], 1, MPI_CHAR, src, tag, cm, status);\
      if (!critter::p2p_blocking_comm_protocol) {critter::internal::_MPI_Recv.start_block();}\
      else {critter::internal::_MPI_Recv.start_synch(); }\
      PMPI_Recv(buf, nelem, t, src, tag, cm, status);\
      if (!critter::p2p_blocking_comm_protocol) {critter::internal::_MPI_Recv.stop_block(false);}\
      else {critter::internal::_MPI_Recv.stop_synch(); }\
    }\
    else{\
      PMPI_Recv(buf, nelem, t, src, tag, cm, status);\
    }\
  } while (0)

#define MPI_Isend(buf, nelem, t, dest, tag, cm, req)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(tag != critter::internal_tag);\
      PMPI_Isend(buf, nelem, t, dest, tag, cm, req);\
      critter::internal::_MPI_Isend.start_nonblock(_critter_curTime_, req, nelem, t, cm, true, dest);\
      critter::internal::computation_timer = MPI_Wtime();\
    }\
    else{\
      PMPI_Isend(buf, nelem, t, dest, tag, cm, req);\
    }\
  } while (0)

#define MPI_Irecv(buf, nelem, t, src, tag, cm, req)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(tag != critter::internal_tag);\
      PMPI_Irecv(buf, nelem, t, src, tag, cm, req);\
      critter::internal::_MPI_Irecv.start_nonblock(_critter_curTime_, req, nelem, t, cm, false, src);\
      critter::internal::computation_timer = MPI_Wtime();\
    }\
    else{\
      PMPI_Irecv(buf, nelem, t, src, tag, cm, req);\
    }\
  } while (0)

#define MPI_Ibcast(buf, nelem, t, root, cm, req)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      PMPI_Ibcast(buf, nelem, t, root, cm, req);\
      critter::internal::_MPI_Ibcast.start_nonblock(_critter_curTime_, req, nelem, t, cm);\
      critter::internal::computation_timer = MPI_Wtime();\
    }\
    else{\
      PMPI_Ibcast(buf, nelem, t, root, cm, req);\
    }\
  } while (0)

#define MPI_Iallreduce(sbuf, rbuf, nelem, t, op, cm, req)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      PMPI_Iallreduce(sbuf, rbuf, nelem, t, op, cm, req);\
      critter::internal::_MPI_Iallreduce.start_nonblock(_critter_curTime_, req, nelem, t, cm);\
      critter::internal::computation_timer = MPI_Wtime();\
    }\
    else{\
      PMPI_Iallreduce(sbuf, rbuf, nelem, t, op, cm, req);\
    }\
  } while (0)

#define MPI_Ireduce(sbuf, rbuf, nelem, t, op, root, cm, req)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      PMPI_Ireduce(buf, nelem, t, root, cm, req);\
      critter::internal::_MPI_Iallreduce.start_nonblock(_critter_curTime_, req, nelem, t, cm);\
      critter::internal::computation_timer = MPI_Wtime();\
    }\
    else{\
      PMPI_Ireduce(sbuf, rbuf, nelem, t, op, root, cm, req);\
    }\
  } while (0)

#define MPI_Igather(sbuf, scount, st, rbuf, rcount, rt, root, cm, req)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      int64_t _critter_recvBufferSize = std::max((int64_t)scount,(int64_t)rcount) * _critter_np;\
      PMPI_Igather(sbuf, scount, st, rbuf, rcount, rt, root, cm, req);\
      critter::internal::_MPI_Igather.start_nonblock(_critter_curTime_, req, _critter_recvBufferSize, st, cm);\
    }\
    else{\
      PMPI_Igather(sbuf, scount, st, rbuf, rcount, rt, root, cm, req);\
    }\
  } while (0)

#define MPI_Igatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm, req)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int64_t _critter_tot_recv=0; int _critter_rank, _critter_npp; MPI_Comm_rank(cm, &_critter_rank); MPI_Comm_size(cm, &_critter_np);\
      if (_critter_rank == root) for (int _critter_i=0; _critter_i<_critter_np; _critter_i++){ _critter_tot_recv += ((int*)rcounts)[_critter_i]; }\
      PMPI_Igatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm, req);\
      critter::internal::_MPI_Igatherv.start_nonblock(_critter_curTime_, req, std::max((int64_t)scount,_critter_tot_recv), st, cm);\
    }\
    else{\
      PMPI_Igatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm, req);\
    }\
  } while (0)

#define MPI_Iallgather(sbuf, scount, st, rbuf, rcount, rt, cm, req)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int _critter_np; MPI_Comm_size(cm, &_critter_np); int64_t _critter_recvBufferSize = std::max((int64_t)scount,(int64_t)rcount) * _critter_np;\
      PMPI_Iallgather(sbuf, scount, st, rbuf, rcount, rt, cm, req);\
      critter::internal::_MPI_Iallgather.start_nonblock(_critter_curTime_, req, _critter_recvBufferSize, st, cm);\
    }\
    else{\
      PMPI_Iallgather(sbuf, scount, st, rbuf, rcount, rt, cm, req);\
    }\
  } while (0)

#define MPI_Iallgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm, req)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int64_t _critter_tot_recv=0; int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      for (int _critter_i=0; _critter_i<_critter_np; _critter_i++){ _critter_tot_recv += rcounts[_critter_i]; }\
      PMPI_Iallgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm, req);\
      critter::internal::_MPI_Iallgatherv.start_nonblock(_critter_curTime_, req, std::max((int64_t)scount,_critter_tot_recv), st, cm);\
    }\
    else{\
      PMPI_Iallgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm, req);\
    }\
  } while (0)

#define MPI_Iscatter(sbuf, scount, st, rbuf, rcount, rt, root, cm, req)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      int64_t _critter_sendBufferSize = std::max((int64_t)scount,(int64_t)rcount) * _critter_np;\
      PMPI_Iscatter(sbuf, scount, st, rbuf, rcount, rt, root, cm, req);\
      critter::internal::_MPI_Iscatter.start_nonblock(_critter_curTime_, req, _critter_sendBufferSize, st, cm);\
    }\
    else{\
      PMPI_Iscatter(sbuf, scount, st, rbuf, rcount, rt, root, cm, req);\
    }\
  } while (0)

#define MPI_Iscatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm, req)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int64_t _critter_tot_send=0;\
      int _critter_rank, _critter_np; MPI_Comm_rank(cm, &_critter_rank); MPI_Comm_size(cm, &_critter_np);\
      if (_critter_rank == root) for (int _critter_i=0; _critter_i<_critter_np; _critter_i++){ _critter_tot_send += ((int*)scounts)[_critter_i]; } \
      PMPI_Iscatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm, req);\
      critter::internal::_MPI_Iscatterv.start_nonblock(_critter_curTime_, req, std::max(_critter_tot_send,(int64_t)rcount), st, cm);\
    }\
    else{\
      PMPI_Iscatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm, req);\
    }\
  } while (0)

#define MPI_Ireduce_scatter(sbuf, rbuf, rcounts, t, op, cm, req)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      int64_t _critter_tot_recv=0;\
      int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      for (int _critter_i=0; _critter_i<_critter_np; _critter_i++){ _critter_tot_recv += rcounts[_critter_i]; }\
      PMPI_Ireduce_scatter(sbuf, rbuf, rcounts, t, op, cm, req);\
      critter::internal::_MPI_Ireduce_scatter.start_nonblock(_critter_curTime_, req, _critter_tot_recv, t, cm);\
    }\
    else{\
      PMPI_Ireduce_scatter(sbuf, rbuf, rcounts, t, op, cm, req);\
    }\
  } while (0)

#define MPI_Ialltoall(sbuf, scount, st, rbuf, rcount, rt, cm, req)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      PMPI_Ialltoall(sbuf, scount, st, rbuf, rcount, rt, cm, req);\
      critter::internal::_MPI_Ialltoall.start_nonblock(_critter_curTime_, req, std::max((int64_t)scount,(int64_t)rcount)*_critter_np, st, cm);\
    }\
    else{\
      PMPI_Ialltoall(sbuf, scount, st, rbuf, rcount, rt, cm, req);\
    }\
  } while (0)

#define MPI_Ialltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm, req)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int64_t _critter_tot_send=0, _critter_tot_recv=0;\
      int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      for (int _critter_i=0; _critter_i<_criter_np; _critter_i++){ _critter_tot_send += scounts[_critter_i]; _critter_tot_recv += rcounts[_critter_i]; }\
      PMPI_Ialltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm, req);\
      critter::internal::_MPI_Ialltoallv.start_nonblock(_critter_curTime_, req, std::max(_critter_tot_send,_critter_tot_recv), st, cm);\
    }\
    else{\
      PMPI_Ialltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm, req);\
    }\
  } while (0)

#define MPI_Wait(req, stat)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime = MPI_Wtime(); double _critter_save_comp_time = _critter_curTime - critter::internal::computation_timer;\
      auto _critter_comm_track_it = critter::internal::internal_comm_track.find(*req);\
      assert(_critter_comm_track_it != critter::internal::internal_comm_track.end());\
      auto _critter_comm_info_it = critter::internal::internal_comm_info.find(*req);\
      MPI_Request _critter_save_request = _critter_comm_info_it->first;\
      volatile double _critter_last_start_time = MPI_Wtime();\
      PMPI_Wait(req, stat);\
      _critter_curTime = MPI_Wtime(); double _critter_save_comm_time = _critter_curTime - _critter_last_start_time;\
      _critter_comm_track_it->second->stop_nonblock(&_critter_save_request, _critter_save_comp_time, _critter_save_comm_time);\
    }\
    else{\
      PMPI_Wait(req, stat);\
    }\
  } while (0)

#define MPI_Waitany(cnt, reqs, indx, stat)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime = MPI_Wtime(); double _critter_save_comp_time = _critter_curTime - critter::internal::computation_timer;\
      std::vector<MPI_Request> _critter_pt(cnt); for (int _critter_i=0;_critter_i<cnt;_critter_i++){_critter_pt[_critter_i]=(reqs)[_critter_i];}\
      volatile double _critter_last_start_time = MPI_Wtime();\
      PMPI_Waitany(cnt, reqs, indx, stat);\
      _critter_curTime = MPI_Wtime(); double _critter_save_comm_time = _critter_curTime - _critter_last_start_time;\
      MPI_Request _critter_request = _critter_pt[*indx];\
      auto _critter_comm_track_it = critter::internal::internal_comm_track.find(_critter_request);\
      assert(_critter_comm_track_it != critter::internal::internal_comm_track.end());\
      _critter_comm_track_it->second->stop_nonblock(&_critter_request, _critter_save_comp_time, _critter_save_comm_time);\
    }\
    else{\
      PMPI_Waitany(cnt, reqs, indx, stat);\
    }\
  } while (0)

#define MPI_Waitsome(incnt, reqs, outcnt, indices, stats)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime = MPI_Wtime(); double _critter_save_comp_time = _critter_curTime - critter::internal::computation_timer;\
      std::vector<MPI_Request> _critter_pt(incnt); for (int _critter_i=0;_critter_i<incnt;_critter_i++){_critter_pt[_critter_i]=(reqs)[_critter_i];}\
      volatile double _critter_last_start_time = MPI_Wtime();\
      PMPI_Waitsome(incnt, reqs, outcnt, indices, stats);\
      _critter_curTime = MPI_Wtime(); double _critter_save_comm_time = _critter_curTime - _critter_last_start_time;\
      for (int _critter_i=0; _critter_i<*outcnt; _critter_i++){\
        MPI_Request _critter_request = _critter_pt[indices[_critter_i]];\
        auto _critter_comm_track_it = critter::internal::internal_comm_track.find(_critter_request);\
        assert(_critter_comm_track_it != critter::internal::internal_comm_track.end());\
        _critter_comm_track_it->second->stop_nonblock(&_critter_request, _critter_save_comp_time, _critter_save_comm_time);\
      }\
    }\
    else{\
      PMPI_Waitsome(incnt, reqs, outcnt, indices, stats);\
    }\
  } while (0)

#define MPI_Waitall(cnt, reqs, stats)\
  do {\
    if (critter::internal::mode>=1){\
      int _critter_indx; MPI_Status _critter_stat; for (int _critter_i=0; _critter_i<cnt; _critter_i++){ MPI_Waitany(cnt, reqs, &_critter_indx, &_critter_stat); if ((MPI_Status*)stats != (MPI_Status*)MPI_STATUSES_IGNORE) ((MPI_Status*)stats)[_critter_indx] = _critter_stat;}\
    }\
    else{\
      PMPI_Waitall(cnt, reqs, stats);\
    }\
  } while (0)

#endif
