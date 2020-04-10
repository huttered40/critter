#ifndef __CRITTER_H__
#define __CRITTER_H__

#include <mpi.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <vector>
#include <stack>
#include <array>
#include <stdint.h>
#include <functional>
#include <map>
#include <unordered_map>
#include <bitset>
#include <cmath>
#include <assert.h>

namespace critter{

// *****************************************************************************************************************************************************************
// Forward declarations needed for 'p2p_type' (see 'User functions' listed below)
namespace internal{
class synchronous; class blocking;
};

// *****************************************************************************************************************************************************************
// User functions
void start(size_t mode = 1);
void stop(size_t mode = 1, size_t factor = 1);

// User variables
// Note: `cost_model_size` must equal `cost_models.count()`. This will not be checked at compile time.
constexpr size_t cost_model_size  = 2;					// must match number of bits set in cost_model (below)
constexpr std::bitset<2> cost_models(0b11);				// alpha-beta butterfly, BSP
// Note: `breakdown_size` must equal `breakdown.count()`. This will not be checked at compile time.
constexpr size_t breakdown_size  			= 5;			// must match number of bits set in breakdown (below)
constexpr std::bitset<5+2*cost_model_size> breakdown(0b110010110);  		// RunTime,CompTime,DataMvtTime,SynchTime,CommTime,ABSynchCost,BSPsynchCost,ABCommCost,BSPcommCost
constexpr int internal_tag                      	= 1669220;		// arbitrary
constexpr int internal_tag1                     	= 1669221;		// arbitrary
constexpr int internal_tag2                     	= 1669222;		// arbitrary
constexpr int internal_tag3                     	= 1669223;		// arbitrary
constexpr int internal_tag4                     	= 1669224;		// arbitrary
using p2p_type 						= internal::blocking;// p2p communication can be tracked as 'synchronous' or 'blocking'
constexpr size_t max_timer_name_length 			= 20;			// max length of a symbol defining a timer
constexpr size_t max_num_symbols       			= 15;			// max number of symbols to be tracked

// *****************************************************************************************************************************************************************
namespace internal{

constexpr auto list_size 				= 32;				// numbers of tracked MPI routines
constexpr auto num_critical_path_measures 		= 5+2*cost_model_size;		// CommCost*, SynchCost*,           CommTime, SynchTime, DataMvtTime, CompTime, RunTime
constexpr auto num_per_process_measures 		= 6+2*cost_model_size;		// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, DataMvtTime, CompTime, RunTime
constexpr auto num_volume_measures 			= 6+2*cost_model_size;		// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, DataMvtTime, CompTime, RunTime
constexpr auto num_tracker_critical_path_measures 	= 3+2*cost_model_size;		// CommCost*, SynchCost*,           CommTime, SynchTime, DataMvtTime
constexpr auto num_tracker_per_process_measures 	= 3+2*cost_model_size;		// CommCost*, SynchCost*,           CommTime, SynchTime, DataMvtTime,
constexpr auto num_tracker_volume_measures 		= 3+2*cost_model_size;		// CommCost*, SynchCost*,           CommTime, SynchTime, DataMvtTime,
constexpr auto num_ftimer_measures                      = 2;				// ExclusiveTime/Cost, InclusiveTime/Cost (NumCalls separate so as to avoid replication)
constexpr auto critical_path_costs_size = num_critical_path_measures+num_tracker_critical_path_measures*breakdown_size*list_size+breakdown_size;
constexpr auto per_process_costs_size = num_per_process_measures+num_tracker_per_process_measures*breakdown_size*list_size+2*breakdown_size;
constexpr auto volume_costs_size = num_volume_measures+num_tracker_volume_measures*list_size;
constexpr auto mode_1_width = 25;
constexpr auto mode_2_width = 15;

void update_critical_path(double* data);
void compute_volume(MPI_Comm cm);
void complete_propagation(MPI_Comm, bool is_sender, int partner1, int partner2);

/* \brief encapsulation of the state of a MPI routine */
class tracker{
  public: 
    /* \brief name of MPI routine */
    std::string name;
    /* \brief integer tag of MPI routine */
    int tag;
    /* \brief local duration of synchronization time */
    double* my_synch_time;
    /* \brief local duration of data mvt time */
    double* my_datamvt_time;
    /* \brief local duration of communication time */
    double* my_comm_time;
    /* \brief local comm cost in #messages */
    double* my_msg_count;
    /* \brief local comm cost in #words */
    double* my_wrd_count;
    /* \brief duration of synchronization time along a critical path */
    double* critical_path_synch_time;
    /* \brief number of time spent moving data about the network along a critical path */
    double* critical_path_datamvt_time;
    /* \brief duration of communication time along a critical path */
    double* critical_path_comm_time;
    /* \brief comm cost in #messages along a critical path */
    double* critical_path_msg_count;
    /* \brief comm cost in #words along a critical path */
    double* critical_path_wrd_count;
    /* \brief function for cost model of MPI routine in bsp cost model, takes (msg_size_in_bytes, number_processors) and returns (latency_cost, bandwidth_cost) */
    std::function< std::pair<double,double>(int64_t,int) > cost_func_bsp;
    /* \brief function for cost model of MPI routine in alpha-beta butterfly cost model, takes (msg_size_in_bytes, number_processors) and returns (latency_cost, bandwidth_cost) */
    std::function< std::pair<double,double>(int64_t,int) > cost_func_alphabeta_butterfly;
    /* \brief duration of computation time for each call made locally, used to save the local computation time between calls to ::start and ::stop variants */
    double save_comp_time;
    /* \brief variable to save time across start and stop tracking */
    volatile double save_time;
    /* \brief time when start() was last called, set to -1.0 initially and after stop() */
    volatile double last_start_time;
    /* \brief time when start() was last called, set to -1.0 initially and after stop() */
    volatile double last_synch_time;
    /* \brief save barrier time across start_synch */
    volatile double last_barrier_time;
    /* \brief cm with which start() was last called */
    MPI_Comm last_cm;
    /* \brief partner with which start() was last called */
    int last_partner1;
    /* \brief partner with which start() was last called */
    int last_partner2;
    /* \brief nbytes with which start() was last called */
    int64_t last_nbytes;
    /* \brief process count with which start() was last called */
    int last_p;
    /* \brief root bool with which start() was last called */
    bool last_is_root;
    /* \brief is_sender bool with which start() was last called */
    bool last_is_sender;
    /** \brief initialization of state called by construtors */
    void init();
    /** */
    void set_header();
    /** */
    void set_critical_path_costs(size_t idx);
    /** */
    void set_per_process_costs(size_t idx);
    /** */
    void set_volume_costs();
    /** \brief sets data members to point into global arrays */
    void set_cost_pointers();
};

class synchronous : public tracker{
public:
    /**
     * \brief constructor
     * \param[in] name symbol name of MPI routine
     * \param[in] tag integer id of MPI routine
     * \param[in] cost_func_bsp function for simple cost model of MPI routine
     * \param[in] cost_func_alphabeta_butterfly function for alpha-beta cost model of MPI routine assuming butterfly algorithms
     */
    synchronous(std::string name, int tag,
            std::function< std::pair<double,double>(int64_t,int)> 
              cost_func_bsp = [](int64_t n, int p){ return std::pair<double,double>(1.,n); },
            std::function< std::pair<double,double>(int64_t,int)> 
              cost_func_alphabeta_butterfly = [](int64_t n, int p){ return std::pair<double,double>(1.,n); }
               );
    /** \brief copy constructor */
    synchronous(synchronous const& t);
    /** \brief starts tracking of synchronous MPI communication of nelem elements of type t over communicator cm */
    void start(volatile double curTime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, bool is_root = false, bool is_sender = false, int partner=-1, int partner2=-1);
    /** \brief starts communication timer for corresponding synchronous MPI communication */
    void intermediate();
    /** \brief completes interception of synchronous communication protocol */
    void stop();
    /** \brief propagates path information */
    void propagate(MPI_Comm cm, bool is_sender=false, int partner1=-1, int partner2=-1);
};

class blocking : public tracker{
public:
    /**
     * \brief constructor
     * \param[in] name symbol name of MPI routine
     * \param[in] tag integer id of MPI routine
     * \param[in] cost_func_bsp function for simple cost model of MPI routine
     * \param[in] cost_func_alphabeta_butterfly function for alpha-beta cost model of MPI routine assuming butterfly algorithms
     */
    blocking(std::string name, int tag,
            std::function< std::pair<double,double>(int64_t,int)> 
              cost_func_bsp = [](int64_t n, int p){ return std::pair<double,double>(1.,n); },
            std::function< std::pair<double,double>(int64_t,int)> 
              cost_func_alphabeta_butterfly = [](int64_t n, int p){ return std::pair<double,double>(1.,n); }
            );
    /** \brief copy constructor */
    blocking(blocking const& t);
    /** \brief starts tracking of blocking MPI communication of nelem elements of type t over communicator cm */
    void start(volatile double curTime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, bool is_root, bool is_sender, int partner1, int partner2=-1);
    /** \brief starts communication timer for corresponding blocking MPI communication */
    void intermediate();
    /** \brief completes interception of blocking communication protocol */
    void stop();
    /** \brief propagates path information */
    void propagate(MPI_Comm cm, bool is_sender, int partner1, int partner2=-1);
};

class nonblocking : public tracker{
public:
    /**
     * \brief constructor
     * \param[in] name symbol name of MPI routine
     * \param[in] tag integer id of MPI routine
     * \param[in] cost_func_bsp function for simple cost model of MPI routine
     * \param[in] cost_func_alphabeta_butterfly function for alpha-beta cost model of MPI routine assuming butterfly algorithms
     */
    nonblocking(std::string name, int tag,
            std::function< std::pair<double,double>(int64_t,int)> 
              cost_func_bsp = [](int64_t n, int p){ return std::pair<double,double>(1.,n); },
            std::function< std::pair<double,double>(int64_t,int)> 
              cost_func_alphabeta_butterfly = [](int64_t n, int p){ return std::pair<double,double>(1.,n); }
               );
    /** \brief copy constructor */
    nonblocking(nonblocking const& t);
    /** \brief starts tracking of nonblocking MPI communication of nelem elements of type t over communicator cm */
    void start(volatile double curTime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, MPI_Request* request, bool is_sender=false, int partner=-1);
    /** \brief starts communication timer for corresponding nonblocking MPI communication */
    void intermediate();
    /** \brief completes interception of nonblocking communication protocol */
    void stop(MPI_Request* request, double comp_time, double comm_time);
    /** \brief propagates path information */
    void propagate(std::pair<int*,int*> envelope_int, std::array<double*,3> envelope_double, char* envelope_char, std::array<MPI_Request,6> internal_requests, MPI_Comm cm, bool is_sender, int partner);
};

extern synchronous
         _MPI_Barrier,
         _MPI_Bcast,
         _MPI_Reduce,
         _MPI_Allreduce,
         _MPI_Gather,
         _MPI_Allgather,
         _MPI_Scatter,
         _MPI_Reduce_scatter,
         _MPI_Alltoall,
         _MPI_Gatherv,
         _MPI_Allgatherv,
         _MPI_Scatterv,
         _MPI_Alltoallv,
         _MPI_Ssend;
extern p2p_type
         _MPI_Send,
         _MPI_Recv,
         _MPI_Sendrecv,
         _MPI_Sendrecv_replace;
extern nonblocking
         _MPI_Isend,
         _MPI_Irecv,
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
         _MPI_Ialltoallv;

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
    void stop(double save_time);
    void start(double save_time);
    bool operator<(const ftimer& w) const ;

    std::string name;
    std::stack<double> start_timer;
    std::array<double,num_critical_path_measures> cp_exclusive_contributions;
    std::array<double,num_per_process_measures> pp_exclusive_contributions;
    std::array<double,num_critical_path_measures> cp_exclusive_measure;
    std::array<double,num_per_process_measures> pp_exclusive_measure;
    double* cp_numcalls; double* pp_numcalls; double* vol_numcalls;
    std::array<double*,num_critical_path_measures> cp_incl_measure;
    std::array<double*,num_critical_path_measures> cp_excl_measure;
    std::array<double*,num_per_process_measures> pp_incl_measure;
    std::array<double*,num_per_process_measures> pp_excl_measure;
    std::array<double*,num_volume_measures> vol_incl_measure;
    std::array<double*,num_volume_measures> vol_excl_measure;
    bool has_been_processed;
};


extern std::string stream_name,file_name;
extern bool flag,is_first_iter,is_world_root,need_new_line,print_volume_symbol;
extern size_t mode;
extern std::ofstream stream;

extern double computation_timer;
extern std::map<MPI_Request,std::pair<std::array<MPI_Request,6>,bool>> internal_comm_info;
extern std::map<MPI_Request,std::pair<MPI_Comm,int>> internal_comm_comm;
extern std::map<MPI_Request,std::pair<int*,int*>> internal_comm_message_int;
extern std::map<MPI_Request,std::array<double*,3>> internal_comm_message_double;
extern std::map<MPI_Request,char*> internal_comm_message_char;
extern std::map<MPI_Request,std::pair<double,double>> internal_comm_data;
extern std::map<MPI_Request,nonblocking*> internal_comm_track;
extern bool decisions[breakdown_size];
extern std::array<double,critical_path_costs_size> critical_path_costs;
extern std::array<double,per_process_costs_size> max_per_process_costs;
extern std::array<double,volume_costs_size> volume_costs;
extern std::map<std::string,std::vector<double>> save_info;
extern double new_cs[critical_path_costs_size];
extern double scratch_pad;
extern std::vector<char> synch_pad_send;
extern std::vector<char> synch_pad_recv;
extern std::array<char,max_timer_name_length*max_num_symbols> symbol_pad;
extern std::array<int,max_num_symbols> symbol_len_pad;
extern std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_local_cp;
extern std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_global_cp;
extern std::array<double,(num_ftimer_measures*num_per_process_measures+1)*max_num_symbols> symbol_timer_pad_local_pp;
extern std::array<double,(num_ftimer_measures*num_per_process_measures+1)*max_num_symbols> symbol_timer_pad_global_pp;
extern std::array<double,(num_ftimer_measures*num_volume_measures+1)*max_num_symbols> symbol_timer_pad_vol;
extern std::unordered_map<std::string,ftimer> symbol_timers;
extern std::stack<std::string> symbol_stack;
extern std::array<std::string,max_num_symbols> symbol_order;
extern double_int timer_info_sender[num_volume_measures];
extern double_int timer_info_receiver[num_volume_measures];
extern bool wait_id;
}
}

// *****************************************************************************************************************************************************************

#define TAU_START(ARG) do {\
  if (critter::internal::mode==2){\
    auto save_time = MPI_Wtime();\
    if (critter::internal::symbol_timers.find(#ARG) == critter::internal::symbol_timers.end()){\
      critter::internal::symbol_timers[#ARG] = critter::internal::ftimer(#ARG);\
      critter::internal::symbol_order[critter::internal::symbol_timers.size()-1] = #ARG;\
      critter::internal::symbol_timers[#ARG].start(save_time);\
    }\
    else{\
      critter::internal::symbol_timers[#ARG].start(save_time);\
    }}}while (0);

#define TAU_STOP(ARG) do {\
  if (critter::internal::mode==2){\
    auto save_time = MPI_Wtime();\
    if (critter::internal::symbol_timers.find(#ARG) == critter::internal::symbol_timers.end()){ assert(0); }\
    else{ critter::internal::symbol_timers[#ARG].stop(save_time); }\
    }}while (0);

// *****************************************************************************************************************************************************************

#define MPI_Init(argc, argv)\
  do {\
     assert(critter::breakdown_size == critter::breakdown.count());\
     assert(critter::cost_model_size == critter::cost_models.count());\
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
     assert(critter::breakdown_size == critter::breakdown.count());\
     assert(critter::cost_model_size == critter::cost_models.count());\
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
      critter::internal::_MPI_Barrier.start(_critter_curTime_, 0, MPI_CHAR, cm);\
      PMPI_Barrier(cm);\
      critter::internal::_MPI_Barrier.intermediate();\
      PMPI_Barrier(cm);\
      critter::internal::_MPI_Barrier.stop();\
    }\
    else{\
      PMPI_Barrier(cm);\
    }\
  } while (0)

#define MPI_Bcast(buf, nelem, t, root, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      int _critter_cm_rank; MPI_Comm_rank(cm,&_critter_cm_rank);\
      critter::internal::_MPI_Bcast.start(_critter_curTime_, nelem, t, cm, _critter_cm_rank==root);\
      PMPI_Bcast(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, root, cm);\
      critter::internal::_MPI_Bcast.intermediate();\
      PMPI_Bcast(buf, nelem, t, root, cm);\
      critter::internal::_MPI_Bcast.stop();\
    }\
    else{\
      PMPI_Bcast(buf, nelem, t, root, cm);\
    }\
  } while (0)

#define MPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      int _critter_cm_rank; MPI_Comm_rank(cm,&_critter_cm_rank);\
      critter::internal::_MPI_Reduce.start(_critter_curTime_, nelem, t, cm, _critter_cm_rank==root);\
      PMPI_Reduce(&critter::internal::synch_pad_send[0], &critter::internal::synch_pad_recv[0], 1, MPI_CHAR, MPI_MAX, root, cm);\
      critter::internal::_MPI_Reduce.intermediate();\
      PMPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm);\
      critter::internal::_MPI_Reduce.stop();\
    }\
    else{\
      PMPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm);\
    }\
  } while (0)

#define MPI_Allreduce(sbuf, rbuf, nelem, t, op, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      int _critter_cm_rank; MPI_Comm_rank(cm,&_critter_cm_rank);\
      critter::internal::_MPI_Allreduce.start(_critter_curTime_, nelem, t, cm,true);\
      PMPI_Allreduce(MPI_IN_PLACE, &critter::internal::synch_pad_send[0], 1, MPI_CHAR, MPI_MAX, cm);\
      critter::internal::_MPI_Allreduce.intermediate();\
      PMPI_Allreduce(sbuf, rbuf, nelem, t, op, cm);\
      critter::internal::_MPI_Allreduce.stop();\
    }\
    else{\
      PMPI_Allreduce(sbuf, rbuf, nelem, t, op, cm);\
    }\
  } while (0)

#define MPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      int _critter_cm_rank; MPI_Comm_rank(cm,&_critter_cm_rank);\
      int64_t _critter_recvBufferSize = std::max((int64_t)scount,(int64_t)rcount) * _critter_np;\
      critter::internal::_MPI_Gather.start(_critter_curTime_, _critter_recvBufferSize, st, cm, _critter_cm_rank==root);\
      PMPI_Gather(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, &critter::internal::synch_pad_recv[0], 1, MPI_CHAR, root, cm);\
      critter::internal::_MPI_Gather.intermediate();\
      PMPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
      critter::internal::_MPI_Gather.stop();\
    }\
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
      critter::internal::_MPI_Allgather.start(_critter_curTime_, _critter_recvBufferSize, st, cm, true);\
      PMPI_Allgather(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, &critter::internal::synch_pad_recv[0], 1, MPI_CHAR, cm);\
      critter::internal::_MPI_Allgather.intermediate();\
      PMPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm);\
      critter::internal::_MPI_Allgather.stop();\
    }\
    else{\
      PMPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm);\
    }\
  } while (0)

#define MPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      int _critter_cm_rank; MPI_Comm_rank(cm,&_critter_cm_rank);\
      int64_t _critter_sendBufferSize = std::max((int64_t)scount,(int64_t)rcount) * _critter_np;\
      critter::internal::_MPI_Scatter.start(_critter_curTime_, _critter_sendBufferSize, st, cm, _critter_cm_rank==root);\
      PMPI_Scatter(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, &critter::internal::synch_pad_recv[0], 1, MPI_CHAR, root, cm);\
      critter::internal::_MPI_Scatter.intermediate();\
      PMPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
      critter::internal::_MPI_Scatter.stop();\
    }\
    else{\
      PMPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm);\
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
      critter::internal::_MPI_Reduce_scatter.start(_critter_curTime_, _critter_tot_recv, t, cm, true);\
      PMPI_Reduce_scatter(&critter::internal::synch_pad_send[0], &critter::internal::synch_pad_recv[0], &_critter_rcounts[0], MPI_CHAR, MPI_MAX, cm);\
      critter::internal::_MPI_Reduce_scatter.intermediate();\
      PMPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm);\
      critter::internal::_MPI_Reduce_scatter.stop();\
    }\
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
      critter::internal::_MPI_Alltoall.start(_critter_curTime_,_critter_recvBufferSize, st, cm, true);\
      PMPI_Alltoall(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, &critter::internal::synch_pad_recv[0], 1, MPI_CHAR, cm);\
      critter::internal::_MPI_Alltoall.intermediate();\
      PMPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm);\
      critter::internal::_MPI_Alltoall.stop();\
    }\
    else{\
      PMPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm);\
    }\
  } while (0)

#define MPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int64_t _critter_tot_recv=0;\
      int _critter_rank, _critter_np; MPI_Comm_rank(cm, &_critter_rank); MPI_Comm_size(cm, &_critter_np);\
      std::vector<int> _critter_rcounts(_critter_np,1); std::vector<int> _critter_rdisp(_critter_np,0);\
      for (int _critter_i=1; _critter_i<_critter_np; _critter_i++) _critter_rdisp[_critter_i]=_critter_rdisp[_critter_i-1]+1;\
      if (_critter_rank == root) for (int _critter_i=0; _critter_i<_critter_np; _critter_i++){ _critter_tot_recv += ((int*)rcounts)[_critter_i]; }\
      critter::internal::_MPI_Gatherv.start(_critter_curTime_, std::max((int64_t)scount,_critter_tot_recv), st, cm, _critter_rank==root);\
      PMPI_Gatherv(&critter::internal::synch_pad_send[0], 1, st, &critter::internal::synch_pad_recv[0], &_critter_rcounts[0], &_critter_rdisp[0], rt, root, cm);\
      critter::internal::_MPI_Gatherv.intermediate();\
      PMPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm);\
      critter::internal::_MPI_Gatherv.stop();\
    }\
    else{\
      PMPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm);\
    }\
  } while (0)

#define MPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(rt==st);\
      int64_t _critter_tot_recv=0; int _critter_np; MPI_Comm_size(cm, &_critter_np);\
      int _critter_cm_rank; MPI_Comm_rank(cm,&_critter_cm_rank);\
      std::vector<int> _critter_rcounts(_critter_np,1); std::vector<int> _critter_rdisp(_critter_np,0);\
      for (int _critter_i=1; _critter_i<_critter_np; _critter_i++) _critter_rdisp[_critter_i]=_critter_rdisp[_critter_i-1]+1;\
      for (int _critter_i=0; _critter_i<_critter_np; _critter_i++){ _critter_tot_recv += rcounts[_critter_i]; }\
      critter::internal::_MPI_Allgatherv.start(_critter_curTime_, std::max((int64_t)scount,_critter_tot_recv), st, cm, true);\
      PMPI_Allgatherv(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, &critter::internal::synch_pad_recv[0], &_critter_rcounts[0], &_critter_rdisp[0], MPI_CHAR, cm);\
      critter::internal::_MPI_Allgatherv.intermediate();\
      PMPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm);\
      critter::internal::_MPI_Allgatherv.stop();\
    }\
    else{\
      PMPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm);\
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
      critter::internal::_MPI_Scatterv.start(_critter_curTime_, std::max(_critter_tot_send,(int64_t)rcount), st, cm, _critter_rank==root);\
      PMPI_Scatterv(&critter::internal::synch_pad_send[0], &_critter_scounts[0], &_critter_sdisp[0], MPI_CHAR, &critter::internal::synch_pad_recv[0], 1, MPI_CHAR, root, cm);\
      critter::internal::_MPI_Scatterv.intermediate();\
      PMPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm);\
      critter::internal::_MPI_Scatterv.stop();\
    }\
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
      critter::internal::_MPI_Alltoallv.start(_critter_curTime_, std::max(_critter_tot_send,_critter_tot_recv), st, cm, true);\
      PMPI_Alltoallv(&critter::internal::synch_pad_send[0], &_critter_counts[0], &_critter_disp[0], MPI_CHAR, &critter::internal::synch_pad_recv[0], &_critter_counts[0], &_critter_disp[0], MPI_CHAR, cm);\
      critter::internal::_MPI_Scatterv.intermediate();\
      PMPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm);\
      critter::internal::_MPI_Alltoallv.stop();\
    }\
    else{\
      PMPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm);\
    }\
  } while (0)


#define MPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(st == rt); assert(stag != critter::internal_tag); assert(rtag != critter::internal_tag);\
      critter::internal::_MPI_Sendrecv.start(_critter_curTime_, std::max(scnt,rcnt), st, cm, false, true, dest, src);\
      if (std::is_same<critter::p2p_type,critter::internal::synchronous>::value){\
        PMPI_Sendrecv(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, dest, critter::internal_tag, &critter::internal::synch_pad_recv[0], 1, MPI_CHAR, src, critter::internal_tag, cm, status);\
      }\
      critter::internal::_MPI_Sendrecv.intermediate();\
      PMPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status);\
      critter::internal::_MPI_Sendrecv.stop();\
    }\
    else{\
      PMPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status);\
    }\
  } while (0)

#define MPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(stag != critter::internal_tag); assert(rtag != critter::internal_tag);\
      critter::internal::_MPI_Sendrecv_replace.start(_critter_curTime_, scnt, st, cm, false, true, dest, src);\
      if (std::is_same<critter::p2p_type,critter::internal::synchronous>::value){\
        PMPI_Sendrecv_replace(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, dest, critter::internal_tag, src, critter::internal_tag, cm, status);\
      }\
      critter::internal::_MPI_Sendrecv_replace.intermediate();\
      PMPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status);\
      critter::internal::_MPI_Sendrecv_replace.stop();\
    }\
    else{\
      PMPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status);\
    }\
  } while (0)

#define MPI_Ssend(buf, nelem, t, dest, tag, cm)\
  do {\
    if (critter::internal::mode>=1){\
      volatile double _critter_curTime_ = MPI_Wtime();\
      assert(tag != critter::internal_tag);\
      critter::internal::_MPI_Ssend.start(_critter_curTime_, nelem, t, cm, false, true, dest);\
      PMPI_Ssend(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, dest, critter::internal_tag, cm);\
      critter::internal::_MPI_Ssend.intermediate();\
      PMPI_Ssend(buf, nelem, t, dest, tag, cm);\
      critter::internal::_MPI_Ssend.stop();\
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
      critter::internal::_MPI_Send.start(_critter_curTime_, nelem, t, cm, false, true, dest);\
      if (std::is_same<critter::p2p_type,critter::internal::synchronous>::value){\
            PMPI_Send(&critter::internal::synch_pad_send[0], 1, MPI_CHAR, dest, critter::internal_tag, cm);\
      }\
      critter::internal::_MPI_Send.intermediate();\
      PMPI_Send(buf, nelem, t, dest, tag, cm);\
      critter::internal::_MPI_Send.stop();\
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
      critter::internal::_MPI_Recv.start(_critter_curTime_, nelem, t, cm, false, false, src);\
      if (std::is_same<critter::p2p_type,critter::internal::synchronous>::value){\
        PMPI_Recv(&critter::internal::synch_pad_recv[0], 1, MPI_CHAR, src, critter::internal_tag, cm, status);\
      }\
      critter::internal::_MPI_Recv.intermediate();\
      PMPI_Recv(buf, nelem, t, src, tag, cm, status);\
      critter::internal::_MPI_Recv.stop();\
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
      critter::internal::_MPI_Isend.start(_critter_curTime_, nelem, t, cm, req, true, dest);\
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
      critter::internal::_MPI_Irecv.start(_critter_curTime_, nelem, t, cm, req, false, src);\
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
      critter::internal::_MPI_Ibcast.start(_critter_curTime_, nelem, t, cm, req);\
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
      critter::internal::_MPI_Iallreduce.start(_critter_curTime_, nelem, t, cm, req);\
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
      critter::internal::_MPI_Iallreduce.start(_critter_curTime_, nelem, t, cm, req);\
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
      critter::internal::_MPI_Igather.start(_critter_curTime_, _critter_recvBufferSize, st, cm, req);\
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
      critter::internal::_MPI_Igatherv.start(_critter_curTime_, std::max((int64_t)scount,_critter_tot_recv), st, cm, req);\
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
      critter::internal::_MPI_Iallgather.start(_critter_curTime_, _critter_recvBufferSize, st, cm, req);\
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
      critter::internal::_MPI_Iallgatherv.start(_critter_curTime_, std::max((int64_t)scount,_critter_tot_recv), st, cm, req);\
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
      critter::internal::_MPI_Iscatter.start(_critter_curTime_, _critter_sendBufferSize, st, cm, req);\
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
      critter::internal::_MPI_Iscatterv.start(_critter_curTime_, std::max(_critter_tot_send,(int64_t)rcount), st, cm, req);\
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
      critter::internal::_MPI_Ireduce_scatter.start(_critter_curTime_, _critter_tot_recv, t, cm, req);\
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
      critter::internal::_MPI_Ialltoall.start(_critter_curTime_, std::max((int64_t)scount,(int64_t)rcount)*_critter_np, st, cm, req);\
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
      critter::internal::_MPI_Ialltoallv.start(_critter_curTime_, std::max(_critter_tot_send,_critter_tot_recv), st, cm, req);\
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
      _critter_comm_track_it->second->stop(&_critter_save_request, _critter_save_comp_time, _critter_save_comm_time);\
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
      _critter_comm_track_it->second->stop(&_critter_request, _critter_save_comp_time, _critter_save_comm_time);\
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
        _critter_comm_track_it->second->stop(_critter_request, _critter_save_comp_time, _critter_save_comm_time);\
      }\
    }\
    else{\
      PMPI_Waitsome(incnt, reqs, outcnt, indices, stats);\
    }\
  } while (0)

#define MPI_Waitall(cnt, reqs, stats)\
  do {\
    if (critter::internal::mode>=1){\
      critter::internal::wait_id=true;\
      int _critter_indx; MPI_Status _critter_stat;\
      for (int _critter_i=0; _critter_i<cnt; _critter_i++){\
        MPI_Waitany(cnt, reqs, &_critter_indx, &_critter_stat);\
        if (_critter_i==0){critter::internal::wait_id=false;}\
        if ((MPI_Status*)stats != (MPI_Status*)MPI_STATUSES_IGNORE) ((MPI_Status*)stats)[_critter_indx] = _critter_stat;\
      }\
      critter::internal::wait_id=true;\
    }\
    else{\
      PMPI_Waitall(cnt, reqs, stats);\
    }\
  } while (0)

#endif /*CRITTER_H_*/
