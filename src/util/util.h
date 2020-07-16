#ifndef CRITTER__UTIL__UTIL_H_
#define CRITTER__UTIL__UTIL_H_

#include <mpi.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <vector>
#include <stack>
#include <stdint.h>
#include <functional>
#include <map>
#include <set>
#include <unordered_map>
#include <cmath>
#include <assert.h>

namespace critter{
namespace internal{

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

struct event{
  event(std::string _kernel, std::vector<double> _measurements){
    tag = -1; measurements = _measurements;
    kernel = _kernel;
  }
  event(std::string _kernel, std::vector<double> _measurements, int _tag, MPI_Comm _comm, int _partner1, int _partner2, bool _is_sender, bool _is_eager){
    measurements = _measurements;
    tag = _tag; comm = _comm; partner1 = _partner1; partner2 = _partner2; is_sender = _is_sender; is_eager = _is_eager;
    kernel = _kernel;
  }
  event(std::string _kernel, std::vector<double> _measurements, int _tag, MPI_Comm _comm, int _partner1, bool _is_sender, bool _is_eager, int _match_id, bool blah){
    measurements = _measurements;
    tag = _tag; comm = _comm; partner1 = _partner1; is_sender = _is_sender; is_eager = _is_eager; match_id = _match_id; is_close = false;;
    kernel = _kernel;
  }
  event(std::string _kernel,std::vector<double> _measurements, std::vector<int> _match_vec){
    measurements = _measurements;
    tag=18; match_vec = _match_vec; match_size = _match_vec.size(); is_close = true;
    kernel = _kernel;
  }
  MPI_Comm comm;
  int tag,partner1,partner2,match_size,match_id;
  bool is_sender,is_eager,is_close;
  std::vector<int> match_vec;
  std::vector<double> measurements;
  std::string kernel;
};

extern size_t cp_symbol_class_count;
extern size_t pp_symbol_class_count;
extern size_t vol_symbol_class_count;
extern size_t mode_1_width;
extern size_t mode_2_width;
extern size_t max_num_symbols;
extern size_t max_timer_name_length;
extern std::string _cost_models_,_symbol_path_select_,_comm_path_select_;
extern size_t cost_model_size;
extern size_t symbol_path_select_size;
extern size_t comm_path_select_size;
extern size_t auto_capture;
extern std::vector<char> cost_models;
extern std::vector<char> symbol_path_select;
extern std::vector<char> comm_path_select;
extern size_t num_critical_path_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime, CompTime, RunTime
extern size_t num_per_process_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, RunTime
extern size_t num_volume_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, RunTime
extern size_t num_tracker_critical_path_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t num_tracker_per_process_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t num_tracker_volume_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t critical_path_costs_size;
extern size_t per_process_costs_size;
extern size_t volume_costs_size;
extern std::string stream_name,file_name;
extern bool flag,is_first_iter,is_world_root,need_new_line,opt;
extern size_t mechanism,mode,stack_id;
extern std::ofstream stream;
extern double computation_timer;
extern std::map<MPI_Request,std::pair<bool,int>> internal_comm_info;
extern std::map<MPI_Request,std::pair<MPI_Comm,int>> internal_comm_comm;
extern std::map<MPI_Request,std::pair<double,double>> internal_comm_data;
extern std::vector<std::pair<double*,int>> internal_comm_prop;
extern std::vector<MPI_Request> internal_comm_prop_req;
extern std::vector<int*> internal_timer_prop_int;
extern std::vector<double*> internal_timer_prop_double;
extern std::vector<double_int*> internal_timer_prop_double_int;
extern std::vector<char*> internal_timer_prop_char;
extern std::vector<MPI_Request> internal_timer_prop_req;
extern std::vector<bool> decisions;
extern std::vector<double> critical_path_costs;
extern std::vector<double> max_per_process_costs;
extern std::vector<double> volume_costs;
extern std::map<std::string,std::vector<double>> save_info;
extern std::vector<double> new_cs;
extern double scratch_pad;
extern std::vector<char> synch_pad_send;
extern std::vector<char> synch_pad_recv;
extern std::vector<char> barrier_pad_send;
extern std::vector<char> barrier_pad_recv;
extern std::vector<char> symbol_pad_cp;
extern std::vector<char> symbol_pad_ncp1;
extern std::vector<char> symbol_pad_ncp2;
extern std::vector<int> symbol_len_pad_cp;
extern std::vector<int> symbol_len_pad_ncp1;
extern std::vector<int> symbol_len_pad_ncp2;
extern std::vector<double> symbol_timer_pad_local_cp;
extern std::vector<double> symbol_timer_pad_global_cp;
extern std::vector<double> symbol_timer_pad_global_cp2;
extern std::vector<double> symbol_timer_pad_local_pp;
extern std::vector<double> symbol_timer_pad_global_pp;
extern std::vector<double> symbol_timer_pad_local_vol;
extern std::vector<double> symbol_timer_pad_global_vol;
extern std::stack<std::string> symbol_stack;
extern std::vector<std::string> symbol_order;
extern std::vector<double_int> info_sender;
extern std::vector<double_int> info_receiver;
extern std::vector<int> symbol_path_select_index;
extern bool wait_id;
extern int internal_tag;
extern int internal_tag1;
extern int internal_tag2;
extern int internal_tag3;
extern int internal_tag4;
extern int internal_tag5;
extern size_t track_collective;
extern size_t track_p2p;
extern size_t track_p2p_idle;
extern size_t eager_p2p;
extern size_t delete_comm;
extern std::vector<char> eager_pad;
extern std::vector<event> event_list;
extern std::vector<int> opt_req_match;
extern std::vector<double> opt_measure_match;
extern size_t event_list_size;
extern size_t opt_max_iter;
extern size_t
         _MPI_Send__id,
         _MPI_Ssend__id,
         _MPI_Bsend__id,
         _MPI_Recv__id,
         _MPI_Sendrecv__id,
         _MPI_Sendrecv_replace__id,
         _MPI_Barrier__id,
         _MPI_Bcast__id,
         _MPI_Reduce__id,
         _MPI_Allreduce__id,
         _MPI_Gather__id,
         _MPI_Allgather__id,
         _MPI_Scatter__id,
         _MPI_Reduce_scatter__id,
         _MPI_Alltoall__id,
         _MPI_Gatherv__id,
         _MPI_Allgatherv__id,
         _MPI_Scatterv__id,
         _MPI_Alltoallv__id,
         _MPI_Isend__id,
         _MPI_Irecv__id,
         _MPI_Ibcast__id,
         _MPI_Iallreduce__id,
         _MPI_Ireduce__id,
         _MPI_Igather__id,
         _MPI_Igatherv__id,
         _MPI_Iallgather__id,
         _MPI_Iallgatherv__id,
         _MPI_Iscatter__id,
         _MPI_Iscatterv__id,
         _MPI_Ireduce_scatter__id,
         _MPI_Ialltoall__id,
         _MPI_Ialltoallv__id;

}
}

#endif /*CRITTER__UTIL__UTIL_H_*/
