#ifndef CRITTER__UTIL_H_
#define CRITTER__UTIL_H_

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

extern size_t num_ftimer_measures;			// ExclusiveTime/Cost, InclusiveTime/Cost (NumCalls separate so as to avoid replication)
extern size_t mode_1_width;
extern size_t mode_2_width;
extern size_t max_num_symbols;
extern size_t max_timer_name_length;
extern std::string _cost_models_,_breakdown_;
extern size_t cost_model_size;
extern size_t breakdown_size;
extern size_t auto_capture;
extern std::vector<char> cost_models;
extern std::vector<char> breakdown;
extern size_t num_critical_path_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime, DataMvtTime, CompTime, RunTime
extern size_t num_per_process_measures;		// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, DataMvtTime, CompTime, RunTime
extern size_t num_volume_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, DataMvtTime, CompTime, RunTime
extern size_t num_tracker_critical_path_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime, DataMvtTime
extern size_t num_tracker_per_process_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime, DataMvtTime
extern size_t num_tracker_volume_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime, DataMvtTime
extern size_t critical_path_costs_size;
extern size_t per_process_costs_size;
extern size_t volume_costs_size;
extern std::string stream_name,pattern_stream_name,file_name;
extern bool flag,is_first_iter,is_world_root,need_new_line;
extern size_t mechanism,mode,stack_id;
extern std::ofstream stream,pattern_stream;
extern double computation_timer;
extern std::map<MPI_Request,bool> internal_comm_info;
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
extern std::vector<char> symbol_pad_ncp;
extern std::vector<int> symbol_len_pad_cp;
extern std::vector<int> symbol_len_pad_ncp;
extern std::vector<double> symbol_timer_pad_local_cp;
extern std::vector<double> symbol_timer_pad_global_cp;
extern std::vector<double> symbol_timer_pad_local_pp;
extern std::vector<double> symbol_timer_pad_global_pp;
extern std::vector<double> symbol_timer_pad_local_vol;
extern std::vector<double> symbol_timer_pad_global_vol;
extern std::stack<std::string> symbol_stack;
extern std::vector<std::string> symbol_order;
extern std::vector<double_int> timer_info_sender;
extern std::vector<double_int> timer_info_receiver;
extern bool wait_id,waitall_id;
extern double waitall_comp_time;
extern std::set<std::pair<int,int>> comm_pattern_table1;
extern std::set<std::pair<int,int>> p2p_table;
extern std::vector<std::pair<std::pair<int,int>,int>> comm_pattern_seq;
extern int internal_tag;
extern int internal_tag1;
extern int internal_tag2;
extern int internal_tag3;
extern int internal_tag4;
extern int internal_tag5;
}
}

#endif /*CRITTER__UTIL_H_*/
