#ifndef CRITTER__DECOMPOSITION__UTIL__UTIL_H_
#define CRITTER__DECOMPOSITION__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace decomposition{

extern int replace_comp;
extern int replace_comm;
extern std::map<comp_kernel_key,std::pair<int,double>> replace_comp_map_local;
extern std::map<comm_kernel_key,std::pair<int,double>> replace_comm_map_local;
extern std::map<comp_kernel_key,std::pair<int,double>> replace_comp_map_global;
extern std::map<comm_kernel_key,std::pair<int,double>> replace_comm_map_global;
extern std::ofstream stream;
extern bool wait_id;
extern size_t mode_1_width;
extern size_t mode_2_width;
extern int internal_tag;
extern int internal_tag1;
extern int internal_tag2;
extern int internal_tag3;
extern int internal_tag4;
extern int internal_tag5;
extern bool is_first_iter;
extern size_t cp_symbol_class_count;
extern size_t pp_symbol_class_count;
extern size_t vol_symbol_class_count;
extern size_t max_num_symbols;
extern size_t max_timer_name_length;
extern std::string _cost_models_,_symbol_path_select_,_comm_path_select_;
extern size_t cost_model_size;
extern size_t symbol_path_select_size;
extern size_t comm_path_select_size;
extern std::vector<char> cost_models;
extern std::vector<char> symbol_path_select;
extern std::vector<char> comm_path_select;
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
extern std::vector<int> symbol_path_select_index;
extern std::vector<double> intercept_overhead;
extern std::vector<double> global_intercept_overhead;
extern size_t num_critical_path_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime, CompTime, CompKernelTime, RunTime
extern size_t num_per_process_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, CompKernelTime, RunTime
extern size_t num_volume_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, CompKernelTime, RunTime
extern size_t num_tracker_critical_path_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t num_tracker_per_process_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t num_tracker_volume_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t critical_path_costs_size;
extern size_t per_process_costs_size;
extern size_t volume_costs_size;
extern std::vector<double> critical_path_costs;
extern std::vector<double> max_per_process_costs;
extern std::vector<double> volume_costs;
extern std::vector<char> synch_pad_send;
extern std::vector<char> synch_pad_recv;
extern std::vector<char> barrier_pad_send;
extern std::vector<char> barrier_pad_recv;
extern std::vector<double_int> info_sender;
extern std::vector<double_int> info_receiver;
extern std::vector<char> eager_pad;
extern double comp_start_time;
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
extern std::map<std::string,std::vector<double>> save_info;
extern std::vector<double> new_cs;

void allocate(MPI_Comm comm);
void reset();
void open_symbol(const char* symbol, double curtime);
void close_symbol(const char* symbol, double curtime);
void final_accumulate(MPI_Comm comm, double last_time);
void clear();
void finalize();

}
}
}

#endif /*CRITTER__DECOMPOSITION__UTIL__UTIL_H_*/
