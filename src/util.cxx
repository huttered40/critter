#include "util.h"

namespace critter{
namespace internal{

size_t num_ftimer_measures;			// ExclusiveTime/Cost, InclusiveTime/Cost (NumCalls separate so as to avoid replication)
size_t mode_1_width;
size_t mode_2_width;
size_t max_num_symbols;
size_t max_timer_name_length;
std::string _cost_models_,_breakdown_;
size_t cost_model_size;
size_t breakdown_size;
size_t auto_capture;
std::vector<char> cost_models;
std::vector<char> breakdown;
size_t num_critical_path_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime, DataMvtTime, CompTime, RunTime
size_t num_per_process_measures;		// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, DataMvtTime, CompTime, RunTime
size_t num_volume_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, DataMvtTime, CompTime, RunTime
size_t num_tracker_critical_path_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime, DataMvtTime
size_t num_tracker_per_process_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime, DataMvtTime
size_t num_tracker_volume_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime, DataMvtTime
size_t critical_path_costs_size;
size_t per_process_costs_size;
size_t volume_costs_size;
std::string stream_name,pattern_stream_name,file_name;
bool flag,is_first_iter,is_world_root,need_new_line;
size_t mechanism,mode,stack_id;
std::ofstream stream,pattern_stream;
double computation_timer;
std::map<MPI_Request,bool> internal_comm_info;
std::map<MPI_Request,std::pair<MPI_Comm,int>> internal_comm_comm;
std::map<MPI_Request,std::pair<double,double>> internal_comm_data;
std::vector<std::pair<double*,int>> internal_comm_prop;
std::vector<MPI_Request> internal_comm_prop_req;
std::vector<int*> internal_timer_prop_int;
std::vector<double*> internal_timer_prop_double;
std::vector<double_int*> internal_timer_prop_double_int;
std::vector<char*> internal_timer_prop_char;
std::vector<MPI_Request> internal_timer_prop_req;
std::vector<bool> decisions;
std::vector<double> critical_path_costs;
std::vector<double> max_per_process_costs;
std::vector<double> volume_costs;
std::map<std::string,std::vector<double>> save_info;
std::vector<double> new_cs;
double scratch_pad;
std::vector<char> synch_pad_send;
std::vector<char> synch_pad_recv;
std::vector<char> barrier_pad_send;
std::vector<char> barrier_pad_recv;
std::vector<char> symbol_pad_cp;
std::vector<char> symbol_pad_ncp;
std::vector<int> symbol_len_pad_cp;
std::vector<int> symbol_len_pad_ncp;
std::vector<double> symbol_timer_pad_local_cp;
std::vector<double> symbol_timer_pad_global_cp;
std::vector<double> symbol_timer_pad_local_pp;
std::vector<double> symbol_timer_pad_global_pp;
std::vector<double> symbol_timer_pad_local_vol;
std::vector<double> symbol_timer_pad_global_vol;
std::stack<std::string> symbol_stack;
std::vector<std::string> symbol_order;
std::vector<double_int> info_sender;
std::vector<double_int> info_receiver;
bool wait_id,waitall_id;
double waitall_comp_time;
std::set<std::pair<int,int>> comm_pattern_table1;
std::set<std::pair<int,int>> p2p_table;
std::vector<std::pair<std::pair<int,int>,int>> comm_pattern_seq;
int internal_tag;
int internal_tag1;
int internal_tag2;
int internal_tag3;
int internal_tag4;
int internal_tag5;
size_t track_collective;
size_t track_p2p;
}
}
