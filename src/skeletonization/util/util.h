#ifndef CRITTER__SKELETONIZATION__UTIL__UTIL_H_
#define CRITTER__SKELETONIZATION__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace skeletonization{

// ****************************************************************************************************************************************************
extern std::map<comm_kernel_key,kernel_key_id> comm_kernel_map;
extern std::map<comp_kernel_key,kernel_key_id> comp_kernel_map;
extern std::vector<std::pair<comm_kernel_key,int>> comm_kernel_select_sort_list;
extern std::vector<std::pair<comp_kernel_key,int>> comp_kernel_select_sort_list;
extern std::vector<int> active_kernels;
extern std::vector<comm_kernel_key> active_comm_kernel_keys;
extern std::vector<comp_kernel_key> active_comp_kernel_keys;
extern volatile double comm_intercept_overhead_stage1;
extern volatile double comm_intercept_overhead_stage2;
extern volatile double comp_intercept_overhead;
extern volatile double comp_start_time;
extern size_t num_critical_path_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime, CompTime, RunTime
extern size_t num_per_process_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, RunTime
extern size_t num_volume_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, RunTime
extern size_t num_tracker_critical_path_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t num_tracker_per_process_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t num_tracker_volume_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t critical_path_costs_size;
extern size_t per_process_costs_size;
extern size_t volume_costs_size;
extern std::vector<double> critical_path_costs;
extern std::vector<double> max_per_process_costs;
extern std::vector<double> volume_costs;
extern std::vector<double_int> info_sender;
extern std::vector<double_int> info_receiver;
extern int internal_tag;
extern int internal_tag1;
extern int internal_tag2;
extern int internal_tag3;
extern int internal_tag4;
extern int internal_tag5;
extern bool is_first_iter;

void allocate(MPI_Comm comm);
void open_symbol(const char* symbol, double curtime);
void close_symbol(const char* symbol, double curtime);
void final_accumulate(MPI_Comm comm, double last_time);
void reset();
void clear();
void finalize();

}
}
}

#endif /*CRITTER__SKELETONIZATION__UTIL__UTIL_H_*/
