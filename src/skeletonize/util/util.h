#ifndef CRITTER__SKELETONIZE__UTIL__UTIL_H_
#define CRITTER__SKELETONIZE__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace skeletonize{

// ****************************************************************************************************************************************************
extern int skeleton_type;
extern std::map<comm_kernel_key,kernel_key_id> comm_kernel_map;
extern std::map<comp_kernel_key,kernel_key_id> comp_kernel_map;
extern std::vector<int> active_kernels;
extern std::vector<comm_kernel_key> active_comm_kernel_keys;
extern std::vector<comp_kernel_key> active_comp_kernel_keys;
extern volatile double comm_intercept_overhead_stage1;
extern volatile double comm_intercept_overhead_stage2;
extern volatile double comp_intercept_overhead;
extern volatile double comp_start_time;
extern size_t num_cp_measures,num_pp_measures;
extern size_t num_vol_measures,num_tracker_cp_measures;
extern size_t num_tracker_pp_measures,num_tracker_vol_measures;
extern size_t cp_costs_size,pp_costs_size,vol_costs_size;
extern std::vector<float> cp_costs;
extern std::vector<float> cp_costs_foreign;
extern std::vector<float> max_pp_costs;
extern std::vector<float> vol_costs;
extern std::vector<float_int> info_sender;
extern std::vector<float_int> info_receiver;
extern int internal_tag,internal_tag1,internal_tag2,internal_tag3;
extern int internal_tag4, internal_tag5;
extern bool is_first_request,is_first_iter;
extern std::vector<char> eager_pad;
extern char barrier_pad_send,barrier_pad_recv;

void allocate(MPI_Comm comm);
void init_symbol(std::vector<std::string>& symbols);
void open_symbol(const char* symbol, double curtime);
void close_symbol(const char* symbol, double curtime);
void final_accumulate(MPI_Comm comm, double last_time);
void reset();
void clear();
void finalize();

}
}
}

#endif /*CRITTER__SKELETONIZE__UTIL__UTIL_H_*/
