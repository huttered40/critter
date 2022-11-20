#ifndef CRITTER__UTIL_H_
#define CRITTER__UTIL_H_

#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <vector>
#include <stack>
#include <stdint.h>
#include <map>
#include <unordered_map>
#include <cmath>
#include <assert.h>
#include <climits>

namespace internal{

extern int save_wildcard_id;
extern int bsp_counter;
extern int propagate_within_timer;
extern int world_rank,debug_rank,debug;
extern volatile double computation_timer;
extern std::vector<double> wall_timer;
extern double _wall_time;
extern size_t auto_capture;
extern bool is_world_root;
extern size_t mode,stack_id;
extern std::vector<float> scratch_pad;
extern size_t profile_blas1;
extern size_t profile_blas2;
extern size_t profile_blas3;
extern size_t profile_lapack;
extern size_t profile_collective;
extern size_t profile_p2p;
extern size_t propagate_collective;
extern size_t propagate_p2p;
extern size_t execute_kernels;
extern size_t execute_kernels_max_message_size;
extern size_t eager_limit;
extern int request_id;

struct float_int{
  float_int(){first=0; second=0;}
  float_int(float one, int two){first=one; second=two;}
  float first; int second;
};

extern std::ofstream stream;
extern int cost_model,path_decomposition;
extern int internal_tag,internal_tag1,internal_tag2;
extern int internal_tag3,internal_tag4,internal_tag5;
extern bool is_first_request,is_first_iter;
extern size_t decomp_text_width,text_width,path_count,path_measure_count;
extern double comp_start_time;
// CommCost, SynchCost, CompCost,
// CommTime, ExecTime
extern size_t num_cp_measures;
// CommCost, SynchCost, CompCost,
// CommTime, ExecTime
extern size_t num_pp_measures;
// CommCost, SynchCost, CompCost,
// CommTime, ExecTime
extern size_t num_vol_measures;
extern size_t num_decomp_cp_measures;
extern size_t num_decomp_pp_measures;
extern size_t num_decomp_vol_measures;
extern size_t cp_costs_size,pp_costs_size,vol_costs_size;
extern size_t num_kernel_ds,exclusive_only,max_num_tracked_kernels;
extern std::string path_select,path_measure_select;
extern std::vector<bool> path_decisions;
extern std::vector<float> cp_costs,max_pp_costs,vol_costs;
extern std::vector<float> cp_costs_foreign;
extern std::vector<char> eager_pad;
extern std::vector<int> path_index,path_measure_index;
extern std::stack<std::string> timer_stack;
extern std::stack<bool> propagate_within_timer_stack;

void initialize(MPI_Comm comm);
void reset();
void register_timer(const char* timer_name);
void __start_timer__(const char* timer_name, double curtime, bool propagate_within, MPI_Comm cm);
void __stop_timer__(const char* timer_name, double curtime, MPI_Comm cm);
void update_time(double last_time);
void collect_volumetric_statistics(MPI_Comm comm);
void clear();
void finalize();
}

#endif /*CRITTER__UTIL_H_*/
