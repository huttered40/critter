#include <limits.h>

#include "util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace skeletonize{

int skeleton_type;
std::map<comm_kernel_key,kernel_key_id> comm_kernel_map;
std::map<comp_kernel_key,kernel_key_id> comp_kernel_map;
std::vector<int> active_kernels;
std::vector<comm_kernel_key> active_comm_kernel_keys;
std::vector<comp_kernel_key> active_comp_kernel_keys;
volatile double comm_intercept_overhead_stage1;
volatile double comm_intercept_overhead_stage2;
volatile double comp_intercept_overhead;
volatile double comp_start_time;
size_t num_cp_measures,num_pp_measures;
size_t num_vol_measures,num_tracker_cp_measures;
size_t num_tracker_pp_measures,num_tracker_vol_measures;
size_t cp_costs_size,pp_costs_size,vol_costs_size;
std::vector<float> cp_costs;
std::vector<float> cp_costs_foreign;
std::vector<float> max_pp_costs;
std::vector<float> vol_costs;
std::vector<float_int> info_sender;
std::vector<float_int> info_receiver;
int internal_tag,internal_tag1,internal_tag2,internal_tag3;
int internal_tag4, internal_tag5;
bool is_first_request,is_first_iter;
std::vector<char> eager_pad;
char barrier_pad_send,barrier_pad_recv;

void allocate(MPI_Comm comm){
  int _world_size; MPI_Comm_size(MPI_COMM_WORLD,&_world_size);
  int _world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&_world_rank);
  internal_tag = 31133;
  internal_tag1 = internal_tag+1;
  internal_tag2 = internal_tag+2;
  internal_tag3 = internal_tag+3;
  internal_tag4 = internal_tag+4;
  internal_tag5 = internal_tag+5;
  is_first_iter = true;
  comm_intercept_overhead_stage1=0;
  comm_intercept_overhead_stage2=0;
  comp_intercept_overhead=0;

  if (std::getenv("CRITTER_SKELETON_TYPE") != NULL){
    skeleton_type = atoi(std::getenv("CRITTER_SKELETON_TYPE"));
    assert(skeleton_type==0 || skeleton_type==1);
  } else{
    skeleton_type = 1;// analytic
  }

  // Estimated execution cost
  num_cp_measures = 1;
  num_pp_measures = 1;
  num_vol_measures = 1;

  // 8 key members in 'comp_kernel' and 8 key members in 'comm_kernel'
  // +1 for each is the schedule count itself
  cp_costs_size = num_cp_measures + 9*comp_kernel_select_count + 9*comm_kernel_select_count;
  pp_costs_size = num_pp_measures;
  vol_costs_size = num_vol_measures;

  cp_costs.resize(cp_costs_size);
  cp_costs_foreign.resize(cp_costs_size);
  max_pp_costs.resize(pp_costs_size);
  vol_costs.resize(vol_costs_size);

  int eager_msg_size[2];
  MPI_Pack_size(1,MPI_CHAR,comm,&eager_msg_size[0]);
  MPI_Pack_size(cp_costs_size,MPI_FLOAT,comm,&eager_msg_size[1]);
  int eager_pad_size = 2*MPI_BSEND_OVERHEAD;
  eager_pad_size += (eager_msg_size[0]+eager_msg_size[1]);
  assert(eager_pad_size <= eager_limit);
  eager_pad.resize(eager_pad_size);
}

void init_symbol(std::vector<std::string>& symbols){}

void open_symbol(const char* symbol, double curtime){}

void close_symbol(const char* symbol, double curtime){}

void final_accumulate(MPI_Comm comm, double last_time){

  if (skeleton_type==0){
    cp_costs[0] += (last_time-computation_timer);
    vol_costs[0] += (last_time-computation_timer);
  }

  max_pp_costs = vol_costs;// copy over the per-process measurements that exist in vol_costs
}

void reset(){
  memset(&cp_costs[0],0,sizeof(float)*cp_costs.size());
  memset(&cp_costs_foreign[0],0,sizeof(float)*cp_costs_foreign.size());
  memset(&max_pp_costs[0],0,sizeof(float)*max_pp_costs.size());
  memset(&vol_costs[0],0,sizeof(float)*vol_costs.size());
  active_kernels.clear();
  active_comm_kernel_keys.clear();
  active_comp_kernel_keys.clear();
  comm_kernel_map.clear();
  comp_kernel_map.clear();
  internal::bsp_counter=0;

  if (std::getenv("CRITTER_MODE") != NULL){
    internal::mode = atoi(std::getenv("CRITTER_MODE"));
  } else{
    internal::mode = 1;
  }
}

void clear(){}

void finalize(){}

}
}
}
