#include <limits.h>

#include "util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace skeletonization{

std::map<comm_kernel_key,kernel_key_id> comm_kernel_map;
std::map<comp_kernel_key,kernel_key_id> comp_kernel_map;
std::vector<std::pair<comm_kernel_key,int>> comm_kernel_select_sort_list;
std::vector<std::pair<comp_kernel_key,int>> comp_kernel_select_sort_list;
std::vector<int> active_kernels;
std::vector<comm_kernel_key> active_comm_kernel_keys;
std::vector<comp_kernel_key> active_comp_kernel_keys;
volatile double comm_intercept_overhead_stage1;
volatile double comm_intercept_overhead_stage2;
volatile double comp_intercept_overhead;
size_t num_critical_path_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime, CompTime, RunTime
size_t num_per_process_measures;		// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, RunTime
size_t num_volume_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, RunTime
size_t num_tracker_critical_path_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime
size_t num_tracker_per_process_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime
size_t num_tracker_volume_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime
size_t critical_path_costs_size;
size_t per_process_costs_size;
size_t volume_costs_size;
std::vector<double> critical_path_costs;
std::vector<double> max_per_process_costs;
std::vector<double> volume_costs;
std::vector<double_int> info_sender;
std::vector<double_int> info_receiver;
int internal_tag;
int internal_tag1;
int internal_tag2;
int internal_tag3;
int internal_tag4;
int internal_tag5;
bool is_first_iter;

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

  // Communication kernel time, computation kernel time, computation time, execution time
  num_critical_path_measures 		= 1;
  num_per_process_measures 		= 1;
  num_volume_measures 			= 1;

  critical_path_costs_size            	= num_critical_path_measures;
  per_process_costs_size              	= num_per_process_measures;
  volume_costs_size                   	= num_volume_measures;

  critical_path_costs.resize(critical_path_costs_size);
  max_per_process_costs.resize(per_process_costs_size);
  volume_costs.resize(volume_costs_size);
  info_sender.resize(num_critical_path_measures);
  info_receiver.resize(num_critical_path_measures);
}

void open_symbol(const char* symbol, double curtime){}

void close_symbol(const char* symbol, double curtime){}

void final_accumulate(MPI_Comm comm, double last_time){
  max_per_process_costs = volume_costs;// copy over the per-process measurements that exist in volume_costs
  double temp_costs[5];
  for (auto i=0; i<critical_path_costs.size(); i++) temp_costs[i] = critical_path_costs[i];
  for (auto i=0; i<max_per_process_costs.size(); i++) temp_costs[critical_path_costs.size()+i] = max_per_process_costs[i];
  temp_costs[critical_path_costs.size()+max_per_process_costs.size()] = comm_intercept_overhead_stage1;
  temp_costs[critical_path_costs.size()+max_per_process_costs.size()+1] = comm_intercept_overhead_stage2;
  temp_costs[critical_path_costs.size()+max_per_process_costs.size()+2] = comp_intercept_overhead;
  PMPI_Allreduce(MPI_IN_PLACE,&temp_costs[0],5,MPI_DOUBLE,MPI_MAX,comm);
  for (auto i=0; i<critical_path_costs.size(); i++) critical_path_costs[i] = temp_costs[i];
  for (auto i=0; i<max_per_process_costs.size(); i++) max_per_process_costs[i] = temp_costs[critical_path_costs.size()+i];
  comm_intercept_overhead_stage1 = temp_costs[critical_path_costs.size()+max_per_process_costs.size()];
  comm_intercept_overhead_stage2 = temp_costs[critical_path_costs.size()+max_per_process_costs.size()+1];
  comp_intercept_overhead = temp_costs[critical_path_costs.size()+max_per_process_costs.size()+2];

  skeletonization::_MPI_Barrier.comm = MPI_COMM_WORLD;
  skeletonization::_MPI_Barrier.partner1 = -1;
  skeletonization::_MPI_Barrier.partner2 = -1;
/*
  Note: the idea below is if one wants to incrementally reduce each kernel's confidence interval length, without having to "cheat"
        and use the per-process count (which usually sufficies).
  After every process owns the same information (the kernel count along the analytical cp),
     iterate over the comp map and the comm map separately and pick out the kernels that have the most kernels or the most cost.
     Perhaps leave that choice for an environment variable?
     Write to two new maps that will be referenced in discretization if scm==3, to propagate kernel information at the reduction,
       with no need to propagate kernel key information because the order is fixed right here.
    scm==1 uses pp information for each kernel, and scm==3 should use cp information for each kernel.
  Just as a judgement call, we only want to propagate the kernel count, and nothing else.
    We want to limit the samples that make up the distributions themselves to per-process ones,
    even if the kernel count along the cp no longer matches the number of kernel schedules local to a processor.
*/
}

void reset(){
  memset(&critical_path_costs[0],0,sizeof(double)*critical_path_costs.size());
  memset(&max_per_process_costs[0],0,sizeof(double)*max_per_process_costs.size());
  memset(&volume_costs[0],0,sizeof(double)*volume_costs.size());
  active_kernels.clear();
  active_comm_kernel_keys.clear();
  active_comp_kernel_keys.clear();
  comm_kernel_map.clear();
  comp_kernel_map.clear();
  internal::bsp_counter=0;

  comm_kernel_select_sort_list.clear();
  comp_kernel_select_sort_list.clear();

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
