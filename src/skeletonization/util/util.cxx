#include <limits.h>

#include "util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace skeletonization{

MPI_Datatype comm_pattern_key_type;
MPI_Datatype comp_pattern_key_type;
std::map<comm_pattern_key,pattern_key_id> comm_pattern_map;
std::map<comp_pattern_key,pattern_key_id> comp_pattern_map;
std::vector<comm_pattern_key> steady_state_comm_pattern_keys;
std::vector<comm_pattern_key> active_comm_pattern_keys;
std::vector<comp_pattern_key> steady_state_comp_pattern_keys;
std::vector<comp_pattern_key> active_comp_pattern_keys;
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

  comp_pattern_key ex_1;
  MPI_Datatype comp_pattern_key_internal_type[2] = { MPI_INT, MPI_DOUBLE };
  int comp_pattern_key_internal_type_block_len[2] = { 7,1 };
  MPI_Aint comp_pattern_key_internal_type_disp[2] = { (char*)&ex_1.tag-(char*)&ex_1, (char*)&ex_1.flops-(char*)&ex_1 };
  PMPI_Type_create_struct(2,comp_pattern_key_internal_type_block_len,comp_pattern_key_internal_type_disp,comp_pattern_key_internal_type,&comp_pattern_key_type);
  PMPI_Type_commit(&comp_pattern_key_type);

  comm_pattern_key ex_2;
  MPI_Datatype comm_pattern_key_internal_type[2] = { MPI_INT, MPI_DOUBLE };
  int comm_pattern_key_internal_type_block_len[2] = { 9,1 };
  MPI_Aint comm_pattern_key_internal_type_disp[2] = { (char*)&ex_2.tag-(char*)&ex_2, (char*)&ex_2.msg_size-(char*)&ex_2 };
  PMPI_Type_create_struct(2,comm_pattern_key_internal_type_block_len,comm_pattern_key_internal_type_disp,comm_pattern_key_internal_type,&comm_pattern_key_type);
  PMPI_Type_commit(&comm_pattern_key_type);

  //TODO: Not a fan of these magic numbers '2' and '9'. Should utilize some error checking for strings that are not of proper length anyways.

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
}

void reset(){
  memset(&critical_path_costs[0],0,sizeof(double)*critical_path_costs.size());
  memset(&max_per_process_costs[0],0,sizeof(double)*max_per_process_costs.size());
  memset(&volume_costs[0],0,sizeof(double)*volume_costs.size());
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
