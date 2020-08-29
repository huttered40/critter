#include "util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace discretization{

int analysis_mode;
int is_optimized;
int autotuning_propagate;
int schedule_kernels;
int update_analysis;
int comm_sample_include_idle;
MPI_Datatype comm_pattern_key_type;
MPI_Datatype comp_pattern_key_type;
MPI_Datatype pattern_type;
size_t pattern_count_limit;
double pattern_time_limit;
double pattern_error_limit;
std::map<MPI_Comm,std::pair<int,int>> communicator_map;
std::map<comm_pattern_key,pattern_key_id> comm_pattern_map;
std::map<comp_pattern_key,pattern_key_id> comp_pattern_map;
std::vector<comm_pattern_key> steady_state_comm_pattern_keys;
std::vector<comm_pattern_key> active_comm_pattern_keys;
std::vector<comp_pattern_key> steady_state_comp_pattern_keys;
std::vector<comp_pattern_key> active_comp_pattern_keys;
std::vector<pattern> steady_state_patterns;
std::vector<pattern> active_patterns;
comm_pattern_key previous_comm_key;
std::map<std::pair<comm_pattern_key,comm_pattern_key>,idle_pattern> comm_pattern_pair_map;

bool is_key_skipable(const comm_pattern_key& key){
  // For now, only barriers cannot be skipped
  if (key.tag == 0){ return false; }
  return true;
}
bool is_key_skipable(const comp_pattern_key& key){
  return true;
}

double get_estimate(const pattern& p, int analysis_param, double unit_count){
  if (analysis_param == 0){// arithmetic mean
    return get_arithmetic_mean(p);
  } else{
    return unit_count*get_harmonic_mean(p);
  }
}
double get_estimate(const pattern_key_id& index, int analysis_param, double unit_count){
  if (analysis_param == 0){// arithmetic mean
    return get_arithmetic_mean(index);
  } else{
    return unit_count*get_harmonic_mean(index);
  }
}

double get_arithmetic_mean(const pattern& p){
  // returns arithmetic mean
  return p.M1;
}
double get_arithmetic_mean(const pattern_key_id& index){
  // returns arithmetic mean
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return pattern_list[index.val_index].M1;
}

double get_harmonic_mean(const pattern& p){
  // returns arithmetic mean
  return 1./p.M1;
}
double get_harmonic_mean(const pattern_key_id& index){
  // returns arithmetic mean
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return 1./pattern_list[index.val_index].M1;
}

double get_variance(const pattern& p, int analysis_param){
  // returns variance
  size_t n = p.num_schedules;
  if (n<=1) return 1000.;
  if (analysis_param == 0){
    return p.M2 / (n-1.);
  } else{
    return 1./p.M2 / (n-1.);
  }
}
double get_variance(const pattern_key_id& index, int analysis_param){
  // returns variance
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  size_t n = pattern_list[index.val_index].num_schedules;
  if (n<=1) return 1000.;
  if (analysis_param == 0){
    return pattern_list[index.val_index].M2 / (n-1.);
  } else{
    return 1./pattern_list[index.val_index].M2 / (n-1.);
  }
}

double get_std_dev(const pattern& p, int analysis_param){
  // returns variance
  return pow(get_variance(p,analysis_param),1./2.);
}
double get_std_dev(const pattern_key_id& index, int analysis_param){
  // returns variance
  return pow(get_variance(index,analysis_param),1./2.);
}

double get_std_error(const pattern& p, int analysis_param){
  // returns standard error
  size_t n = p.num_schedules;
  return get_std_dev(p,analysis_param) / pow(n*1.,1./2.);
}
double get_std_error(const pattern_key_id& index, int analysis_param){
  // returns standard error
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  size_t n = pattern_list[index.val_index].num_schedules;
  return get_std_dev(index,analysis_param) / pow(n*1.,1./2.);
}

double get_confidence_interval(const pattern& p, int analysis_param, double level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(p,analysis_param);
}
double get_confidence_interval(const pattern_key_id& index, int analysis_param, double level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(index,analysis_param);
}

bool is_steady(const pattern& p, int analysis_param){
  return ((get_confidence_interval(p,analysis_param) / (2.*get_estimate(p,analysis_param))) < pattern_error_limit) &&
          (p.num_schedules >= pattern_count_limit) &&
          (p.total_exec_time >= pattern_time_limit);
}
bool is_steady(const pattern_key_id& index, int analysis_param){
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return ((get_confidence_interval(index,analysis_param) / (2.*get_estimate(index,analysis_param))) < pattern_error_limit) &&
          (pattern_list[index.val_index].num_schedules >= pattern_count_limit) &&
          (pattern_list[index.val_index].total_exec_time >= pattern_time_limit);
}

bool steady_test(const comm_pattern_key& key, const pattern& p, int analysis_param){
  if (!is_key_skipable(key)) return false;
  return is_steady(p,analysis_param);
}
bool steady_test(const comm_pattern_key& key, const pattern_key_id& index, int analysis_param){
  if (!is_key_skipable(key)) return false;
  return is_steady(index,analysis_param);
}
bool steady_test(const comp_pattern_key& key, const pattern& p, int analysis_param){
  if (!is_key_skipable(key)) return false;
  return is_steady(p,analysis_param);
}
bool steady_test(const comp_pattern_key& key, const pattern_key_id& index, int analysis_param){
  if (!is_key_skipable(key)) return false;
  return is_steady(index,analysis_param);
}

void update_kernel_stats(pattern& p, int analysis_param, volatile double exec_time, double unit_count){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  if (exec_time == 0) { exec_time=1.e-9; }
  if (p.steady_state == 0){
    p.num_schedules++;
    p.num_scheduled_units += unit_count;
    p.num_propagations++;
    p.total_exec_time += exec_time;
    // Online computation of up to 4th-order central moments using compunication time samples
    size_t n1 = p.num_schedules-1;
    size_t n = p.num_schedules;
    double x;
    if (analysis_param == 0){x = exec_time; }	// prep for arithmetic mean
    else                   {x = (unit_count>0 ? unit_count : 1.)/exec_time; }	// prep for harmonic mean
    double delta = x - p.M1;
    double delta_n = delta / n;
    double delta_n2 = delta_n*delta_n;
    double term1 = delta*delta_n*n1;
    p.M1 += delta_n;
    p.M2 += term1;
  }
  else{
    p.num_non_schedules++;
    p.num_non_scheduled_units += unit_count;
    p.num_non_propagations++;
  }
}
void update_kernel_stats(const pattern_key_id& index, int analysis_param, volatile double exec_time, double unit_count){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  if (exec_time == 0) { exec_time=1.e-9; }
  if (pattern_list[index.val_index].steady_state == 0){
    pattern_list[index.val_index].num_schedules++;
    pattern_list[index.val_index].num_scheduled_units += unit_count;
    pattern_list[index.val_index].num_propagations++;
    pattern_list[index.val_index].total_exec_time += exec_time;
    // Online computation of up to 4th-order central moments using compunication time samples
    size_t n1 = pattern_list[index.val_index].num_schedules-1;
    size_t n = pattern_list[index.val_index].num_schedules;
    double x;
    if (analysis_param == 0){x = exec_time; }	// prep for arithmetic mean
    else                   {x = (unit_count>0 ? unit_count : 1.)/exec_time; }	// prep for harmonic mean
    double delta = x - pattern_list[index.val_index].M1;
    double delta_n = delta / n;
    double delta_n2 = delta_n*delta_n;
    double term1 = delta*delta_n*n1;
    pattern_list[index.val_index].M1 += delta_n;
    pattern_list[index.val_index].M2 += term1;
  }
  else{
    pattern_list[index.val_index].num_non_schedules++;
    pattern_list[index.val_index].num_non_scheduled_units += unit_count;
    pattern_list[index.val_index].num_non_propagations++;
  }
}
void update_kernel_stats(pattern& dest, const pattern& src, int analysis_param){
  // This function will implement the parallel algorithm computing the mean and variance of two partitions
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  // Online computation of up to 4th-order central moments using compunication time samples
  size_t n1 = dest.num_schedules;
  size_t n2 = src.num_schedules;
  double delta = dest.M1 - src.M1;
  dest.M1 = (n1*dest.M1 + n2*src.M1)/(n1+n2);
  dest.M2 = dest.M2 + src.M2 + delta/(n1+n2)*delta*(n1*n2);
  dest.num_schedules += src.num_schedules;
  dest.num_scheduled_units += src.num_scheduled_units;
  dest.num_non_schedules += src.num_non_schedules;
  dest.num_non_scheduled_units += src.num_non_scheduled_units;
  dest.num_propagations += src.num_propagations;
  dest.num_non_propagations += src.num_non_propagations;
  dest.total_exec_time += src.total_exec_time;
}

void update_kernel_stats(idle_pattern& p, bool is_global_steady_state, volatile double exec_time){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  if (!is_global_steady_state){
    p.num_schedules++;
    // Online computation of up to 4th-order central moments using compunication time samples
    size_t n1 = p.num_schedules-1;
    size_t n = p.num_schedules;
    double x = exec_time;
    double delta = x - p.M1;
    double delta_n = delta / n;
    double delta_n2 = delta_n*delta_n;
    double term1 = delta*delta_n*n1;
    p.M1 += delta_n;
    p.M2 += term1;
  }
  else{
    p.num_non_schedules++;
  }
}

int should_schedule(const pattern& p){
  return (p.steady_state==1 ? 0 : 1);
}
int should_schedule(const pattern_key_id& index){
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return (pattern_list[index.val_index].steady_state==1 ? 0 : 1);
}

int should_schedule_global(const pattern& p){
  return (p.global_steady_state==1 ? 0 : 1);
}
int should_schedule_global(const pattern_key_id& index){
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return (pattern_list[index.val_index].global_steady_state==1 ? 0 : 1);
}

void set_kernel_state(pattern& p, bool schedule_decision){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  p.steady_state = (schedule_decision==true ? 0 : 1);
}
void set_kernel_state(const pattern_key_id& index, bool schedule_decision){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  pattern_list[index.val_index].steady_state = (schedule_decision==true ? 0 : 1);
}

void set_kernel_state_global(pattern& p, bool schedule_decision){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  p.global_steady_state = (schedule_decision==true ? 0 : 1);
}
void set_kernel_state_global(const pattern_key_id& index, bool schedule_decision){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  pattern_list[index.val_index].global_steady_state = (schedule_decision==true ? 0 : 1);
}


void allocate(MPI_Comm comm){
  int _world_size; MPI_Comm_size(MPI_COMM_WORLD,&_world_size);
  int _world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&_world_rank);
  mode_1_width = 25;
  mode_2_width = 15;

  communicator_map[MPI_COMM_WORLD] = std::make_pair(_world_size,0);

  comp_pattern_key ex_1;
  MPI_Datatype comp_pattern_key_internal_type[2] = { MPI_INT, MPI_DOUBLE };
  int comp_pattern_key_internal_type_block_len[2] = { 7,1 };
  MPI_Aint comp_pattern_key_internal_type_disp[2] = { (char*)&ex_1.tag-(char*)&ex_1, (char*)&ex_1.flops-(char*)&ex_1 };
  PMPI_Type_create_struct(2,comp_pattern_key_internal_type_block_len,comp_pattern_key_internal_type_disp,comp_pattern_key_internal_type,&comp_pattern_key_type);
  PMPI_Type_commit(&comp_pattern_key_type);

  comm_pattern_key ex_2;
  MPI_Datatype comm_pattern_key_internal_type[2] = { MPI_INT, MPI_DOUBLE };
  int comm_pattern_key_internal_type_block_len[2] = { 5,1 };
  MPI_Aint comm_pattern_key_internal_type_disp[2] = { (char*)&ex_2.tag-(char*)&ex_2, (char*)&ex_2.msg_size-(char*)&ex_2 };
  PMPI_Type_create_struct(2,comm_pattern_key_internal_type_block_len,comm_pattern_key_internal_type_disp,comm_pattern_key_internal_type,&comm_pattern_key_type);
  PMPI_Type_commit(&comm_pattern_key_type);

  pattern ex_3;
  MPI_Datatype pattern_internal_type[2] = { MPI_INT, MPI_DOUBLE };
  int pattern_internal_block_len[2] = { 6,5 };
  MPI_Aint pattern_internal_disp[2] = { (char*)&ex_3.steady_state-(char*)&ex_3, (char*)&ex_3.num_scheduled_units-(char*)&ex_3 };
  PMPI_Type_create_struct(2,pattern_internal_block_len,pattern_internal_disp,pattern_internal_type,&pattern_type);
  PMPI_Type_commit(&pattern_type);

  //TODO: Not a fan of these magic numbers '2' and '9'. Should utilize some error checking for strings that are not of proper length anyways.

  // Communication kernel time, computation kernel time, computation time, execution time
  num_critical_path_measures 		= 4;
  num_per_process_measures 		= 4;
  num_volume_measures 			= 4;

  critical_path_costs_size            	= num_critical_path_measures;
  per_process_costs_size              	= num_per_process_measures;
  volume_costs_size                   	= num_volume_measures;

  critical_path_costs.resize(critical_path_costs_size);
  max_per_process_costs.resize(per_process_costs_size);
  volume_costs.resize(volume_costs_size);
  new_cs.resize(critical_path_costs_size);
  info_sender.resize(num_critical_path_measures);
  info_receiver.resize(num_critical_path_measures);

  if (eager_p2p){
    int eager_msg_sizes[5];
    MPI_Pack_size(1,MPI_CHAR,comm,&eager_msg_sizes[0]);
    MPI_Pack_size(1,MPI_CHAR,comm,&eager_msg_sizes[1]);
    MPI_Pack_size(num_critical_path_measures,MPI_DOUBLE_INT,comm,&eager_msg_sizes[2]);
    MPI_Pack_size(critical_path_costs_size,MPI_DOUBLE,comm,&eager_msg_sizes[3]);
    MPI_Pack_size(1,MPI_INT,comm,&eager_msg_sizes[4]);
    int eager_pad_size = 5*MPI_BSEND_OVERHEAD;
    for (int i=0; i<8; i++) { eager_pad_size += eager_msg_sizes[i]; }
    eager_pad.resize(eager_pad_size);
  }
}

void open_symbol(const char* symbol, double curtime){}

void close_symbol(const char* symbol, double curtime){}

void final_accumulate(MPI_Comm comm, double last_time){
  critical_path_costs[num_critical_path_measures-1]+=(last_time-computation_timer);	// update critical path runtime
  volume_costs[num_volume_measures-1]+=(last_time-computation_timer);			// update runtime volume

  PMPI_Allreduce(MPI_IN_PLACE,&critical_path_costs[0],critical_path_costs.size(),MPI_DOUBLE,MPI_MAX,comm);
  PMPI_Allreduce(MPI_IN_PLACE,&volume_costs[0],volume_costs.size(),MPI_DOUBLE,MPI_MAX,comm);
  // Find the max per-process overhead
  double intercept_overhead[3] = {comm_intercept_overhead_stage1,comm_intercept_overhead_stage2,
                                  comp_intercept_overhead};
  PMPI_Allreduce(MPI_IN_PLACE,&intercept_overhead,3,MPI_DOUBLE,MPI_MAX,comm);
  comm_intercept_overhead_stage1 = intercept_overhead[0];
  comm_intercept_overhead_stage2 = intercept_overhead[1];
  comp_intercept_overhead = intercept_overhead[2];
}

void reset(bool track_statistical_data_override, bool schedule_kernels_override, bool force_steady_statistical_data_overide, bool update_statistical_data_overide){
  for (auto i=0; i<list_size; i++){ list[i]->init(); }
  memset(&critical_path_costs[0],0,sizeof(double)*critical_path_costs.size());
  memset(&max_per_process_costs[0],0,sizeof(double)*max_per_process_costs.size());
  memset(&volume_costs[0],0,sizeof(double)*volume_costs.size());

  // Reset these global variables, as some are updated by function arguments for convenience
  if (std::getenv("CRITTER_AUTOTUNING_MODE") != NULL){
    analysis_mode = atoi(std::getenv("CRITTER_AUTOTUNING_MODE"));
    assert(analysis_mode>0 && analysis_mode<=3);
  } else{
    analysis_mode = 0;
  }
  if (std::getenv("CRITTER_AUTOTUNING_OPTIMIZE") != NULL){
    is_optimized = atoi(std::getenv("CRITTER_AUTOTUNING_OPTIMIZE"));
    assert(is_optimized>=0 && is_optimized<=1);
  } else{
    is_optimized = 0;
  }
  if (std::getenv("CRITTER_COMM_ENVELOPE_PARAM") != NULL){
    comm_envelope_param = atoi(std::getenv("CRITTER_COMM_ENVELOPE_PARAM"));
  } else{
    comm_envelope_param = 0;
  }
  if (std::getenv("CRITTER_COMM_UNIT_PARAM") != NULL){
    comm_unit_param = atoi(std::getenv("CRITTER_COMM_UNIT_PARAM"));
  } else{
    comm_unit_param = 0;
  }
  if (std::getenv("CRITTER_COMM_ANALYSIS_PARAM") != NULL){
    comm_analysis_param = atoi(std::getenv("CRITTER_COMM_ANALYSIS_PARAM"));
  } else{
    comm_analysis_param = 0;
  }
  if (std::getenv("CRITTER_COMP_ENVELOPE_PARAM") != NULL){
    comp_envelope_param = atoi(std::getenv("CRITTER_COMP_ENVELOPE_PARAM"));
  } else{
    comp_envelope_param = 0;
  }
  if (std::getenv("CRITTER_COMP_UNIT_PARAM") != NULL){
    comp_unit_param = atoi(std::getenv("CRITTER_COMP_UNIT_PARAM"));
  } else{
    comp_unit_param = 0;
  }
  if (std::getenv("CRITTER_COMP_ANALYSIS_PARAM") != NULL){
    comp_analysis_param = atoi(std::getenv("CRITTER_COMP_ANALYSIS_PARAM"));
  } else{
    comp_analysis_param = 0;
  }
  if (std::getenv("CRITTER_PATTERN_COUNT_LIMIT") != NULL){
    pattern_count_limit = atoi(std::getenv("CRITTER_PATTERN_COUNT_LIMIT"));
  } else{
    pattern_count_limit = 1;
  }
  if (std::getenv("CRITTER_PATTERN_TIME_LIMIT") != NULL){
    pattern_time_limit = atof(std::getenv("CRITTER_PATTERN_TIME_LIMIT"));
  } else{
    pattern_time_limit = .00001;
  }
  if (std::getenv("CRITTER_PATTERN_ERROR_LIMIT") != NULL){
    pattern_error_limit = atof(std::getenv("CRITTER_PATTERN_ERROR_LIMIT"));
  } else{
    pattern_error_limit = .5;
  }
  if (std::getenv("CRITTER_SCHEDULE_KERNELS") != NULL){
    schedule_kernels = atoi(std::getenv("CRITTER_SCHEDULE_KERNELS"));
  } else{
    schedule_kernels = 1;
  }
  if (std::getenv("CRITTER_COMM_SAMPLE_INCLUDE_IDLE") != NULL){
    comm_sample_include_idle = atoi(std::getenv("CRITTER_COMM_SAMPLE_INCLUDE_IDLE"));
  } else{
    comm_sample_include_idle = 0;
  }
  if (analysis_mode>0){ analysis_mode = (track_statistical_data_override ? analysis_mode : 0); }
  if (schedule_kernels==1){ schedule_kernels = (schedule_kernels_override ? schedule_kernels : 0); }
  update_analysis = (update_statistical_data_overide ? 1 : 0);
  autotuning_propagate=1;// means nothing if analysis_mode != 3
  if (force_steady_statistical_data_overide && analysis_mode<3){// DO NOT SET TO STEADY_STATE WHEN PERFORMING CRITICAL PATH ANALYSIS
    // This branch is to be entered only after tuning a space of algorithmic parameterizations, in which the expectation is that all kernels,
    //   both comm and comp, have reached a sufficiently-predictable state (steady state). It is also only valid for per-process or volumetric analysis.
/*
    .. my problem with this now is that it breaks the invariant that kernels in global steady state must be both inactive and in the steady state buffers.
    ..   I have recently modified the routines in path.cxx to expect this, which I did not do before (and thus is the reason the code below worked for pp/vol analysis).
    .. I cannot simply flush each pattern that is not aleady in global steady state, as 
    // set all kernels into global steady state -- note this is reasonable for now,
    //   particularly as a debugging technique, but not sure if this makes sense permanently
    for (auto it : comm_pattern_map){
      set_kernel_state(it.second,false);
      set_kernel_state_global(it.second,false);
    }
    for (auto it : comp_pattern_map){
      set_kernel_state(it.second,false);
      set_kernel_schedule_global(it.second,false);
    }
*/
  }
  else{
    update_analysis=1;// Its too risky to allow no updating to pattern members when performing critical path analysis.
/*
    .. yes but how do we protect the statistics generated when tuning the parameterization space from being corrupted, even if most kernels are already steady?
    ..   because I think they would still update the "non" pattern members, right?
*/
  }
  previous_comm_key = comm_pattern_key();
}

void clear(){
  // I don't see any reason to clear the communicator map. In fact, doing so would be harmful
  comm_pattern_map.clear();
  comp_pattern_map.clear();
  steady_state_comm_pattern_keys.clear();
  active_comm_pattern_keys.clear();
  steady_state_comp_pattern_keys.clear();
  active_comp_pattern_keys.clear();
  steady_state_patterns.clear();
  active_patterns.clear();
  comm_pattern_pair_map.clear();
  comm_intercept_overhead_stage1=0;
  comm_intercept_overhead_stage2=0;
  comp_intercept_overhead=0;
  previous_comm_key = comm_pattern_key();
}

}
}
}
