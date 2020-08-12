#include "util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace discretization{

double get_arithmetic_mean(const pattern_key_id& index){
  // returns arithmetic mean
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return pattern_list[index.val_index].M1;
}
double get_variance(const pattern_key_id& index){
  // returns variance
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  size_t n = pattern_list[index.val_index].num_schedules;
  if (n<=1) return 1000.;
  return pattern_list[index.val_index].M2 / (n-1.);
}
double get_std_dev(const pattern_key_id& index){
  // returns variance
  return pow(get_variance(index),1./2.);
}
double get_std_error(const pattern_key_id& index){
  // returns standard error
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  size_t n = pattern_list[index.val_index].num_schedules;
  return get_std_dev(index) / pow(n*1.,1./2.);
}
double get_confidence_interval(const pattern_key_id& index, double level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(index);
}
/*
double get_skewness(pattern_key_id index){
  // returns skewness
  size_t n = pattern_list[index].num_schedules;
  if (n<=2) return 1000.;
  return std::pow(n*1.,1./2.) * pattern_list[index].M3 / std::pow(pattern_list[index].M2,1.5);
}
double get_kurtosis(pattern_key_id index){
  // returns the excess kurtosis
  size_t n = pattern_list[index].num_schedules;
  if (n<=3) return 1000.;
  return (n*pattern_list[index].M4) / (pattern_list[index].M2*pattern_list[index].M2) - 3.;
}
double get_jacque_barra(pattern_key_id index){
  // Note this test may be too sensitive for our purposes. But we can try this first.
  size_t n = pattern_list[index].num_schedules;
  double k = get_kurtosis();
  double s = get_skewness();
  double jb_test_stat = (n/6.)*(s*s + (1./4.)*(k*k));
  return jb_test_stat;
}
*/
void error_test(const pattern_key_id& index){
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  //assert(index.is_active);
  bool decision = ((get_confidence_interval(index) / (2.*get_arithmetic_mean(index))) < pattern_error_limit) &&
                  (pattern_list[index.val_index].num_schedules >= pattern_count_limit) &&
                  (pattern_list[index.val_index].total_exec_time >= pattern_time_limit);
  if (decision){
    pattern_list[index.val_index].steady_state = 1;
  } else{
    pattern_list[index.val_index].steady_state = 0;
  }
}
int should_schedule(const pattern_key_id& index){
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return (pattern_list[index.val_index].steady_state==1 ? 0 : 1);
}
int should_schedule_global(const pattern_key_id& index){
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return (pattern_list[index.val_index].global_steady_state==1 ? 0 : 1);
}
void set_schedule(const pattern_key_id& index, bool schedule_decision){
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  //assert(index.is_active);
  pattern_list[index.val_index].steady_state = (schedule_decision==true ? 0 : 1);
  pattern_list[index.val_index].global_steady_state = (schedule_decision==true ? 0 : 1);
}
void update(const pattern_key_id& index, volatile double exec_time, double unit_count){
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  if (pattern_list[index.val_index].steady_state == 0){
    pattern_list[index.val_index].num_schedules++;
    pattern_list[index.val_index].num_scheduled_units += unit_count;
    pattern_list[index.val_index].total_exec_time += exec_time;
//    pattern_list[index].save_exec_times.push_back((double)exec_time);	// only temporary
    // Online computation of up to 4th-order central moments using compunication time samples
    size_t n1 = pattern_list[index.val_index].num_schedules-1;
    size_t n = pattern_list[index.val_index].num_schedules;
    double x = exec_time;
    double delta = x - pattern_list[index.val_index].M1;
    double delta_n = delta / n;
    double delta_n2 = delta_n*delta_n;
    double term1 = delta*delta_n*n1;
    pattern_list[index.val_index].M1 += delta_n;
//    this->M4 = this->M4 + term1*delta_n2*(n*n - 3*n + 3) + 6*delta_n2*this->M2-4*delta_n*this->M3;
//    this->M3 = this->M3 + term1*delta_n*(n-2)-3*delta_n*this->M2;
    pattern_list[index.val_index].M2 += term1;
    error_test(index);
  }
  else{
    pattern_list[index.val_index].num_non_schedules++;
    pattern_list[index.val_index].num_non_scheduled_units += unit_count;
  }
}

void allocate(MPI_Comm comm){
  int _world_size; MPI_Comm_size(MPI_COMM_WORLD,&_world_size);
  int _world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&_world_rank);
  mode_1_width = 25;
  mode_2_width = 15;

  communicator_map[MPI_COMM_WORLD] = std::make_pair(_world_size,0);

  comp_pattern_param1_key ex_1;
  MPI_Datatype comp_pattern_key_internal_type[2] = { MPI_INT, MPI_DOUBLE };
  int comp_pattern_key_internal_type_block_len[2] = { 7,1 };
  MPI_Aint comp_pattern_key_internal_type_disp[2] = { (char*)&ex_1.tag-(char*)&ex_1, (char*)&ex_1.flops-(char*)&ex_1 };
  PMPI_Type_create_struct(2,comp_pattern_key_internal_type_block_len,comp_pattern_key_internal_type_disp,comp_pattern_key_internal_type,&comp_pattern_key_type);
  PMPI_Type_commit(&comp_pattern_key_type);

  comm_pattern_param1_key ex_2;
  MPI_Datatype comm_pattern_key_internal_type[2] = { MPI_INT, MPI_DOUBLE };
  int comm_pattern_key_internal_type_block_len[2] = { 5,1 };
  MPI_Aint comm_pattern_key_internal_type_disp[2] = { (char*)&ex_2.tag-(char*)&ex_2, (char*)&ex_2.msg_size-(char*)&ex_2 };
  PMPI_Type_create_struct(2,comm_pattern_key_internal_type_block_len,comm_pattern_key_internal_type_disp,comm_pattern_key_internal_type,&comm_pattern_key_type);
  PMPI_Type_commit(&comm_pattern_key_type);

  pattern_param1 ex_3;
  MPI_Datatype pattern_internal_type[2] = { MPI_INT, MPI_DOUBLE };
  int pattern_internal_block_len[2] = { 4,5 };
  MPI_Aint pattern_internal_disp[2] = { (char*)&ex_3.steady_state-(char*)&ex_3, (char*)&ex_3.num_scheduled_units-(char*)&ex_3 };
  PMPI_Type_create_struct(2,pattern_internal_block_len,pattern_internal_disp,pattern_internal_type,&pattern_type);
  PMPI_Type_commit(&pattern_type);

  cost_model_size=0; symbol_path_select_size=0; comm_path_select_size=0;
  //TODO: Not a fan of these magic numbers '2' and '9'. Should utilize some error checking for strings that are not of proper length anyways.

  num_critical_path_measures 		= 5;
  num_per_process_measures 		= 6;
  num_volume_measures 			= 6;

  // The '3*comm_path_select_size' used below are used to track {computation cost, computation time, idle time} along each of the 'comm_path_select_size' paths.
  critical_path_costs_size            	= num_critical_path_measures;
  per_process_costs_size              	= num_per_process_measures;
  volume_costs_size                   	= num_volume_measures;

  decisions.resize(comm_path_select_size);
  critical_path_costs.resize(critical_path_costs_size);
  max_per_process_costs.resize(per_process_costs_size);
  volume_costs.resize(volume_costs_size);
  new_cs.resize(critical_path_costs_size);
  info_sender.resize(num_critical_path_measures);
  info_receiver.resize(num_critical_path_measures);

  if (eager_p2p){
    int eager_msg_sizes[8];
    MPI_Pack_size(1,MPI_CHAR,comm,&eager_msg_sizes[0]);
    MPI_Pack_size(1,MPI_CHAR,comm,&eager_msg_sizes[1]);
    MPI_Pack_size(num_critical_path_measures,MPI_DOUBLE_INT,comm,&eager_msg_sizes[2]);
    MPI_Pack_size(critical_path_costs_size,MPI_DOUBLE,comm,&eager_msg_sizes[3]);
    MPI_Pack_size(1,MPI_INT,comm,&eager_msg_sizes[4]);
    MPI_Pack_size(max_num_symbols,MPI_INT,comm,&eager_msg_sizes[5]);
    MPI_Pack_size(max_num_symbols*max_timer_name_length,MPI_CHAR,comm,&eager_msg_sizes[6]);
    MPI_Pack_size(symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*max_num_symbols,MPI_DOUBLE,comm,&eager_msg_sizes[7]);
    int eager_pad_size = 8*MPI_BSEND_OVERHEAD;
    for (int i=0; i<8; i++) { eager_pad_size += eager_msg_sizes[i]; }
    eager_pad.resize(eager_pad_size);
  }
}



void open_symbol(const char* symbol, double curtime){}

void close_symbol(const char* symbol, double curtime){}

void final_accumulate(double last_time){
  critical_path_costs[num_critical_path_measures-2]+=(last_time-computation_timer);	// update critical path computation time
  critical_path_costs[num_critical_path_measures-1]+=(last_time-computation_timer);	// update critical path runtime
  volume_costs[num_volume_measures-2]+=(last_time-computation_timer);			// update computation time volume
  volume_costs[num_volume_measures-1]+=(last_time-computation_timer);			// update runtime volume
}

void reset(){
  for (auto i=0; i<list_size; i++){ list[i]->init(); }
  memset(&critical_path_costs[0],0,sizeof(double)*critical_path_costs.size());
  memset(&max_per_process_costs[0],0,sizeof(double)*max_per_process_costs.size());
  memset(&volume_costs[0],0,sizeof(double)*volume_costs.size());
  autotuning_propagate=1;// means nothing if autotuning_mode != 3
}

void clear(){}

}
}
}
