#include "util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace decomposition{

size_t cp_symbol_class_count;
size_t pp_symbol_class_count;
size_t vol_symbol_class_count;
size_t max_num_symbols;
size_t max_timer_name_length;
std::string _cost_models_,_symbol_path_select_,_comm_path_select_;
size_t cost_model_size;
size_t symbol_path_select_size;
size_t comm_path_select_size;
std::vector<char> cost_models;
std::vector<char> symbol_path_select;
std::vector<char> comm_path_select;
std::vector<char> symbol_pad_cp;
std::vector<char> symbol_pad_ncp1;
std::vector<char> symbol_pad_ncp2;
std::vector<int> symbol_len_pad_cp;
std::vector<int> symbol_len_pad_ncp1;
std::vector<int> symbol_len_pad_ncp2;
std::vector<double> symbol_timer_pad_local_cp;
std::vector<double> symbol_timer_pad_global_cp;
std::vector<double> symbol_timer_pad_global_cp2;
std::vector<double> symbol_timer_pad_local_pp;
std::vector<double> symbol_timer_pad_global_pp;
std::vector<double> symbol_timer_pad_local_vol;
std::vector<double> symbol_timer_pad_global_vol;
std::stack<std::string> symbol_stack;
std::vector<std::string> symbol_order;
std::vector<int> symbol_path_select_index;
std::vector<event> event_list;
std::vector<int> opt_req_match;
std::vector<double> opt_measure_match;
size_t event_list_size;
size_t opt_max_iter;
size_t gradient_jump_size;
size_t num_gradient_points;

void allocate(MPI_Comm comm){
  int _world_size; MPI_Comm_size(MPI_COMM_WORLD,&_world_size);
  int _world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&_world_rank);
  cp_symbol_class_count = 4;
  pp_symbol_class_count = 4;
  vol_symbol_class_count = 4;// should truly be 2, but set to 4 to conform to pp_symbol_class_count
  mode_1_width = 25;
  mode_2_width = 15;
  event_list_size = 0;

/*
  if (std::getenv("CRITTER_OPT") != NULL){
    opt = atoi(std::getenv("CRITTER_OPT"));
    delete_comm = 0;
  } else{
    opt = 0;
  }
  if (std::getenv("CRITTER_OPT_MAX_ITER") != NULL){
    opt_max_iter = atoi(std::getenv("CRITTER_OPT_MAX_ITER"));
  } else{
    opt_max_iter = 5;
  }
  if (std::getenv("CRITTER_OPT_NUM_GRADIENT_POINTS") != NULL){
    num_gradient_points = atoi(std::getenv("CRITTER_OPT_NUM_GRADIENT_POINTS"));
  } else{
    num_gradient_points = 10;
  }
  if (std::getenv("CRITTER_OPT_GRADIENT_JUMP_SIZE") != NULL){
    gradient_jump_size = atoi(std::getenv("CRITTER_OPT_GRADIENT_JUMP_SIZE"));
  } else{
    gradient_jump_size = 5;// signifies 5%
  }
*/
  if (std::getenv("CRITTER_MODEL_SELECT") != NULL){
    _cost_models_ = std::getenv("CRITTER_MODEL_SELECT");
  } else{
    _cost_models_ = "11";
  }
  if (std::getenv("CRITTER_SYMBOL_PATH_SELECT") != NULL){
    _symbol_path_select_ = std::getenv("CRITTER_SYMBOL_PATH_SELECT");
  } else{
    _symbol_path_select_ = "000000000";
  }
  if (std::getenv("CRITTER_COMM_PATH_SELECT") != NULL){
    _comm_path_select_ = std::getenv("CRITTER_COMM_PATH_SELECT");
  } else{
    _comm_path_select_ = "000000000";
  }
  if (std::getenv("CRITTER_VIZ_FILE") != NULL){
    flag = 1;
    file_name = std::getenv("CRITTER_VIZ_FILE");
    stream_name = file_name + ".txt";
  }
  if (std::getenv("CRITTER_MAX_NUM_SYMBOLS") != NULL){
    max_num_symbols = atoi(std::getenv("CRITTER_MAX_NUM_SYMBOLS"));
  } else{
    max_num_symbols = 15;
  }
  if (std::getenv("CRITTER_MAX_SYMBOL_LENGTH") != NULL){
    max_timer_name_length = atoi(std::getenv("CRITTER_MAX_SYMBOL_LENGTH"));
  } else{
    max_timer_name_length = 25;
  }
  assert(_cost_models_.size()==2);
  assert(_comm_path_select_.size()==9);
  assert(_symbol_path_select_.size()==9);

  cost_model_size=0; symbol_path_select_size=0; comm_path_select_size=0;
  //TODO: Not a fan of these magic numbers '2' and '9'. Should utilize some error checking for strings that are not of proper length anyways.
  for (auto i=0; i<2; i++){
    if (_cost_models_[i] == '1'){ cost_model_size++; }
    cost_models.push_back(_cost_models_[i]);
  } 
  for (auto i=0; i<9; i++){
    if (_symbol_path_select_[i] == '1'){ symbol_path_select_size++; symbol_path_select_index.push_back(i);}
    symbol_path_select.push_back(_symbol_path_select_[i]);
  } 
  for (auto i=0; i<9; i++){
    if (_comm_path_select_[i] == '1'){ comm_path_select_size++; }
    comm_path_select.push_back(_comm_path_select_[i]);
  } 

  num_critical_path_measures 		= 5+2*cost_model_size;// Reason for '5' instead of '6' is because we are not interested in the critical-path idle time.
  num_per_process_measures 		= 6+2*cost_model_size;
  num_volume_measures 			= 6+2*cost_model_size;
  num_tracker_critical_path_measures 	= 2+2*cost_model_size;
  num_tracker_per_process_measures 	= 2+2*cost_model_size;
  num_tracker_volume_measures 		= 2+2*cost_model_size;

  // The '3*comm_path_select_size' used below are used to track {computation cost, computation time, idle time} along each of the 'comm_path_select_size' paths.
  critical_path_costs_size            	= num_critical_path_measures+num_tracker_critical_path_measures*comm_path_select_size*list_size+3*comm_path_select_size;
  per_process_costs_size              	= num_per_process_measures+num_tracker_per_process_measures*comm_path_select_size*list_size+3*comm_path_select_size;
  volume_costs_size                   	= num_volume_measures+num_tracker_volume_measures*list_size;

  synch_pad_send.resize(_world_size);
  synch_pad_recv.resize(_world_size);
  barrier_pad_send.resize(_world_size);
  barrier_pad_recv.resize(_world_size);

  decisions.resize(comm_path_select_size);
  critical_path_costs.resize(critical_path_costs_size);

  max_per_process_costs.resize(per_process_costs_size);
  volume_costs.resize(volume_costs_size);
  new_cs.resize(critical_path_costs_size);
  // The reason 'symbol_pad_cp' and 'symbol_len_pad_cp' are a factor 'symbol_path_select_size' larger than the 'ncp*'
  //   variants is because those variants are used solely for p2p, in which we simply transfer a process's path data, rather than reduce it using a special multi-root trick.
  symbol_pad_cp.resize(symbol_path_select_size*max_timer_name_length*max_num_symbols);
  symbol_pad_ncp1.resize(max_timer_name_length*max_num_symbols);
  symbol_pad_ncp2.resize(max_timer_name_length*max_num_symbols);
  symbol_len_pad_cp.resize(symbol_path_select_size*max_num_symbols);
  symbol_len_pad_ncp1.resize(max_num_symbols);
  symbol_len_pad_ncp2.resize(max_num_symbols);
  // Note: we use 'num_per_process_measures' rather than 'num_critical_path_measures' for specifying the
  //   length of 'symbol_timer_pad_*_cp' because we want to track idle time contribution of each symbol along a path.
  symbol_timer_pad_local_cp.resize(symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*max_num_symbols,0.);
  symbol_timer_pad_global_cp.resize(symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*max_num_symbols,0.);
  symbol_timer_pad_local_pp.resize((pp_symbol_class_count*num_per_process_measures+1)*max_num_symbols,0.);
  symbol_timer_pad_global_pp.resize((pp_symbol_class_count*num_per_process_measures+1)*max_num_symbols,0.);
  symbol_timer_pad_local_vol.resize((vol_symbol_class_count*num_volume_measures+1)*max_num_symbols,0.);
  symbol_timer_pad_global_vol.resize((vol_symbol_class_count*num_volume_measures+1)*max_num_symbols,0.);
  symbol_order.resize(max_num_symbols);
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

void reset(){
  for (auto i=0; i<list_size; i++){ list[i]->init(); }
  memset(&critical_path_costs[0],0,sizeof(double)*critical_path_costs.size());
  memset(&max_per_process_costs[0],0,sizeof(double)*max_per_process_costs.size());
  memset(&volume_costs[0],0,sizeof(double)*volume_costs.size());
  memset(&symbol_timer_pad_local_cp[0],0,sizeof(double)*symbol_timer_pad_local_cp.size());
  memset(&symbol_timer_pad_local_pp[0],0,sizeof(double)*symbol_timer_pad_local_pp.size());
  memset(&symbol_timer_pad_local_vol[0],0,sizeof(double)*symbol_timer_pad_local_vol.size());

  comm_intercept_overhead_stage1=0;
  comm_intercept_overhead_stage2=0;
  comm_intercept_overhead_stage3=0;
  comm_intercept_overhead_stage4=0;
  comp_intercept_overhead=0;
}

void open_symbol(const char* symbol, double curtime){
  if (symbol_timers.find(symbol) == symbol_timers.end()){
    symbol_timers[symbol] = symbol_tracker(symbol);
    symbol_order[symbol_timers.size()-1] = symbol;
    symbol_timers[symbol].start(curtime);
  }
  else{
    symbol_timers[symbol].start(curtime);
  }
}

void close_symbol(const char* symbol, double curtime){
  if (symbol_timers.find(symbol) == symbol_timers.end()){
    assert(0);
  }
  else{
    symbol_timers[symbol].stop(curtime);
  }
}


void final_accumulate(MPI_Comm comm, double last_time){
  critical_path_costs[num_critical_path_measures-2]+=(last_time-computation_timer);	// update critical path computation time
  critical_path_costs[num_critical_path_measures-1]+=(last_time-computation_timer);	// update critical path runtime
  volume_costs[num_volume_measures-2]+=(last_time-computation_timer);			// update computation time volume
  volume_costs[num_volume_measures-1]+=(last_time-computation_timer);			// update runtime volume
  // update the computation time (i.e. time between last MPI synchronization point and this function invocation) along all paths decomposed by MPI communication routine
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-1-i] += (last_time-computation_timer); }
  // Save the communication pattern
  if (opt){
    //TODO: we will assume both costs models are chosen.
    std::vector<double> measurements(num_per_process_measures,0.);
    measurements[num_per_process_measures-1]=(last_time-computation_timer);
    measurements[num_per_process_measures-2]=(last_time-computation_timer);
    event_list.push_back(event("",std::move(measurements)));
  }
}


void clear(){
  symbol_timers.clear();
}

}
}
}
