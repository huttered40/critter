#include "util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace decomposition{

void allocate(MPI_Comm comm){
  int _world_size; MPI_Comm_size(MPI_COMM_WORLD,&_world_size);
  cp_symbol_class_count = 4;
  pp_symbol_class_count = 4;
  vol_symbol_class_count = 4;// should truly be 2, but set to 4 to conform to pp_symbol_class_count
  mode_1_width = 25;
  mode_2_width = 15;
  event_list_size = 0;

  cost_model_size=0; symbol_path_select_size=0; comm_path_select_size=0;
  //TODO: Not a fan of these magic numbers '2' and '9'. Should utilize some error checking for strings that are not of proper length anyways.
  for (auto i=0; i<2; i++){
    if (_cost_models_[i] == '1'){ cost_model_size++; }
    cost_models.push_back(_cost_models_[i]);
  } 
  for (auto i=0; i<8; i++){
    if (_symbol_path_select_[i] == '1'){ symbol_path_select_size++; symbol_path_select_index.push_back(i);}
    symbol_path_select.push_back(_symbol_path_select_[i]);
  } 
  for (auto i=0; i<8; i++){
    if (_comm_path_select_[i] == '1'){ comm_path_select_size++; }
    comm_path_select.push_back(_comm_path_select_[i]);
  } 

  num_critical_path_measures 		= 4+2*cost_model_size;// Reason for '4' instead of '5' is because we are not interested in the critical-path idle time.
  num_per_process_measures 		= 5+2*cost_model_size;
  num_volume_measures 			= 5+2*cost_model_size;
  num_tracker_critical_path_measures 	= 2+2*cost_model_size;
  num_tracker_per_process_measures 	= 2+2*cost_model_size;
  num_tracker_volume_measures 		= 2+2*cost_model_size;

  // The '2*comm_path_select_size' used below are used to track the computation time and idle time along each of the 'comm_path_select_size' paths.
  critical_path_costs_size            	= num_critical_path_measures+num_tracker_critical_path_measures*comm_path_select_size*list_size+2*comm_path_select_size;
  per_process_costs_size              	= num_per_process_measures+num_tracker_per_process_measures*comm_path_select_size*list_size+2*comm_path_select_size;
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


void final_accumulate(double last_time){
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
