#include <limits.h>

#include "util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace decomposition{

int invoke_max_barrier;
int track_synchronization;
std::map<comp_kernel_key,std::pair<int,double>> comp_kernel_info;
std::map<comm_kernel_key,std::pair<int,double>> comm_kernel_info;
std::ofstream stream;
bool is_first_request;
int internal_tag;
int internal_tag1;
int internal_tag2;
int internal_tag3;
int internal_tag4;
int internal_tag5;
size_t mode_1_width;
size_t mode_2_width;
bool is_first_iter;
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
std::vector<double> intercept_overhead;
std::vector<double> global_intercept_overhead;
size_t num_critical_path_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime, CompTime, CompKernelTime, RunTime
size_t num_per_process_measures;		// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, CompKernelTime, RunTime
size_t num_volume_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, CompKernelTime, RunTime
size_t num_tracker_critical_path_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime
size_t num_tracker_per_process_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime
size_t num_tracker_volume_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime
size_t critical_path_costs_size;
size_t per_process_costs_size;
size_t volume_costs_size;
std::vector<double> critical_path_costs;
std::vector<double> max_per_process_costs;
std::vector<double> volume_costs;
std::vector<char> synch_pad_send;
std::vector<char> synch_pad_recv;
std::vector<char> barrier_pad_send;
std::vector<char> barrier_pad_recv;
std::vector<double_int> info_sender;
std::vector<double_int> info_receiver;
std::vector<char> eager_pad;
double comp_start_time;

std::vector<std::pair<double*,int>> internal_comm_prop;
std::vector<MPI_Request> internal_comm_prop_req;

std::vector<int*> internal_timer_prop_int;
std::vector<double*> internal_timer_prop_double;
std::vector<double_int*> internal_timer_prop_double_int;
std::vector<char*> internal_timer_prop_char;
std::vector<MPI_Request> internal_timer_prop_req;

std::vector<bool> decisions;
std::map<std::string,std::vector<double>> save_info;
std::vector<double> new_cs;

void allocate(MPI_Comm comm){
  int _world_size; MPI_Comm_size(MPI_COMM_WORLD,&_world_size);
  int _world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&_world_rank);
  cp_symbol_class_count = 4;
  pp_symbol_class_count = 4;
  vol_symbol_class_count = 4;// should truly be 2, but set to 4 to conform to pp_symbol_class_count
  mode_1_width = 25;
  mode_2_width = 15;
  internal_tag = 31133;
  internal_tag1 = internal_tag+1;
  internal_tag2 = internal_tag+2;
  internal_tag3 = internal_tag+3;
  internal_tag4 = internal_tag+4;
  internal_tag5 = internal_tag+5;
  is_first_iter = true;

  if (std::getenv("CRITTER_INVOKE_MAX_BARRIER") != NULL){
    invoke_max_barrier = atoi(std::getenv("CRITTER_INVOKE_MAX_BARRIER"));
  } else{
    invoke_max_barrier = 0;
  }
  if (std::getenv("CRITTER_TRACK_SYNCHRONIZATION") != NULL){
    track_synchronization = atoi(std::getenv("CRITTER_TRACK_SYNCHRONIZATION"));
  } else{
    track_synchronization = 0;
  }
  if (std::getenv("CRITTER_MODEL_SELECT") != NULL){
    _cost_models_ = std::getenv("CRITTER_MODEL_SELECT");
  } else{
    _cost_models_ = "11";
  }
  if (std::getenv("CRITTER_SYMBOL_PATH_SELECT") != NULL){
    _symbol_path_select_ = std::getenv("CRITTER_SYMBOL_PATH_SELECT");
  } else{
    _symbol_path_select_ = "0000000000";
  }
  if (std::getenv("CRITTER_COMM_PATH_SELECT") != NULL){
    _comm_path_select_ = std::getenv("CRITTER_COMM_PATH_SELECT");
  } else{
    _comm_path_select_ = "0000000000";
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
  if (std::getenv("CRITTER_VIZ_FILE") != NULL){
    std::string stream_name = std::getenv("CRITTER_VIZ_FILE");
    stream_name += "_decomposition.txt";
    if (is_world_root){
      stream.open(stream_name.c_str(),std::ofstream::out);
      assert(stream.good());
    }
  }
  assert(_cost_models_.size()==2);
  assert(_comm_path_select_.size()==10);
  assert(_symbol_path_select_.size()==10);

  cost_model_size=0; symbol_path_select_size=0; comm_path_select_size=0;
  //TODO: Not a fan of these magic numbers '2' and '9'. Should utilize some error checking for strings that are not of proper length anyways.
  for (auto i=0; i<2; i++){
    if (_cost_models_[i] == '1'){ cost_model_size++; }
    cost_models.push_back(_cost_models_[i]);
  } 
  for (auto i=0; i<10; i++){
    if (_symbol_path_select_[i] == '1'){ symbol_path_select_size++; symbol_path_select_index.push_back(i);}
    symbol_path_select.push_back(_symbol_path_select_[i]);
  } 
  for (auto i=0; i<10; i++){
    if (_comm_path_select_[i] == '1'){ comm_path_select_size++; }
    comm_path_select.push_back(_comm_path_select_[i]);
  } 

  num_critical_path_measures 		= 6+2*cost_model_size;// Reason for '5' instead of '6' is because we are not interested in the critical-path idle time.
  num_per_process_measures 		= 7+2*cost_model_size;
  num_volume_measures 			= 7+2*cost_model_size;
  num_tracker_critical_path_measures 	= 2+2*cost_model_size;
  num_tracker_per_process_measures 	= 2+2*cost_model_size;
  num_tracker_volume_measures 		= 2+2*cost_model_size;

  // The '3*comm_path_select_size' used below are used to track {computation cost, computation time, computational kernel time, idle time} along each of the 'comm_path_select_size' paths.
  critical_path_costs_size            	= num_critical_path_measures+num_tracker_critical_path_measures*comm_path_select_size*list_size+4*comm_path_select_size;
  per_process_costs_size              	= num_per_process_measures+num_tracker_per_process_measures*comm_path_select_size*list_size+4*comm_path_select_size;
  volume_costs_size                   	= num_volume_measures+num_tracker_volume_measures*list_size;

  synch_pad_send.resize(_world_size);
  synch_pad_recv.resize(_world_size);
  barrier_pad_send.resize(_world_size);
  barrier_pad_recv.resize(_world_size);

  decisions.resize(comm_path_select_size);
  critical_path_costs.resize(critical_path_costs_size);

  intercept_overhead.resize(3,0);
  global_intercept_overhead.resize(3,0);

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

  int eager_msg_sizes[8];
  MPI_Pack_size(2,MPI_CHAR,comm,&eager_msg_sizes[0]);
  MPI_Pack_size(1,MPI_DOUBLE,comm,&eager_msg_sizes[1]);
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

void reset(){
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  if (std::getenv("CRITTER_MODE") != NULL){
    internal::mode = atoi(std::getenv("CRITTER_MODE"));
  } else{
    internal::mode = 1;
  }

  for (auto i=0; i<list_size; i++){ list[i]->init(); }
  save_info.clear();
  memset(&critical_path_costs[0],0,sizeof(double)*critical_path_costs.size());
  memset(&max_per_process_costs[0],0,sizeof(double)*max_per_process_costs.size());
  memset(&volume_costs[0],0,sizeof(double)*volume_costs.size());
  memset(&symbol_timer_pad_local_cp[0],0,sizeof(double)*symbol_timer_pad_local_cp.size());
  memset(&symbol_timer_pad_local_pp[0],0,sizeof(double)*symbol_timer_pad_local_pp.size());
  memset(&symbol_timer_pad_local_vol[0],0,sizeof(double)*symbol_timer_pad_local_vol.size());
  memset(&intercept_overhead[0],0,sizeof(double)*intercept_overhead.size());
  internal::bsp_counter=0;

  comp_kernel_info.clear();
  comm_kernel_info.clear();
  
  is_first_request=true;
}

void open_symbol(const char* symbol, double curtime){
  if (symbol_timers.find(symbol) == symbol_timers.end()){
    symbol_timers[symbol] = symbol_tracker(symbol);
    if (symbol_timers.size() < max_num_symbols) symbol_order[symbol_timers.size()-1] = std::string(symbol);
    symbol_timers[symbol].start(curtime);
  }
  else{
    symbol_timers[symbol].start(curtime);
  }
}

void close_symbol(const char* symbol, double curtime){
  if (symbol_timers.find(symbol) == symbol_timers.end()){ assert(0); }
  else{
    symbol_timers[symbol].stop(curtime);
  }
}


void final_accumulate(MPI_Comm comm, double last_time){
  assert(nonblocking_internal_info.size() == 0);
  critical_path_costs[num_critical_path_measures-3]+=(last_time-computation_timer);	// update critical path computation time
  critical_path_costs[num_critical_path_measures-1]+=(last_time-computation_timer);	// update critical path runtime
  volume_costs[num_volume_measures-3]+=(last_time-computation_timer);			// update computation time volume
  volume_costs[num_volume_measures-1]+=(last_time-computation_timer);			// update runtime volume
  // update the computation time (i.e. time between last MPI synchronization point or comp kernel and this function invocation) along all paths decomposed by MPI communication routine
  // TODO: Why don't I apply the same update to max_per_process_costs?
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-comm_path_select_size-1-i] += (last_time-computation_timer); }
}


void clear(){
  symbol_timers.clear();
}

void finalize(){
  if (std::getenv("CRITTER_VIZ_FILE") != NULL){
    if (is_world_root){
      stream.close();
    }
  }
}

}
}
}
