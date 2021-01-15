#include <limits.h>

#include "util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace decomposition{

std::ofstream stream;
int track_synchronization,cost_model,path_decomposition;
int include_barrier_time;
int internal_tag,internal_tag1,internal_tag2;
int internal_tag3,internal_tag4,internal_tag5;
bool is_first_request,is_first_iter;
size_t decomp_text_width,text_width,path_count,path_measure_count;
char barrier_pad_send,barrier_pad_recv;
double comp_start_time;
// CommCost, SynchCost, CompCost,
// CommTime, SynchTime, CompTime, CompKernelTime, ExecTime
size_t num_cp_measures;
// CommCost, SynchCost, CompCost,
// IdleTime, CommTime, SynchTime, CompTime, CompKernelTime, ExecTime
size_t num_pp_measures;
// CommCost, SynchCost, CompCost,
// IdleTime, CommTime, SynchTime, CompTime, CompKernelTime, ExecTime
size_t num_vol_measures;
size_t num_decomp_cp_measures;
size_t num_decomp_pp_measures;
size_t num_decomp_vol_measures;
size_t cp_costs_size,pp_costs_size,vol_costs_size;
std::string path_select,path_measure_select;
std::vector<bool> path_decisions;
std::vector<float> intercept_overhead,global_intercept_overhead;
std::vector<float> cp_costs,max_pp_costs,vol_costs;
std::vector<char> synch_pad_send,synch_pad_recv,eager_pad;
std::vector<int> path_index,path_measure_index;
std::vector<float> cp_costs_foreign;
std::stack<std::string> symbol_stack;
std::vector<std::string> symbol_order;
std::map<comp_kernel_key,std::pair<int,float>> comp_kernel_info;
std::map<comm_kernel_key,std::pair<int,float>> comm_kernel_info;

void allocate(MPI_Comm comm){
  int _world_size; MPI_Comm_size(MPI_COMM_WORLD,&_world_size);
  int _world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&_world_rank);
  decomp_text_width = 15; text_width = 15; is_first_iter = true;
  internal_tag = 31133; internal_tag1 = internal_tag+1;
  internal_tag2 = internal_tag+2; internal_tag3 = internal_tag+3;
  internal_tag4 = internal_tag+4; internal_tag5 = internal_tag+5;
  intercept_overhead.resize(3,0);
  global_intercept_overhead.resize(3,0);

  if (std::getenv("CRITTER_TRACK_SYNCHRONIZATION") != NULL){
    track_synchronization = atoi(std::getenv("CRITTER_TRACK_SYNCHRONIZATION"));
    assert(track_synchronization >= 0 && track_synchronization <= 1);
  } else{ track_synchronization = 0; }
  if (std::getenv("CRITTER_INCLUDE_BARRIER_TIME") != NULL){
    include_barrier_time = atoi(std::getenv("CRITTER_INCLUDE_BARRIER_TIME"));
    assert(include_barrier_time >= 0 && include_barrier_time <= 1);
  } else{ include_barrier_time = 0; }
  if (std::getenv("CRITTER_COST_MODEL") != NULL){
    cost_model = atoi(std::getenv("CRITTER_COST_MODEL"));
    assert(cost_model >= 0 && cost_model <= 1);
  } else cost_model = 0;// BSP
  if (std::getenv("CRITTER_PATH_DECOMPOSITION") != NULL){
    path_decomposition = atoi(std::getenv("CRITTER_PATH_DECOMPOSITION"));
    assert(path_decomposition >=0 && path_decomposition <= 2);
  } else path_decomposition = 0;// Neither symbols nor MPI path breakdown
  if (std::getenv("CRITTER_PATH_SELECT") != NULL){
    path_select = std::getenv("CRITTER_PATH_SELECT");
    assert(path_select.size() == 6);
  } else path_select = "000000";
  // SynchCost,CommCost,CompCost,CommTime,CompTime,ExecTime
  if (std::getenv("CRITTER_PATH_MEASURE_SELECT") != NULL){
    path_measure_select = std::getenv("CRITTER_PATH_MEASURE_SELECT");
    if (path_decomposition<=1) assert(path_measure_select.size() == 4);
    else if (path_decomposition==2) assert(path_measure_select.size() == 7);
  } else{
    if (path_decomposition == 2){
      // SynchCost,CommCost,CompCost,CommTime,SynchTime,CompTime,ExecTime
      path_measure_select = "0000000";
    } else if (path_decomposition == 1){
      // SynchCost,CommCost,CommTime,SynchTime
      path_measure_select = "0000";
    } else{
      path_measure_select = "0000";
    }
  }
  if (std::getenv("CRITTER_VIZ_FILE") != NULL){
    std::string stream_name = std::getenv("CRITTER_VIZ_FILE");
    stream_name += "_decomposition.txt";
    if (is_world_root){
      stream.open(stream_name.c_str(),std::ofstream::out);
      assert(stream.good());
    }
  }

  path_count=0;
  for (auto i=0; i<path_select.size(); i++){
    if (path_select[i] == '1'){ path_count++;
                                path_index.push_back(i); }
  } 
  path_measure_count=0;
  for (auto i=0; i<path_measure_select.size(); i++){
    if (path_measure_select[i] == '1'){ path_measure_count++;
                                        path_measure_index.push_back(i); }
  } 
  path_decisions.resize(path_count);
  if (track_synchronization){ synch_pad_send.resize(_world_size);
                              synch_pad_recv.resize(_world_size); }

  // 1st term - number of timers
  // 2nd term - number of costs
  // Each of these measures a single critical path, with no decomposition
  num_cp_measures = 5+3;
  num_pp_measures = 6+3;
  num_vol_measures = 6+3;
  // These are the measures we attribute to each tracker/symbol
  num_decomp_cp_measures = path_measure_count;
  num_decomp_pp_measures = path_measure_select.size();
  num_decomp_vol_measures = path_measure_select.size();

  // The '4*path_count' used below measure along each tracked path:
  // {computation cost, idle time, computation time, computational kernel time}
  cp_costs_size = num_cp_measures;
  pp_costs_size = num_pp_measures;
  vol_costs_size = num_vol_measures;
  if (path_decomposition <= 1){
    cp_costs_size += num_decomp_cp_measures*path_count*list_size;
    cp_costs_size += 4*path_count;
    pp_costs_size += num_decomp_pp_measures*path_count*list_size;
    pp_costs_size += 4*path_count;
    vol_costs_size += num_decomp_vol_measures*list_size;
    cp_costs.resize(cp_costs_size);
    cp_costs_foreign.resize(cp_costs_size);
    max_pp_costs.resize(pp_costs_size);
    vol_costs.resize(vol_costs_size);

    int eager_msg_sizes[2];
    MPI_Pack_size(1,MPI_CHAR,comm,&eager_msg_sizes[0]);
    MPI_Pack_size(cp_costs_size,MPI_FLOAT,comm,&eager_msg_sizes[1]);
    int eager_pad_size = 2*MPI_BSEND_OVERHEAD;
    for (int i=0; i<2; i++) { eager_pad_size += eager_msg_sizes[i]; }
    assert(eager_pad_size <= eager_limit);
    eager_pad.resize(eager_pad_size);
  }
}

void reset(){
  if (std::getenv("CRITTER_MODE") != NULL){
    mode = atoi(std::getenv("CRITTER_MODE"));
    assert(mode >=0 && mode <=1);
  } else mode = 1;

  for (auto i=0; i<list_size; i++){ list[i]->init(); }
  memset(&cp_costs[0],0,sizeof(float)*cp_costs.size());
  memset(&cp_costs_foreign[0],0,sizeof(float)*cp_costs.size());
  memset(&max_pp_costs[0],0,sizeof(float)*max_pp_costs.size());
  memset(&vol_costs[0],0,sizeof(float)*vol_costs.size());
  memset(&intercept_overhead[0],0,sizeof(float)*intercept_overhead.size());
  bsp_counter=0; is_first_request=true;
  comp_kernel_info.clear(); comm_kernel_info.clear();
}

void init_symbol(std::vector<std::string>& symbols){
  if (path_decomposition != 2) return;
  // We only want to initialize once, but do not know the number of symbols before this point.
  size_t cp_symbol_select_size = 4*num_decomp_cp_measures+1;
  size_t pp_symbol_select_size = 4*num_decomp_pp_measures+1;
  size_t vol_symbol_select_size = 4*num_decomp_vol_measures+1;
  cp_costs_size = num_cp_measures + cp_symbol_select_size*path_count*symbols.size();
  pp_costs_size = num_pp_measures + pp_symbol_select_size*symbols.size();
  vol_costs_size = num_vol_measures + vol_symbol_select_size*symbols.size();
  cp_costs.resize(cp_costs_size);
  cp_costs_foreign.resize(cp_costs_size);
  max_pp_costs.resize(pp_costs_size);
  vol_costs.resize(vol_costs_size);
  if (eager_pad.size()==0){
    int eager_msg_sizes[2];
    MPI_Pack_size(1,MPI_CHAR,MPI_COMM_WORLD,&eager_msg_sizes[0]);
    MPI_Pack_size(cp_costs_size,MPI_FLOAT,MPI_COMM_WORLD,&eager_msg_sizes[1]);
    int eager_pad_size = 2*MPI_BSEND_OVERHEAD;
    for (int i=0; i<2; i++) { eager_pad_size += eager_msg_sizes[i]; }
    assert(eager_pad_size <= eager_limit);
    eager_pad.resize(eager_pad_size);
  }

  for (auto& it :symbols){
    if (symbol_timers.find(it) == symbol_timers.end()){
      symbol_timers[it] = symbol_tracker(it);
      symbol_order.push_back(it);
      decomp_text_width = std::max(decomp_text_width,it.size());
    }
  }
}

void open_symbol(const char* symbol, double curtime){
  if (path_decomposition != 2) return;
  // If symbol not found already, ignore the symbol.
  // The new instructions are for user to specify tracked symbols apriori
  if (symbol_timers.find(symbol) != symbol_timers.end()){
    symbol_timers[symbol].start(curtime); }
}

void close_symbol(const char* symbol, double curtime){
  if (path_decomposition != 2) return;
  // If symbol not found already, ignore the symbol.
  // The new instructions are for user to specify tracked symbols apriori
  if (symbol_timers.find(symbol) != symbol_timers.end()){
    symbol_timers[symbol].stop(curtime); }
}

void final_accumulate(MPI_Comm comm, double last_time){
  assert(nonblocking_internal_info.size() == 0);
  auto last_comp_time = last_time-computation_timer;
  cp_costs[num_cp_measures-3] += last_comp_time;
  cp_costs[num_cp_measures-1] += last_comp_time;
  vol_costs[num_vol_measures-3] += last_comp_time;
  vol_costs[num_vol_measures-1] += last_comp_time;
  vol_costs[num_vol_measures-3] = std::min(vol_costs[num_vol_measures-3],cp_costs[num_cp_measures-3]);
  vol_costs[num_vol_measures-1] = std::min(vol_costs[num_vol_measures-1],cp_costs[num_cp_measures-1]);
  if (path_decomposition==1){
    for (size_t i=0; i<path_count; i++){
      cp_costs[cp_costs_size-path_count-1-i] += last_comp_time; }
  }
}

void clear(){ symbol_timers.clear(); }

void finalize(){
  if (std::getenv("CRITTER_VIZ_FILE") != NULL){
    if (is_world_root){ stream.close(); }
  }
}

}
}
}
