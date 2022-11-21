#include "util.h"
#include "container/comm_tracker.h"
#include "container/kernel_tracker.h"

namespace internal{

int save_wildcard_id;
int bsp_counter;
int propagate_within_timer;
int world_rank,debug_rank,debug;
volatile double computation_timer;
std::vector<double> wall_timer;
double _wall_time;
size_t auto_capture;
bool is_world_root;
size_t mode,stack_id;
std::vector<float> scratch_pad;
size_t profile_blas1;
size_t profile_blas2;
size_t profile_blas3;
size_t profile_lapack;
size_t profile_collective;
size_t profile_p2p;
size_t propagate_collective;
size_t propagate_p2p;
size_t execute_kernels;
size_t execute_kernels_max_message_size;
size_t eager_limit;
int request_id;

std::ofstream stream;
int cost_model,path_decomposition;
int internal_tag,internal_tag1,internal_tag2;
int internal_tag3,internal_tag4,internal_tag5;
bool is_first_request,is_first_iter;
size_t decomp_text_width,text_width,path_count,path_measure_count;
double comp_start_time;
// CommCost, SynchCost, CompCost,
// CommTime, ExecTime
size_t num_cp_measures;
// CommCost, SynchCost, CompCost,
// CommTime, ExecTime
size_t num_pp_measures;
// CommCost, SynchCost, CompCost,
// CommTime, ExecTime
size_t num_vol_measures;
size_t num_decomp_cp_measures;
size_t num_decomp_pp_measures;
size_t num_decomp_vol_measures;
size_t cp_costs_size,pp_costs_size,vol_costs_size;
size_t num_kernel_ds,exclusive_only,max_num_tracked_kernels;
std::string path_select,path_measure_select;
std::vector<bool> path_decisions;
std::vector<float> cp_costs,max_pp_costs,vol_costs;
std::vector<float> cp_costs_foreign;
std::vector<char> eager_pad;
std::vector<int> path_index,path_measure_index;
std::stack<std::string> timer_stack;
std::stack<bool> propagate_within_timer_stack;

void initialize(MPI_Comm comm){
  // Allocate is always invoked when MPI_Init* is intercepted.
  int _world_size; MPI_Comm_size(MPI_COMM_WORLD,&_world_size);
  int _world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&_world_rank);
  // Set up print-statement details. These should ideally be user-specified.
  decomp_text_width = 40; text_width = 15; is_first_iter = true;
  // Select arbitrary tag for internal point-to-point
  //   messages that will (hopefully) not conflict with user messages
  internal_tag = 31133; internal_tag1 = internal_tag+1;
  internal_tag2 = internal_tag+2; internal_tag3 = internal_tag+3;
  internal_tag4 = internal_tag+4; internal_tag5 = internal_tag+5;

  // Two cost-model options: alpha-beta (1) and BSP (0)
  if (std::getenv("CRITTER_COST_MODEL") != NULL){
    cost_model = atoi(std::getenv("CRITTER_COST_MODEL"));
    assert(cost_model >= 0 && cost_model <= 1);
  } else cost_model = 0;// BSP
  // Three path decompositions: total measures (0),
  //                            communication routine breakdown along specific critical paths (1),
  //                            kernel breakdown along specific critical paths (2).
  if (std::getenv("CRITTER_PATH_PROFILE") != NULL){
    path_decomposition = atoi(std::getenv("CRITTER_PATH_PROFILE"));
    assert(path_decomposition >=0 && path_decomposition <= 2);
  } else path_decomposition = 0;// Neither symbols nor MPI path breakdown
  // Select which critical-path to decompose (in bit-string order):
  //   {SynchCost,CommCost,CompCost,CommTime,ExecTime}
  if (std::getenv("CRITTER_PATH_SELECT") != NULL){
    path_select = std::getenv("CRITTER_PATH_SELECT");
    assert(path_select.size() == 5);
  } else path_select = "00000";// Decompose down none of them (lowest overhead)
  // Select which measurements to track along the specified (see above) critical paths.
  //   If path_decomposition==1, choose among (in bit-string order): {CommCost, SynchCost, CommTime}.
  //   If path_decomposition==2, choose among (in bit-string order): {CommCost, SynchCost, CompCost, CommTime, ExecTime}.
  if (std::getenv("CRITTER_PATH_MEASURE_SELECT") != NULL){
    path_measure_select = std::getenv("CRITTER_PATH_MEASURE_SELECT");
    // The reason we currently force specification of 3 path measures for a particular path
    //   when setting path_decomposition==1 is because the execution time and computational cost
    //   along the particular path can be tracked usng just 2 doubles, whereas the other 3 measures
    //   are tracked per MPI routine (and thus take up more memory and incur more overhead).
    // I will now change that so that a bit-string of size 5 is always expected.
    // TODO: There is one limitation: if we want to track communication time, but not specific
    //   to a particular MPI routine, we cannot support this without tracking a distinct
    //   bucket to store the contribution from each MPI routine and add it up at the end.
    assert(path_measure_select.size() == 5);
  } else{// {SynchCost,CommCost,CompCost,CommTime,ExecTime}
      path_measure_select = "00000";
  }
  // Specify file-name prefix
  if (std::getenv("CRITTER_VIZ_FILE") != NULL){
    std::string stream_name = std::getenv("CRITTER_VIZ_FILE");
    stream_name += "_decomposition.txt";
    if (is_world_root){
      stream.open(stream_name.c_str(),std::ofstream::out);
      assert(stream.good());
    }
  }
  // Specify whether to profile kernel exclusive time only or
  //   both exclusive time and inclusive time.
  // The benefit of the former is less communication overhead.
  if (std::getenv("CRITTER_PROFILE_EXCLUSIVE_TIME_ONLY") != NULL){
    exclusive_only = atoi(std::getenv("CRITTER_PROFILE_EXCLUSIVE_TIME_ONLY"));
    assert(exclusive_only >=0 && exclusive_only <= 2);
  } else exclusive_only = 0;
  if (exclusive_only == 1){
    num_kernel_ds = 1;
  } else{
    num_kernel_ds = 4;
  }
  // Specify whether to profile kernel exclusive time only or
  //   both exclusive time and inclusive time.
  // The benefit of the former is less communication overhead.
  if (std::getenv("CRITTER_PROFILE_MAX_NUM_KERNELS") != NULL){
    max_num_tracked_kernels = atoi(std::getenv("CRITTER_PROFILE_MAX_NUM_KERNELS"));
  } else max_num_tracked_kernels = 10;

  path_count=0;
  for (auto i=0; i<path_select.size(); i++){
    if (path_select[i] == '1'){ path_count++;
                                path_index.push_back(i); }
  } 
  path_measure_count=0;
  for (auto i=0; i<path_measure_select.size(); i++){
    if (path_measure_select[i] == '1'){ if ((path_decomposition==1 && i<3) || (path_decomposition==2)){// Extra logic not needed otherwise.
                                         path_measure_count++;
                                         path_measure_index.push_back(i); } }
  } 
  path_decisions.resize(path_count);// Used in internal reductions of path profiles
  scratch_pad.resize(100);	// necessary (100 perhaps overkill) to avoid overwriting past bounds of global variable in kernel_tracker

  // 1st term - number of timers
  // 2nd term - number of costs
  // Each of these measures a single critical path, with no decomposition
  num_cp_measures = 2+3;
  num_pp_measures = 2+3;
  num_vol_measures = 2+3;
  // These are the measures we attribute to each tracker/symbol
  // Note that, below, these aren't set to the same value
  //   because there is negligable overhead to tracking
  //   more max-per-process or volumetric measures, but there
  //   certainly is more overhead to tracking more critical-path
  //   measures.
  num_decomp_cp_measures = path_measure_count;
  num_decomp_pp_measures = path_measure_select.size();
  num_decomp_vol_measures = path_measure_select.size();
  if (path_decomposition <= 1){
    num_decomp_pp_measures -= 2;// Only need to track {CommCost,SynchCost,CommTime} per MPI routine
    num_decomp_vol_measures -= 2;// Only need to track {CommCost,SynchCost,CommTime} per MPI routine
  }

  // The '2*path_count' used below measure along each tracked path:
  // {computation cost, execution time}
  cp_costs_size = num_cp_measures;
  pp_costs_size = num_pp_measures;
  vol_costs_size = num_vol_measures;
  if (path_decomposition <= 1){
    cp_costs_size += num_decomp_cp_measures*path_count*list_size;
    cp_costs_size += 2*path_count;
    pp_costs_size += num_decomp_pp_measures*path_count*list_size;
    pp_costs_size += 2*path_count;
    // Note: cp_costs and pp_costs are stored differently.
    // cp_costs puts the per-path measures at the end in groups,
    //   whereas pp_costs puts the per-path measures at the end of a
    //   collection of contributions from all MPI routines per path.
    vol_costs_size += num_decomp_vol_measures*list_size;
    cp_costs.resize(cp_costs_size);
    cp_costs_foreign.resize(cp_costs_size);
    max_pp_costs.resize(pp_costs_size);
    vol_costs.resize(vol_costs_size);
  } else if (path_decomposition == 2){
    // We only want to initialize once, but do not know the number of symbols before this point.
    size_t cp_symbol_select_size = num_kernel_ds*num_decomp_cp_measures+1;
    size_t pp_symbol_select_size = num_kernel_ds*num_decomp_pp_measures+1;
    size_t vol_symbol_select_size = num_kernel_ds*num_decomp_vol_measures+1;
    cp_costs_size = num_cp_measures + cp_symbol_select_size*path_count*max_num_tracked_kernels;
    pp_costs_size = num_pp_measures + pp_symbol_select_size*max_num_tracked_kernels;
    vol_costs_size = num_vol_measures + vol_symbol_select_size*max_num_tracked_kernels;
    cp_costs.resize(cp_costs_size);
    cp_costs_foreign.resize(cp_costs_size);
    max_pp_costs.resize(pp_costs_size);
    vol_costs.resize(vol_costs_size);
  }
  if (eager_pad.size()==0){
    int eager_msg_sizes[2];
    MPI_Pack_size(1,MPI_CHAR,MPI_COMM_WORLD,&eager_msg_sizes[0]);
    MPI_Pack_size(cp_costs_size,MPI_FLOAT,MPI_COMM_WORLD,&eager_msg_sizes[1]);
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
  bsp_counter=0; is_first_request=true;
}

void register_timer(const char* timer_name){
  assert(path_decomposition == 2);
  if (timers.find(timer_name) == timers.end()){
    timers[timer_name] = kernel_tracker(timer_name);
  }
}

void __start_timer__(const char* timer_name, double curtime, bool propagate_within){
  // If timer_name not found already, ignore the timer.
  // The new instructions are for user to specify tracked kernels apriori
  if (timers.find(timer_name) == timers.end() && timers.size() < max_num_tracked_kernels){
    timers[timer_name] = kernel_tracker(timer_name);
  }
  if (timers.find(timer_name) != timers.end()) {
    //NOTE: If propagate_within_timer=1, critical path information will not be propagated
    //      If propagate_within_timer==1, this cannot be overwritten by timers invoked within the call stack,
    //        so state must be saved internally.
    propagate_within_timer = propagate_within && (propagate_within_timer_stack.size()==0 ? true : propagate_within_timer_stack.top());
    propagate_within_timer_stack.push(propagate_within_timer);
    timers[timer_name].start(curtime);
  }
}

void __stop_timer__(const char* timer_name, double curtime){
  // If timer_name not found already, ignore the timer.
  // The new instructions are for user to specify tracked timers apriori
  if (timers.find(timer_name) != timers.end()){
    timers[timer_name].stop(curtime);
    if (propagate_within_timer_stack.size()>0) propagate_within_timer_stack.pop();
    propagate_within_timer = (propagate_within_timer_stack.size()==0 ? true : propagate_within_timer_stack.top());
  }
}

void update_time(double last_time){
  auto last_comp_time = last_time-computation_timer;
  // Update ExecTime.
  cp_costs[num_cp_measures-1] += last_comp_time;
  vol_costs[num_vol_measures-1] += last_comp_time;
  // Due to granularity of timing measurements, force time along critical-path to be
  //   greater than max per-process time (which doesn't include idle time) if necessary.
  //vol_costs[num_vol_measures-1] = std::min(vol_costs[num_vol_measures-1],cp_costs[num_cp_measures-1]);
  if (path_decomposition==1){
    for (size_t i=0; i<path_count; i++){
      cp_costs[cp_costs_size-1-i] += last_comp_time; }
  }
}

void collect_volumetric_statistics(MPI_Comm cm){
  if (mode==0) return;
  assert(nonblocking_internal_info.size() == 0);
  // First compute per-process max
  int rank; MPI_Comm_rank(cm,&rank);
  float_int buffer[num_pp_measures];
  for (size_t i=0; i<num_pp_measures; i++){
    max_pp_costs[i] = vol_costs[i];
    buffer[i].first = vol_costs[i];
    buffer[i].second = rank;
  }
  PMPI_Allreduce(MPI_IN_PLACE, &max_pp_costs[0], num_pp_measures, MPI_FLOAT, MPI_MAX, cm);
  PMPI_Allreduce(MPI_IN_PLACE, &buffer[0], num_pp_measures, MPI_FLOAT_INT, MPI_MAXLOC, cm);
  if (path_decomposition == 2){
    // First 'num_pp_measures' will be overwritten by the same data
    std::memcpy(&vol_costs[num_vol_measures], &max_pp_costs[num_pp_measures], (max_pp_costs.size()-num_pp_measures)*sizeof(float));
  }
  if (path_decomposition == 1 && path_count>0){
    size_t save=0;
    for (size_t i=0; i<path_index.size(); i++){
      size_t z = path_index[i];
      // If I am the processor that incurs the largest per-process time along any one path, then
      //   save the data and communicate it via reduction.
      if (rank == buffer[z].second){
        for (size_t j=0; j<num_decomp_pp_measures*list_size; j++){
          // The magic number 2 is for tracking the CompCost and ExecutionTime along each path
          max_pp_costs[num_pp_measures+save*(num_decomp_pp_measures*list_size+2)+j] = vol_costs[num_vol_measures+j];
        }
        // Set the two final measurements
        max_pp_costs[num_pp_measures+(save+1)*(num_decomp_pp_measures*list_size+2)-2] = vol_costs[num_vol_measures-3];// CompCost
        max_pp_costs[num_pp_measures+(save+1)*(num_decomp_pp_measures*list_size+2)-1] = vol_costs[num_vol_measures-1];// ExecTime
      }
      else{
        for (size_t j=0; j<num_decomp_pp_measures*list_size+2; j++){
          max_pp_costs[num_pp_measures+save*(num_decomp_pp_measures*list_size+2)+j] = 0.;
        }
      }
      PMPI_Allreduce(MPI_IN_PLACE, &max_pp_costs[num_pp_measures+save*(num_decomp_pp_measures*list_size+2)], num_decomp_pp_measures*list_size+2, MPI_FLOAT, MPI_MAX, cm);
      save++;
    }
  }
  else if (path_decomposition == 2 && path_count>0){
    size_t path_select_offset = num_kernel_ds*num_decomp_pp_measures+1;
    for (size_t i=0; i<path_index.size(); i++){
      size_t z = path_index[i];
      if (rank == buffer[z].second){
        for (size_t j=0; j<timers.size(); j++){
          std::memcpy(&max_pp_costs[num_pp_measures+j*path_select_offset*path_count+i*path_select_offset],
                      &vol_costs[num_pp_measures+j*path_select_offset*path_count+i*path_select_offset],
                      path_select_offset*sizeof(float));
        }
      }
      else{
        for (size_t j=0; j<timers.size(); j++){
          std::memset(&max_pp_costs[num_pp_measures+j*path_select_offset*path_count+i*path_select_offset],
                      (float)0., path_select_offset*sizeof(float));
        }
      }
    }
    PMPI_Allreduce(MPI_IN_PLACE, &max_pp_costs[0], (max_pp_costs.size()-num_pp_measures), MPI_FLOAT, MPI_MAX, cm);
  }
  // Now compute volumetric average
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  PMPI_Allreduce(MPI_IN_PLACE, &vol_costs[0], vol_costs.size(), MPI_FLOAT, MPI_SUM, cm);
  for (int i=0; i<vol_costs.size(); i++){ vol_costs[i] /= world_size; }
/*
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  if (mode && symbol_path_select_size>0){
    size_t active_size = world_size;
    size_t active_rank = world_rank;
    size_t active_mult = 1;
    while (active_size>1){
      if (active_rank % 2 == 1){
        int partner = (active_rank-1)*active_mult;
        PMPI_Send(...)
        break;
      }
      else if ((active_rank % 2 == 0) && (active_rank < (active_size-1))){
        int partner = (active_rank+1)*active_mult;
        PMPI_Recv(...)
      }
      active_size = active_size/2 + active_size%2;
      active_rank /= 2;
      active_mult *= 2;
    }
  }
*/
}

void finalize(){
  if (std::getenv("CRITTER_VIZ_FILE") != NULL){
    if (is_world_root){ stream.close(); }
  }
}

}
