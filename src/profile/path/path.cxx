#include "path.h"
#include "../container/symbol_tracker.h"
#include "../../accelerate/util/util.h"
#include "../util/util.h"
#include "../../util/util.h"

namespace critter{
namespace internal{
namespace profile{

void path::exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm){
  assert(0); return;
  // Accumulate the computation time between last communication routine
  //   as both execution-time and computation time into both
  //   the execution-time critical path data structures and the
  //   per-process data structures.
  volatile auto curtime = MPI_Wtime();
  auto save_comp_time = curtime - computation_timer;
  if (mode==1){
    cp_costs[num_cp_measures-1] += save_comp_time;
    vol_costs[num_vol_measures-1] += save_comp_time;
    if (path_decomposition == 1){
      for (size_t i=0; i<path_count; i++){
        cp_costs[cp_costs_size-path_count-1-i] += save_comp_time; }
    }
    if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
      // Get the current symbol's execution-time since last
      //   communication routine or its inception.
      // Accumulate as both execution-time and computation time
      //   into both the execution-time critical path data structures
      //     and the per-process data structures.
      auto last_symbol_time = curtime - symbol_timers[symbol_stack.top()].start_timer.top();
      for (auto i=0; i<path_count; i++){
        for (auto j=0; j<path_measure_index.size(); j++){
          if (path_measure_index[j] == path_measure_select.size()-1){
            symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += last_symbol_time;
            symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += last_symbol_time;
          }
        }
      }
      symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_decomp_pp_measures-1] += last_symbol_time;
      symbol_timers[symbol_stack.top()].pp_excl_measure[num_decomp_pp_measures-1] += last_symbol_time;
    }
  }
  if (mode==1){
    computation_timer = MPI_Wtime();
    if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
      symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
  }
}

bool path::initiate_comp(size_t id, volatile double curtime, float flop_count,
                         int param1, int param2, int param3, int param4, int param5){
  auto save_comp_time = curtime - computation_timer;
  // Update {CompTime, ExecTime}
  cp_costs[num_cp_measures-1] += save_comp_time;
  vol_costs[num_vol_measures-1] += save_comp_time;
  if (path_decomposition == 1){
    for (size_t i=0; i<path_count; i++){
      cp_costs[cp_costs_size-1-i] += save_comp_time; }
  }
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    // Get the current symbol's execution-time since last
    //   communication routine or its inception.
    // Accumulate as both execution-time and computation time
    //   into both the execution-time critical path data structures
    //     and the per-process data structures.
    auto last_symbol_time = curtime - symbol_timers[symbol_stack.top()].start_timer.top();
    for (auto i=0; i<path_count; i++){
      for (auto j=0; j<path_measure_index.size(); j++){
        if (path_measure_index[j] == path_measure_select.size()-1){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += last_symbol_time;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += last_symbol_time;
        }
      }
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_decomp_pp_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_decomp_pp_measures-1] += last_symbol_time;
  }

  // start compunication timer for compunication routine
  comp_start_time = MPI_Wtime();
  return true;
}

void path::complete_comp(double errtime, size_t id, float flop_count,
                         int param1, int param2, int param3, int param4, int param5){
  volatile auto comp_time = MPI_Wtime() - comp_start_time - errtime;	// complete computation time
/*
  // Save kernel information
  if (autotuning_debug == 1){
    comp_kernel_key key(-1,id,flop_count,param1,param2,param3,param4,param5);
    if (comp_kernel_info.find(key) == comp_kernel_info.end()){
      comp_kernel_info[key] = std::make_pair(1,comp_time);
    } else{
      comp_kernel_info[key].first++;
      comp_kernel_info[key].second += comp_time;
    }
  }
*/
  // Decompose measurements along multiple paths by symbol
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    for (auto i=0; i<path_count; i++){
      for (auto j=0; j<path_measure_index.size(); j++){
        if (path_measure_index[j] == path_measure_select.size()-3){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += flop_count;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += flop_count;
        } else if (path_measure_index[j] == path_measure_select.size()-1){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += comp_time;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += comp_time;
        }
      }
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_decomp_pp_measures-3] += flop_count;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_decomp_pp_measures-3] += flop_count;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_decomp_pp_measures-1] += comp_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_decomp_pp_measures-1] += comp_time;
  }
  cp_costs[num_cp_measures-3] += flop_count;
  cp_costs[num_cp_measures-1] += comp_time;
  if (path_decomposition == 1){
    for (size_t i=0; i<path_count; i++){ cp_costs[cp_costs_size-1-i] += comp_time; }
    for (size_t i=0; i<path_count; i++){ cp_costs[cp_costs_size-path_count-1-i] += flop_count; }
  }
  vol_costs[num_vol_measures-3] += flop_count;
  vol_costs[num_vol_measures-1] += comp_time;

  computation_timer = MPI_Wtime();
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
}

static void update_cp_decomp1(float* in, float* inout, size_t len){
  assert(len == cp_costs_size);	// this assert prevents user from obtaining wrong output if MPI implementation cuts up the message.
  // We do not track paths for synchronization time nor comp-kernel time, so check for those
  for (int i=0; i<path_index.size(); i++){
    int j=path_index[i];
    path_decisions[i] = inout[j] > in[j];
  }
  for (int i=0; i<num_cp_measures; i++){
    inout[i] = std::max(inout[i],in[i]);
  }
  for (int i=num_cp_measures; i<cp_costs_size; i++){
    int idx = (i-num_cp_measures)%path_count;
    inout[i] = (path_decisions[idx] ? inout[i] : in[i]);
  }
}
static void update_cp_decomp2(float* in, float* inout, size_t len){
  assert(len == cp_costs_size);	// this assert prevents user from obtaining wrong output if MPI implementation cuts up the message.
  assert(path_count>0);
  // We do not track paths for synchronization time nor comp-kernel time, so check for those
  for (int i=0; i<path_index.size(); i++){
    int j=path_index[i];
    path_decisions[i] = inout[j] > in[j];
  }
  for (int i=0; i<num_cp_measures; i++){
    inout[i] = std::max(inout[i],in[i]);
  }
  size_t path_select_offset = 4*num_decomp_cp_measures+1;
  for (int i=num_cp_measures; i<cp_costs_size; i++){
    int idx = (i-num_cp_measures)%(path_count*path_select_offset);// restarts for each symbol
    idx /= path_select_offset;
    inout[i] = (path_decisions[idx] ? inout[i] : in[i]);
  }
}


static void propagate_cp_decomp1_op(float* in, float* inout, int* len, MPI_Datatype* dtype){
  update_cp_decomp1(in,inout,static_cast<size_t>(*len));
}
static void propagate_cp_decomp2_op(float* in, float* inout, int* len, MPI_Datatype* dtype){
  update_cp_decomp2(in,inout,static_cast<size_t>(*len));
}

bool path::initiate_comm(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                            bool is_sender, int partner1, int user_tag1, int partner2, int user_tag2){
  // Check for conflicting communication tag(s)
  if (partner1 != -1) assert(!(user_tag1 >= internal_tag && user_tag1 <= internal_tag5));
  if (partner2 != -1) assert(!(user_tag2 >= internal_tag && user_tag2 <= internal_tag5));
  // Save and accumulate the computation time between last communication routine as both execution-time and computation time
  //   into both the execution-time critical path data structures and the per-process data structures.
  tracker.comp_time = curtime - computation_timer;
  int rank; MPI_Comm_rank(comm, &rank);

  // Save caller communication attributes into reference object for use in corresponding static method 'complete'
  int word_size,np; MPI_Type_size(t, &word_size);
  int64_t nbytes = word_size * nelem;
  MPI_Comm_size(comm, &np);
  tracker.nbytes = nbytes;
  tracker.comm = comm;
  tracker.comm_size = np;
  tracker.is_sender = is_sender;
  tracker.partner1 = partner1;
  tracker.partner2 = partner2 != -1 ? partner2 : partner1;// Useful in propagation
  if (partner1 == MPI_ANY_SOURCE){// only possible for 'MPI_Recv'
    MPI_Status st;
    PMPI_Probe(partner1,user_tag1,comm,&st);
    save_wildcard_id = st.MPI_SOURCE;
    tracker.partner1 = save_wildcard_id;
  } else if (partner2 == MPI_ANY_SOURCE){// only possible for 'MPI_Sendrecv' / 'MPI_Sendrecv_replace'
    assert((tracker.tag == 13) || (tracker.tag == 14));
    MPI_Status st;
    PMPI_Probe(partner2,user_tag2,comm,&st);
    save_wildcard_id = st.MPI_SOURCE;
    tracker.partner2 = save_wildcard_id;
  }

  cp_costs[num_cp_measures-1] += tracker.comp_time;
  vol_costs[num_vol_measures-1] += tracker.comp_time;
  if (path_decomposition == 1){
    for (size_t i=0; i<path_count; i++){
      cp_costs[cp_costs_size-1-i] += tracker.comp_time;
    }
  }
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    auto last_symbol_time = curtime - symbol_timers[symbol_stack.top()].start_timer.top();
    for (auto i=0; i<path_count; i++){
      for (auto j=0; j<path_measure_index.size(); j++){
        if (path_measure_index[j] == path_measure_select.size()-1){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += last_symbol_time;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += last_symbol_time;
        }
      }
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_decomp_pp_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_decomp_pp_measures-1] += last_symbol_time;
  }

  volatile auto init_time = MPI_Wtime();
  propagate(tracker);

  // Do as much work as possible post-propagate and pre-barrier (which should be "dead" time that shouldn't affect correctness).
  // Update measurements that define the critical path for each metric.
  std::pair<float,float> cost = cost_model == 0 ?
    tracker.cost_func_bsp(tracker.nbytes, tracker.comm_size)
    : tracker.cost_func_alphabeta(tracker.nbytes, tracker.comm_size);
  // Decompose measurements along multiple paths by MPI routine.
  *tracker.my_msg_count += cost.first;
  *tracker.my_wrd_count += cost.second;
  for (size_t i=0; i<path_count; i++){
    *(tracker.cp_msg_count+i) += cost.first;
    *(tracker.cp_wrd_count+i) += cost.second;
  }
  // Decompose measurements along multiple paths by symbol
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    for (auto i=0; i<path_count; i++){
      // update all communication-related measures for the top symbol in stack
      for (auto j=0; j<path_measure_index.size(); j++){
        if (path_measure_index[j] == 0){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += cost.second;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += cost.second;
        } else if (path_measure_index[j] == 1){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += cost.second;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += cost.second;
        }
      }
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[0] += cost.second;
    symbol_timers[symbol_stack.top()].pp_excl_measure[0] += cost.second;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[1] += cost.first;
    symbol_timers[symbol_stack.top()].pp_excl_measure[1] += cost.first;
  }

  cp_costs[num_cp_measures-5] += cost.second;
  cp_costs[num_cp_measures-4] += cost.first;

  vol_costs[num_vol_measures-5] += cost.second;
  vol_costs[num_vol_measures-4] += cost.first;
  void* temp_buf; int temp_size;
  MPI_Buffer_detach(&temp_buf,&temp_size);
  bsp_counter++;

  // To accurately measure communication time of a collective, must re-"synchronize" (blocking barrier, not synchronous).
  //   The need for this ensuing barrier has been validated by empirical studies.
  if (tracker.partner1 != -1) PMPI_Barrier(tracker.comm);
  // start communication timer for communication routine
  tracker.start_time = MPI_Wtime();
  return true;
}

void path::complete_comm(blocking& tracker){
  volatile auto comm_time = MPI_Wtime() - tracker.start_time;	// complete communication time
  int rank; MPI_Comm_rank(tracker.comm, &rank);

/*
  // Save kernel information
  if (autotuning_debug == 1){
    assert(comm_channel_map.find(tracker.comm) != comm_channel_map.end());
    int comm_sizes[2]={0,0}; int comm_strides[2]={0,0};
    for (auto i=0; i<comm_channel_map[tracker.comm]->id.size(); i++){
      comm_sizes[i]=comm_channel_map[tracker.comm]->id[i].first;
      comm_strides[i]=comm_channel_map[tracker.comm]->id[i].second;
    }
    comm_kernel_key key(rank,-1,tracker.tag,comm_sizes,comm_strides,tracker.nbytes,tracker.partner1);
    if (comm_kernel_info.find(key) == comm_kernel_info.end()){
      comm_kernel_info[key] = std::make_pair(1,comm_time);
    } else{
      comm_kernel_info[key].first++;
      comm_kernel_info[key].second += comm_time;
    }
  }
*/

  // Decompose measurements along multiple paths by MPI routine.
  for (size_t i=0; i<path_count; i++){
    *(tracker.cp_comm_time+i) += comm_time;
  }

  // Decompose measurements along multiple paths by symbol
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    for (auto i=0; i<path_count; i++){
      // update all communication-related measures for the top symbol in stack
      for (auto j=0; j<path_measure_index.size(); j++){
        if (path_measure_index[j] == path_measure_select.size()-2){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += comm_time;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += comm_time;
        } else if (path_measure_index[j] == path_measure_select.size()-1){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += comm_time;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += comm_time;
        }
      }
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_decomp_pp_measures-2] += comm_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_decomp_pp_measures-2] += comm_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_decomp_pp_measures-1] += comm_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_decomp_pp_measures-1] += comm_time;
  }

  // Update measurements that define the critical path for each metric.
  cp_costs[num_cp_measures-2] += comm_time;
  cp_costs[num_cp_measures-1] += comm_time;
  vol_costs[num_vol_measures-2] += comm_time;
  vol_costs[num_vol_measures-1] += comm_time;

  // Note that this block of code below is left in solely for blocking communication to avoid over-counting the idle time
  // If per-process execution-time gets larger than execution-time along the execution-time critical path,
  //   subtract out the difference from idle time.
  //auto path_diff = vol_costs[num_vol_measures-1]-cp_costs[num_cp_measures-1];
  //vol_costs[num_vol_measures-65] -= std::max((float)0.,path_diff);

/* Commented out for performance reasons. Often not necessary anyways.
  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
  for (int i=1; i<=5; i++){
    vol_costs[num_vol_measures-i] = std::min(vol_costs[num_vol_measures-i],cp_costs[num_cp_measures-i]);
  }
*/

  // Prepare to leave interception and re-enter user code by restarting computation timers.
  computation_timer = MPI_Wtime();
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
}

// Called by both nonblocking p2p and nonblocking collectives
bool path::initiate_comm(nonblocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm, bool is_sender, int partner, int user_tag){
  // Check for conflicting communication tag(s)
  if (partner != -1) assert(!(user_tag >= internal_tag && user_tag <= internal_tag5));
  tracker.comp_time = curtime - computation_timer;
  cp_costs[num_cp_measures-1] += tracker.comp_time;
  vol_costs[num_vol_measures-1] += tracker.comp_time;
  if (path_decomposition == 1){
    for (size_t i=0; i<path_count; i++){
      cp_costs[cp_costs_size-1-i] += tracker.comp_time;
    }
  }
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    auto last_symbol_time = curtime - symbol_timers[symbol_stack.top()].start_timer.top();
    for (auto i=0; i<path_count; i++){
      for (auto j=0; j<path_measure_index.size(); j++){
        if (path_measure_index[j] == path_measure_select.size()-1){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += last_symbol_time;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += last_symbol_time;
        }
      }
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_decomp_pp_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_decomp_pp_measures-1] += last_symbol_time;
  }
  // Note: routine below will be called immediately afterward.
  return true;
}

// Called by both nonblocking p2p and nonblocking collectives
void path::initiate_comm(nonblocking& tracker, volatile double itime, int64_t nelem,
                         MPI_Datatype t, MPI_Comm comm, MPI_Request* request, bool is_sender, int partner, int user_tag){
  // Check for conflicting communication tag(s)
  if (partner != -1) assert(!(user_tag >= internal_tag && user_tag <= internal_tag5));
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  tracker.comp_time = itime;
  cp_costs[num_cp_measures-1] += tracker.comp_time;
  vol_costs[num_vol_measures-1] += tracker.comp_time;
  if (path_decomposition == 1){
    for (size_t i=0; i<path_count; i++){
     cp_costs[cp_costs_size-1-i] += tracker.comp_time;
    }
  }
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    auto last_symbol_time = itime;
    for (auto i=0; i<path_count; i++){
      for (auto j=0; j<path_measure_index.size(); j++){
        if (path_measure_index[j] == path_measure_select.size()-1){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += last_symbol_time;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += last_symbol_time;
        }
      }
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_decomp_pp_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_decomp_pp_measures-1] += last_symbol_time;
  }

  int rank; MPI_Comm_rank(comm,&rank); 
  int el_size,p;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(comm, &p);
  MPI_Buffer_attach(&eager_pad[0],eager_pad.size());
  tracker.comm = comm;
  tracker.is_sender = is_sender;
  tracker.partner1 = partner;
  tracker.partner2 = -1;
  tracker.comm_size = p;
  tracker.nbytes = nbytes;

  MPI_Request prop_req = MPI_REQUEST_NULL;
  float* path_data = nullptr;
  if (propagate_p2p) propagate(tracker,path_data,&prop_req);
  nonblocking_info msg_info(path_data,prop_req,is_sender,partner,comm,(float)nbytes,(float)p,&tracker);
  nonblocking_internal_info[*request] = msg_info;

  void* temp_buf; int temp_size;
  MPI_Buffer_detach(&temp_buf,&temp_size);
  computation_timer = MPI_Wtime();
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
}

void path::complete_comm(double comp_time){
  cp_costs[num_cp_measures-1] += comp_time;
  if (path_decomposition == 1){
    for (size_t i=0; i<path_count; i++){
      cp_costs[cp_costs_size-1-i] += comp_time;
    }
  }
  vol_costs[num_vol_measures-1] += comp_time;
  // Decompose measurements along multiple paths by symbol
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    for (auto i=0; i<path_count; i++){
      for (auto j=0; j<path_measure_index.size(); j++){
        if (path_measure_index[j] == path_measure_select.size()-1){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += comp_time;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += comp_time;
        }
      }
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_decomp_pp_measures-1] += comp_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_decomp_pp_measures-1] += comp_time;
  }
}

void path::complete_comm(nonblocking& tracker, MPI_Request* request, double comp_time, double comm_time){
  auto info_it = nonblocking_internal_info.find(*request);
  assert(info_it != nonblocking_internal_info.end());

  tracker.is_sender = info_it->second.is_sender;
  tracker.comm = info_it->second.comm;
  tracker.partner1 = info_it->second.partner;
  tracker.partner2 = -1;
  tracker.nbytes = info_it->second.nbytes;
  tracker.comm_size = info_it->second.comm_size;

  int rank; MPI_Comm_rank(tracker.comm, &rank);
/*
  if (autotuning_debug == 1){
    assert(comm_channel_map.find(tracker.comm) != comm_channel_map.end());
    int comm_sizes[2]={0,0}; int comm_strides[2]={0,0};
    for (auto i=0; i<comm_channel_map[tracker.comm]->id.size(); i++){
      comm_sizes[i]=comm_channel_map[tracker.comm]->id[i].first;
      comm_strides[i]=comm_channel_map[tracker.comm]->id[i].second;
    }
    comm_kernel_key key(rank,-1,tracker.tag,comm_sizes,comm_strides,tracker.nbytes,tracker.partner1);
    if (comm_kernel_info.find(key) == comm_kernel_info.end()){
      comm_kernel_info[key] = std::make_pair(1,comm_time);
    } else{
      comm_kernel_info[key].first++;
      comm_kernel_info[key].second += comm_time;
    }
  }
*/

  std::pair<float,float> cost = cost_model == 0 ?
    tracker.cost_func_bsp(tracker.nbytes, tracker.comm_size)
    : tracker.cost_func_alphabeta(tracker.nbytes, tracker.comm_size);
  if (!is_first_request && cost_model==0) cost.first = 0;

  // Update measurements that define the critical path for each metric.
  cp_costs[0] += cost.second;
  cp_costs[1] += cost.first;
  cp_costs[num_cp_measures-2] += comm_time;
  cp_costs[num_cp_measures-1] += comp_time+comm_time;
  if (path_decomposition == 1){
    for (size_t i=0; i<path_count; i++){
      cp_costs[cp_costs_size-1-i] += comp_time;
    }
  }

  vol_costs[0] += cost.second;
  vol_costs[1] += cost.first;
  vol_costs[num_vol_measures-2] += comm_time;
  vol_costs[num_vol_measures-1] += comp_time+comm_time;
  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
  for (int i=1; i<=5; i++){
    vol_costs[num_vol_measures-i] = std::min(vol_costs[num_vol_measures-i],cp_costs[num_cp_measures-i]);
  }

  // Decompose measurements along multiple paths by MPI routine.
  // Accumuate MPI routine-local measurements. The "my_..." members will never modify the accumulations, while the "cp_..." will first accumulate before path propagation.
  *tracker.my_comm_time += comm_time;
  *tracker.my_msg_count += cost.first;
  *tracker.my_wrd_count += cost.second;
  for (size_t i=0; i<path_count; i++){
    *(tracker.cp_comm_time+i) += comm_time;
    *(tracker.cp_msg_count+i) += cost.first;
    *(tracker.cp_wrd_count+i) += cost.second;
  }

  // Decompose measurements along multiple paths by symbol
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    for (auto i=0; i<path_count; i++){
      for (auto j=0; j<path_measure_index.size(); j++){
        if (path_measure_index[j] == 0){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += cost.second;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += cost.second;
        } else if (path_measure_index[j] == 1){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += cost.first;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += cost.first;
        } else if (path_measure_index[j] == path_measure_select.size()-2){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += comm_time;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += comm_time;
        } else if (path_measure_index[j] == path_measure_select.size()-1){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += (comp_time+comm_time);
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += (comp_time+comm_time);
        }
      }
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[0] += cost.second;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[1] += cost.first;
    symbol_timers[symbol_stack.top()].pp_excl_measure[0] += cost.second;
    symbol_timers[symbol_stack.top()].pp_excl_measure[1] += cost.first;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_decomp_pp_measures-2] += comm_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_decomp_pp_measures-1] += (comp_time+comm_time);
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_decomp_pp_measures-2] += comm_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_decomp_pp_measures-1] += (comp_time+comm_time);
  }

  nonblocking_internal_info.erase(*request);
}

int path::complete_comm(double curtime, MPI_Request* request, MPI_Status* status, int is_test, int* flag){
  auto comp_time = curtime - computation_timer;
  int ret = MPI_SUCCESS;
  auto info_it = nonblocking_internal_info.find(*request);
  assert(info_it != nonblocking_internal_info.end());
  MPI_Request save_request = info_it->first;
  volatile auto last_start_time = MPI_Wtime();
  if (is_test == 0) ret = PMPI_Wait(request, status);
  else{
    ret = PMPI_Test(request,flag,status);
    if (*flag == 0){
      auto save_comm_time = MPI_Wtime() - last_start_time;
      complete_comm(save_comm_time+comp_time);
      computation_timer = MPI_Wtime();
      if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
        symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
      return ret;
    }
  }
  auto save_comm_time = MPI_Wtime() - last_start_time;
  int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
  if (rank != info_it->second.partner){
    // If receiver or collective, complete the barrier and the path data propagation
    if (!info_it->second.is_sender || info_it->second.partner==-1){
      assert(info_it->second.path_data != nullptr);
      MPI_Request req_array[] = {info_it->second.prop_req};
      PMPI_Wait(&req_array[0], MPI_STATUS_IGNORE);
      if (path_decomposition <= 1) update_cp_decomp1(info_it->second.path_data,&cp_costs[0],cp_costs_size);
      else if (path_decomposition == 2) update_cp_decomp2(info_it->second.path_data,&cp_costs[0],cp_costs_size);
      if (info_it->second.path_data != nullptr) free(info_it->second.path_data);
      if (info_it->second.partner == MPI_ANY_SOURCE) { info_it->second.track->partner1 = status->MPI_SOURCE; }
    }
  }
  complete_comm(*info_it->second.track, &save_request, comp_time, save_comm_time);
  computation_timer = MPI_Wtime();
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
  return ret;
}

int path::complete_comm(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status, int is_test, int* flag){

  auto comp_time = curtime - computation_timer;
  int ret = MPI_SUCCESS;
  // We must save the requests before the completition of a request by the MPI implementation
  //   because its tag is set to MPI_REQUEST_NULL and lost forever
  std::vector<MPI_Request> pt(count); for (int i=0;i<count;i++){pt[i]=(array_of_requests)[i];}
  volatile auto last_start_time = MPI_Wtime();
  if (is_test == 0) ret = PMPI_Waitany(count,array_of_requests,indx,status);
  else{
    ret = PMPI_Testany(count,array_of_requests,indx,flag,status);
    if (*flag == 0){
      auto save_comm_time = MPI_Wtime() - last_start_time;
      complete_comm(save_comm_time+comp_time);
      computation_timer = MPI_Wtime();
      if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
        symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
      return ret;
    }
  }
  auto waitany_comm_time = MPI_Wtime() - last_start_time;
  MPI_Request request = pt[*indx];
  auto info_it = nonblocking_internal_info.find(request);
  assert(info_it != nonblocking_internal_info.end());
  if (propagate_p2p){
    int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
    if (rank != info_it->second.partner){
      // If receiver, complete the barrier and the path data propagation
      if (!info_it->second.is_sender || info_it->second.partner==-1){
        assert(info_it->second.path_data != nullptr);
        MPI_Request req_array[] = {info_it->second.prop_req};
        //PMPI_Waitall(2, &req_array[0], MPI_STATUSES_IGNORE);
        PMPI_Wait(&req_array[0], MPI_STATUS_IGNORE);
        if (path_decomposition <= 1) update_cp_decomp1(info_it->second.path_data,&cp_costs[0],cp_costs_size);
        else if (path_decomposition == 2) update_cp_decomp2(info_it->second.path_data,&cp_costs[0],cp_costs_size);
        if (info_it->second.path_data != nullptr) free(info_it->second.path_data);
        if (info_it->second.partner == MPI_ANY_SOURCE) { info_it->second.track->partner1 = status->MPI_SOURCE; }
      }
    }
  }
  complete_comm(*info_it->second.track, &request, comp_time, waitany_comm_time);
  computation_timer = MPI_Wtime();
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
  return ret;
}

int path::complete_comm(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[],
                        MPI_Status array_of_statuses[], int is_test){

  auto waitsome_comp_time = curtime - computation_timer;
  int ret = MPI_SUCCESS;
  is_first_request=true;
  // We must save the requests before the completition of a request by the MPI implementation
  //   because its tag is set to MPI_REQUEST_NULL and lost forever
  std::vector<MPI_Request> pt(incount); for (int i=0;i<incount;i++){pt[i]=(array_of_requests)[i];}
  volatile auto last_start_time = MPI_Wtime();
  if (is_test == 0) ret = PMPI_Waitsome(incount,array_of_requests,outcount,array_of_indices,array_of_statuses);
  else{
    ret = PMPI_Testsome(incount,array_of_requests,outcount,array_of_indices,array_of_statuses);
    if (*outcount == 0){
      auto save_comm_time = MPI_Wtime() - last_start_time;
      complete_comm(save_comm_time+waitsome_comp_time);
      computation_timer = MPI_Wtime();
      if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
        symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
      return ret;
    }
  }
  auto waitsome_comm_time = MPI_Wtime() - last_start_time;
  for (int i=0; i<*outcount; i++){
    MPI_Request request = pt[(array_of_indices)[i]];
    auto info_it = nonblocking_internal_info.find(request);
    assert(info_it != nonblocking_internal_info.end());
    if (propagate_p2p){
      int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
      if (rank != info_it->second.partner){
        // If receiver, complete the barrier and the path data propagation
        if (!info_it->second.is_sender || info_it->second.partner==-1){
          assert(info_it->second.path_data != nullptr);
          MPI_Request req_array[] = {info_it->second.prop_req};
          //PMPI_Waitall(2, &req_array[0], MPI_STATUSES_IGNORE);
          PMPI_Wait(&req_array[0], MPI_STATUS_IGNORE);
          if (path_decomposition <= 1) update_cp_decomp1(info_it->second.path_data,&cp_costs[0],cp_costs_size);
          else if (path_decomposition == 2) update_cp_decomp2(info_it->second.path_data,&cp_costs[0],cp_costs_size);
          if (info_it->second.path_data != nullptr) free(info_it->second.path_data);
          if (info_it->second.partner == MPI_ANY_SOURCE) { info_it->second.track->partner1 = (array_of_statuses)[i].MPI_SOURCE; }
        }
      }
    }
    complete_comm(*info_it->second.track, &request, waitsome_comp_time, waitsome_comm_time);
    waitsome_comp_time=0; waitsome_comm_time=0; is_first_request=false;
  }
  computation_timer = MPI_Wtime();
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
  return ret;
}

int path::complete_comm(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[], int is_test, int* flag){
  auto waitall_comp_time = curtime - computation_timer;
  int ret = MPI_SUCCESS;
  is_first_request=true;
  // We must save the requests before the completition of a request by the MPI implementation
  //   because its tag is set to MPI_REQUEST_NULL and lost forever
  std::vector<MPI_Request> pt(count); for (int i=0;i<count;i++){pt[i]=(array_of_requests)[i];}
  volatile auto last_start_time = MPI_Wtime();
  if (is_test == 0) ret = PMPI_Waitall(count,array_of_requests,array_of_statuses);
  else{
    ret = PMPI_Testall(count,array_of_requests,flag,array_of_statuses);
    if (*flag == 0){
      auto save_comm_time = MPI_Wtime() - last_start_time;
      complete_comm(save_comm_time+waitall_comp_time);
      computation_timer = MPI_Wtime();
      if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
        symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
      return ret;
    }
  }
  auto waitall_comm_time = MPI_Wtime() - last_start_time;
  for (int i=0; i<count; i++){
    MPI_Request request = pt[i];
    auto info_it = nonblocking_internal_info.find(request);
    assert(info_it != nonblocking_internal_info.end());
    if (propagate_p2p){
      int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
      if (rank != info_it->second.partner){
        // If receiver, complete the barrier and the path data propagation
        if (!info_it->second.is_sender || info_it->second.partner==-1){
          assert(info_it->second.path_data != nullptr);
          MPI_Request req_array[] = {info_it->second.prop_req};
          //PMPI_Waitall(2, &req_array[0], MPI_STATUSES_IGNORE);
          PMPI_Wait(&req_array[0], MPI_STATUS_IGNORE);
          if (path_decomposition <= 1) update_cp_decomp1(info_it->second.path_data,&cp_costs[0],cp_costs_size);
          else if (path_decomposition == 2) update_cp_decomp2(info_it->second.path_data,&cp_costs[0],cp_costs_size);
          if (info_it->second.path_data != nullptr) free(info_it->second.path_data);
          if (info_it->second.partner == MPI_ANY_SOURCE) { info_it->second.track->partner1 = (array_of_statuses)[i].MPI_SOURCE; }
        }
      }
    }
    complete_comm(*info_it->second.track, &request, waitall_comp_time, waitall_comm_time);
    // Although we have to exchange the path data for each request, we do not want to float-count the computation time nor the communicaion time
    waitall_comp_time=0; waitall_comm_time=0; is_first_request=false;
  }
  computation_timer = MPI_Wtime();
  if (path_decomposition == 2 && path_count>0 && symbol_stack.size()>0){
    symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
  return ret;
}

void path::propagate(blocking& tracker){
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  if ((rank == tracker.partner1) && (rank == tracker.partner2)) { return; } 
  // Exchange the tracked routine critical path data
  if (tracker.partner1 == -1){
    if (propagate_collective==0) return;
    MPI_Op op = MPI_MAX;
    if (path_decomposition <= 1 && path_count>0) MPI_Op_create((MPI_User_function*) propagate_cp_decomp1_op,0,&op);
    else if (path_decomposition == 2) MPI_Op_create((MPI_User_function*) propagate_cp_decomp2_op,0,&op);
    PMPI_Allreduce(MPI_IN_PLACE, &cp_costs[0], cp_costs.size(), MPI_FLOAT, op, tracker.comm);
    //PMPI_Allreduce(MPI_IN_PLACE, &cp_costs[0], 5+51*5, MPI_FLOAT, MPI_MAX, tracker.comm);
    if (path_count>0 || path_decomposition >= 2) MPI_Op_free(&op);
/*
      size_t active_size = tracker.comm_size;
      size_t active_rank = rank;
      size_t active_mult = 1;
      while (active_size>1){
	if (active_rank % 2 == 1){
          MPI_Buffer_attach(&eager_pad[0],eager_pad.size());
	  int partner = (active_rank-1)*active_mult;
          PMPI_Bsend(&cp_costs[0], cp_costs.size(), MPI_FLOAT, partner, internal_tag2, tracker.comm);
          void* temp_buf; int temp_size;
          MPI_Buffer_detach(&temp_buf,&temp_size);
	  break;
	}
	else if ((active_rank % 2 == 0) && (active_rank < (active_size-1))){
	  int partner = (active_rank+1)*active_mult;
          PMPI_Recv(&cp_costs_foreign[0], cp_costs.size(), MPI_FLOAT, partner, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
          if (path_decomposition <= 1) update_cp_decomp1(&cp_costs_foreign[0],&cp_costs[0],cp_costs_size);
          else if (path_decomposition == 2) update_cp_decomp2(&cp_costs_foreign[0],&cp_costs[0],cp_costs_size);
	}
	active_size = active_size/2 + active_size%2;
	active_rank /= 2;
	active_mult *= 2;
      }
    PMPI_Barrier(tracker.comm);// necessary to track communication time accuracy (else some processes will leave early and wait).
*/
  } else{
    if (propagate_p2p==0) return;
    MPI_Buffer_attach(&eager_pad[0],eager_pad.size());
    if (tracker.is_sender){
      PMPI_Bsend(&cp_costs[0], cp_costs.size(), MPI_FLOAT, tracker.partner1, internal_tag2, tracker.comm);
    } else{
      PMPI_Recv(&cp_costs_foreign[0], cp_costs.size(), MPI_FLOAT, tracker.partner1, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
      if (path_decomposition <= 1) update_cp_decomp1(&cp_costs_foreign[0],&cp_costs[0],cp_costs_size);
      else if (path_decomposition == 2) update_cp_decomp2(&cp_costs_foreign[0],&cp_costs[0],cp_costs_size);
    }
    if (tracker.partner2 != tracker.partner1){
      PMPI_Recv(&cp_costs_foreign[0], cp_costs.size(), MPI_FLOAT, tracker.partner2, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
      if (path_decomposition <= 1) update_cp_decomp1(&cp_costs_foreign[0],&cp_costs[0],cp_costs_size);
      else if (path_decomposition == 2) update_cp_decomp2(&cp_costs_foreign[0],&cp_costs[0],cp_costs_size);
    }
  }
}

void path::propagate(nonblocking& tracker, float*& path_data, MPI_Request* prop_req){
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  if (rank == tracker.partner1) { return; } 
  // Exchange the tracked routine critical path data
  if (tracker.partner1 == -1){
    MPI_Op op;
    if (path_decomposition <= 1) MPI_Op_create((MPI_User_function*) propagate_cp_decomp1_op,0,&op);
    else if (path_decomposition == 2) MPI_Op_create((MPI_User_function*) propagate_cp_decomp2_op,0,&op);
    path_data = (float*)malloc(cp_costs.size()*sizeof(float));
    std::memcpy(path_data, &cp_costs[0], cp_costs.size()*sizeof(float));
    PMPI_Iallreduce(MPI_IN_PLACE,path_data,cp_costs.size(),MPI_FLOAT,op,tracker.comm,prop_req);
    //MPI_Op_free(&op);
  }
  else{
    if (tracker.is_sender){
      PMPI_Bsend(&cp_costs[0], cp_costs.size(), MPI_FLOAT, tracker.partner1, internal_tag2, tracker.comm);
    }
    else{
      path_data = (float*)malloc(cp_costs.size()*sizeof(float));
      PMPI_Irecv(path_data, cp_costs.size(), MPI_FLOAT, tracker.partner1, internal_tag2, tracker.comm, prop_req);
    }
  }
}

}
}
}
