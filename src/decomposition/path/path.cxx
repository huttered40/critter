#include "path.h"
#include "../container/symbol_tracker.h"
#include "../../discretization/util/util.h"
#include "../util/util.h"
#include "../../util/util.h"

namespace critter{
namespace internal{
namespace decomposition{

void path::exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm){
  // Save and accumulate the computation time between last communication routine as both execution-time and computation time
  //   into both the execution-time critical path data structures and the per-process data structures.
  float save_comp_time = MPI_Wtime() - computation_timer;
  critical_path_costs[num_critical_path_measures-1] += save_comp_time;	// update critical path execution time
  critical_path_costs[num_critical_path_measures-3] += save_comp_time;	// update critical path computation time
  volume_costs[num_volume_measures-1]        += save_comp_time;		// update local execution time
  volume_costs[num_volume_measures-3]        += save_comp_time;		// update local computation time
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-comm_path_select_size-1-i] += save_comp_time; }// update each metric's critical path's computation time

  generate_aggregate_channels(oldcomm,newcomm);
  PMPI_Barrier(oldcomm);
  computation_timer = MPI_Wtime();
}

bool path::initiate_comp(size_t id, volatile float curtime, float flop_count, int param1, int param2, int param3, int param4, int param5){
  // accumulate computation time
  float save_comp_time = curtime - computation_timer;
  critical_path_costs[num_critical_path_measures-3] += save_comp_time;	// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += save_comp_time;	// update critical path runtime
  volume_costs[num_volume_measures-3]        += save_comp_time;		// update local computation time
  volume_costs[num_volume_measures-1]        += save_comp_time;		// update local runtime
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-comm_path_select_size-1-i] += save_comp_time; }// update each metric's critical path's computation time
  if (symbol_path_select_size>0 && symbol_stack.size()>0){
    // Get the current symbol's execution-time since last communication routine or its inception.
    // Accumulate as both execution-time and computation time into both the execution-time critical path data structures and the per-process data structures.
    auto last_symbol_time = curtime - symbol_timers[symbol_stack.top()].start_timer.top();
    for (auto i=0; i<symbol_path_select_size; i++){
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-1] += last_symbol_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-3] += last_symbol_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-1] += last_symbol_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-3] += last_symbol_time;
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-3] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-3] += last_symbol_time;
  }

  // start compunication timer for compunication routine
  comp_start_time = MPI_Wtime();
  return true;
}

void path::complete_comp(float errtime, size_t id, float flop_count, int param1, int param2, int param3, int param4, int param5){
  volatile float comp_time = MPI_Wtime() - comp_start_time - errtime;	// complete computation time

  // Save kernel information
  if (autotuning_debug == 1){
    comp_kernel_key key(-1,id,flop_count,param1,param2,param3,param4,param5);// '-1' argument is arbitrary, does not influence overloaded operators
    if (comp_kernel_info.find(key) == comp_kernel_info.end()){
      comp_kernel_info[key] = std::make_pair(1,comp_time);
    } else{
      comp_kernel_info[key].first++;
      comp_kernel_info[key].second += comp_time;
    }
  }

  // Decompose measurements along multiple paths by symbol
  if (symbol_path_select_size>0 && symbol_stack.size()>0){
    for (auto i=0; i<symbol_path_select_size; i++){
      // update all computation-related measures for the top symbol in stack
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-7] += flop_count;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-3] += comp_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-2] += comp_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-1] += comp_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-7] += flop_count;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-3] += comp_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-2] += comp_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-1] += comp_time;
    }
    // update all computation-related measures for the top symbol in stack
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-7] += flop_count;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-3] += comp_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-2] += comp_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += comp_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-7] += flop_count;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-3] += comp_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-2] += comp_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += comp_time;
  }
  critical_path_costs[num_critical_path_measures-6] += flop_count;//'6', not '7' because 'critical_path_costs' does not store idle time
  critical_path_costs[num_critical_path_measures-3] += comp_time;
  critical_path_costs[num_critical_path_measures-2] += comp_time;
  critical_path_costs[num_critical_path_measures-1] += comp_time;
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-1-i] += comp_time; }// update each metric's critical path's computation kernel time
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-comm_path_select_size-1-i] += comp_time; }// update each metric's critical path's computation time
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-3*comm_path_select_size-1-i] += flop_count; }// update each metric's critical path's computation time

  volume_costs[num_volume_measures-7] += flop_count;
  volume_costs[num_volume_measures-3] += comp_time;
  volume_costs[num_volume_measures-2] += comp_time;
  volume_costs[num_volume_measures-1] += comp_time;

  computation_timer = MPI_Wtime();
  if (symbol_path_select_size>0 && symbol_stack.size()>0){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
}

static void update_critical_path(float* in, float* inout, size_t len){
  assert(len == critical_path_costs_size);	// this assert prevents user from obtaining wrong output if MPI implementation cuts up the message.
  if (comm_path_select_size > 0){
    size_t breakdown_idx=0;
    for (int i=0; i<num_critical_path_measures; i++){
      if (comm_path_select[i]=='1'){ decisions[breakdown_idx++] = inout[i] > in[i]; }
    }
    for (int i=0; i<num_critical_path_measures; i++){
      inout[i] = std::max(inout[i],in[i]);
    }
    for (int i=num_critical_path_measures; i<critical_path_costs_size; i++){
      int idx = (i-num_critical_path_measures)%comm_path_select_size;
      inout[i] = (decisions[idx] ? inout[i] : in[i]);
    }
  } else{
    for (int i=0; i<num_critical_path_measures; i++){
      inout[i] = std::max(inout[i],in[i]);
    }
  }
}

static void propagate_critical_path_op(float* in, float* inout, int* len, MPI_Datatype* dtype){
  update_critical_path(in,inout,static_cast<size_t>(*len));
}

static void complete_timers(float* remote_path_data, size_t msg_id){
  // Note that only receives will enter this function anyways, so no need to branch
  int* envelope_int[2] = { internal_timer_prop_int[2*msg_id], internal_timer_prop_int[2*msg_id+1] };
  float* envelope_float[2] = { remote_path_data, internal_timer_prop_float[msg_id] };
  char* envelope_char = internal_timer_prop_char[msg_id];
  for (auto k=0; k<symbol_path_select_size; k++){
    // Up until this very point, we had no idea whether we, or our partner rank, determined the path for a specific metric.
    if (envelope_float[0][symbol_path_select_index[k]] > critical_path_costs[symbol_path_select_index[k]]){
      int ftimer_size = *envelope_int[0];
      int symbol_offset = 0;
      for (int i=0; i<ftimer_size; i++){
        auto reconstructed_symbol = std::string(envelope_char+symbol_offset,envelope_char+symbol_offset+envelope_int[1][i]);
        if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
          symbol_timers[reconstructed_symbol] = symbol_tracker(reconstructed_symbol);
          symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
        }
        std::memcpy(symbol_timers[reconstructed_symbol].cp_numcalls[k],
                    &envelope_float[1][(i*symbol_path_select_size+k)*(cp_symbol_class_count*num_per_process_measures+1)],
                    sizeof(float)*(cp_symbol_class_count*num_per_process_measures+1));
        symbol_timers[reconstructed_symbol].has_been_processed = true;
        symbol_offset += envelope_int[1][i];
      }
      // Now cycle through and find the symbols that were not processed and set their accumulated measures to 0
      for (auto& it : symbol_timers){
        if (it.second.has_been_processed){ it.second.has_been_processed = false; }
        else{
          it.second.cp_numcalls[k][0] = 0;
          for (int j=0; j<num_per_process_measures; j++){
            it.second.cp_incl_measure[k][j] = 0;
            it.second.cp_excl_measure[k][j] = 0;
          }
        }
      }
    }
  }
}

static void complete_path_update(){
  PMPI_Waitall(internal_comm_prop_req.size(), &internal_comm_prop_req[0], MPI_STATUSES_IGNORE);
  if (symbol_path_select_size>0) { PMPI_Waitall(internal_timer_prop_req.size(), &internal_timer_prop_req[0], MPI_STATUSES_IGNORE); }
  size_t msg_id=0;
  for (auto& it : internal_comm_prop){
    if (!it.second){
      if (symbol_path_select_size>0) complete_timers(it.first,msg_id++);
      update_critical_path(it.first,&critical_path_costs[0],critical_path_costs_size);
    }
    free(it.first);
  }
  internal_comm_prop.clear(); internal_comm_prop_req.clear();
  for (auto& it : internal_timer_prop_int){ free(it); }
  for (auto& it : internal_timer_prop_float){ free(it); }
  for (auto& it : internal_timer_prop_float_int){ free(it); }
  for (auto& it : internal_timer_prop_char){ free(it); }
  internal_timer_prop_int.clear(); internal_timer_prop_float.clear(); internal_timer_prop_float_int.clear();
  internal_timer_prop_char.clear(); internal_timer_prop_req.clear();
}


bool path::initiate_comm(blocking& tracker, volatile float curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                            bool is_sender, int partner1, int partner2){
  // Save and accumulate the computation time between last communication routine as both execution-time and computation time
  //   into both the execution-time critical path data structures and the per-process data structures.
  tracker.comp_time = curtime - computation_timer;
  // Note that use of MPI_ANY_SOURCE may in fact be ok, but it has not been tested.
  assert(partner1 != MPI_ANY_SOURCE); if ((tracker.tag == 13) || (tracker.tag == 14)){ assert(partner2 != MPI_ANY_SOURCE); }
  int rank; MPI_Comm_rank(comm, &rank);
  MPI_Buffer_attach(&eager_pad[0],eager_pad.size());

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
  tracker.synch_time = 0.;// might get updated below
  tracker.barrier_time=0.;// might get updated below
  // Below: tricky situation in which the receiver assumes eager protocol due to message size, but sender is synchronous.
  if (tracker.tag == 15) assert(tracker.nbytes > eager_limit);
  bool eager = ((tracker.partner1 != -1) && (tracker.tag!=13) && (tracker.tag!=14) && (tracker.tag!=15) && (tracker.nbytes <= eager_limit));

  // Use a barrier or synchronous (rendezvous protocol) send/recv to track idle time (i.e. one process will be the latest to arrive at this segment of code, thus all other processes directly wait for it)
  // This avoids corruption of communication time when processes are still waiting for the initial synchronization to proceed with the communication.
  // Note that unlike the execution-time critical path, critical paths defined by other metrics besides execution-time can incur idle time.

  // Note that the only reason we separate out the blocking p2p idle time communication is due to the possibility that the other side being nonblocking.
  // Nonblocking sends can not utilize sendrecv because they may be embedded within a Waitall that does not guarantee ordering.
  // To account for this for blocking sends, we use the synchronous send to avoid automatic eager-protocol (which would have allowed early exit from p2p barrier).
  // The only reason we do not always issue a sendrecv with partners 1&2 is because of the possibility in which the other side issued a nonblocking communication request and our handling for idle time there
  //   is for a nonblocking sender to issue a send and a nonblocking receiver to issue a recv, rather than both. This also is more aligned with a future goal of utilizing eager sends, which might be constitute a different mechanism.
  // Note that we favor {Issend,Irecv} rather than {Ssend,Recv,Sendrecv} only because it simplifies logic when handling multiple possible two-sided p2p communication patterns.
  //   A user Sendrecv cannot be handled with separate Send+recv because a Sendrecv is implemented via nonblocking p2p in all MPI implementations as it is a construct used in part to prevent deadlock.

  volatile float init_time = MPI_Wtime();
  if (partner1 == -1){
    PMPI_Barrier(tracker.comm);
    tracker.barrier_time = MPI_Wtime() - init_time;
  }
  else {
    MPI_Request barrier_reqs[3]; int barrier_count=0;
    char sbuf='H'; char rbuf='H';
    if ((is_sender) && (rank != partner1)){
      if (eager) { PMPI_Bsend(&sbuf, 1, MPI_CHAR, partner1, internal_tag3, tracker.comm); }
      else       { PMPI_Issend(&sbuf, 1, MPI_CHAR, partner1, internal_tag3, tracker.comm, &barrier_reqs[barrier_count]); barrier_count++; }
    }
    if ((!is_sender) && (rank != partner1)){
      PMPI_Irecv(&rbuf, 1, MPI_CHAR, partner1, internal_tag3, tracker.comm, &barrier_reqs[barrier_count]); barrier_count++;
    }
    if ((partner2 != -1) && (rank != partner2)){
      PMPI_Irecv(&rbuf, 1, MPI_CHAR, partner2, internal_tag3, tracker.comm, &barrier_reqs[barrier_count]); barrier_count++;
    }
    PMPI_Waitall(barrier_count,&barrier_reqs[0],MPI_STATUSES_IGNORE);
    if (barrier_count>0) tracker.barrier_time = MPI_Wtime() - init_time;
  }
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-1-i-2*comm_path_select_size] += tracker.barrier_time; }

  critical_path_costs[num_critical_path_measures-3] += tracker.comp_time;	// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += tracker.comp_time;	// update critical path runtime
  volume_costs[num_volume_measures-3]        += tracker.comp_time;		// update local computation time
  volume_costs[num_volume_measures-1]        += tracker.comp_time;		// update local runtime
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-1-i-comm_path_select_size] += tracker.comp_time; }// update each metric's critical path's computation time
  if (symbol_path_select_size>0 && symbol_stack.size()>0){
    // Get the current symbol's execution-time since last communication routine or its inception.
    // Accumulate as both execution-time and computation time into both the execution-time critical path data structures and the per-process data structures.
    auto last_symbol_time = curtime - symbol_timers[symbol_stack.top()].start_timer.top();
    for (auto i=0; i<symbol_path_select_size; i++){
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-1] += last_symbol_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-3] += last_symbol_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-1] += last_symbol_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-3] += last_symbol_time;
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-3] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-3] += last_symbol_time;
  }

  if (track_synchronization && tracker.partner1==-1){
    // Use the user communication routine to measre synchronization time.
    // Note the following consequences of using a tiny 1-byte message (note that 0-byte is trivially handled by most MPI implementations) on measuring synchronization time:
    // 	1) The collective communication algorithm is likely different for small messages than large messages.
    // 	2) The eager sending protocol will be utilized, which would incur a potentially significant difference in synchronization time than if rendezvous protocol was invoked.

    // Special arrays for use in the collective -v routines as well as Reduce_scatter.
    std::vector<int> counts(np,1); std::vector<int> disp(np,0);
    if (tracker.tag>=9 && tracker.tag<=12) for (int i=1; i<np; i++) disp[i]=disp[i-1]+1;
    // start synchronization timer for communication routine
    tracker.start_time = MPI_Wtime();
    switch (tracker.tag){
      case 0:
        PMPI_Barrier(tracker.comm);
        break;
      case 1:
        PMPI_Bcast(&synch_pad_send[0], 1, MPI_CHAR, 0, tracker.comm);// arbitrary root 0
        break;
      case 2:
        PMPI_Reduce(&synch_pad_send[0], &synch_pad_recv[0], 1, MPI_CHAR, MPI_MAX, 0, tracker.comm);// arbitrary root 0
        break;
      case 3:
        PMPI_Allreduce(MPI_IN_PLACE, &synch_pad_send[0], 1, MPI_CHAR, MPI_MAX, tracker.comm);
        break;
      case 4:
        PMPI_Gather(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, 0, tracker.comm);// arbitrary root 0
        break;
      case 5:
        PMPI_Allgather(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, tracker.comm);
        break;
      case 6:
        PMPI_Scatter(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, 0, tracker.comm);// arbitrary root 0
        break;
      case 7:
        PMPI_Reduce_scatter(&synch_pad_send[0], &synch_pad_recv[0], &counts[0], MPI_CHAR, MPI_MAX, tracker.comm);
        break;
      case 8:
        PMPI_Alltoall(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, tracker.comm);
        break;
      case 9:
        PMPI_Gatherv(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], &counts[0], &disp[0], MPI_CHAR, 0, tracker.comm);// arbitrary root 0
        break;
      case 10:
        PMPI_Allgatherv(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], &counts[0], &disp[0], MPI_CHAR, tracker.comm);
        break;
      case 11:
        PMPI_Scatterv(&synch_pad_send[0], &counts[0], &disp[0], MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, 0, tracker.comm);// arbitrary root 0
        break;
      case 12:
        PMPI_Alltoallv(&synch_pad_send[0], &counts[0], &disp[0], MPI_CHAR, &synch_pad_recv[0], &counts[0], &disp[0], MPI_CHAR, tracker.comm);
        break;
    }
    tracker.synch_time = MPI_Wtime()-tracker.start_time;
  }

  void* temp_buf; int temp_size;
  MPI_Buffer_detach(&temp_buf,&temp_size);
  // start communication timer for communication routine
  tracker.start_time = MPI_Wtime();
  return true;
}

void path::complete_comm(blocking& tracker, int recv_source){
  volatile float comm_time = MPI_Wtime() - tracker.start_time;	// complete communication time
  int rank; MPI_Comm_rank(tracker.comm, &rank);
  MPI_Buffer_attach(&eager_pad[0],eager_pad.size());

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

  std::pair<float,float> cost_bsp = tracker.cost_func_bsp(tracker.nbytes, tracker.comm_size);
  std::pair<float,float> cost_alphabeta = tracker.cost_func_alphabeta(tracker.nbytes, tracker.comm_size);
  std::vector<std::pair<float,float>> costs = {cost_bsp,cost_alphabeta};

  // We handle wildcard sources (for MPI_Recv variants) only after the user communication.
  if (recv_source != -1){
    if ((tracker.tag == 13) || (tracker.tag == 14)){ tracker.partner2=recv_source; }
    else{ assert(tracker.tag==17); tracker.partner1=recv_source; }
  }

  // Decompose measurements along multiple paths by MPI routine.
  // Accumuate MPI routine-local measurements. The "my_..." members will never modify the accumulations, while the "critical_path_..." will first accumulate before path propagation.
  *tracker.my_synch_time   += tracker.synch_time;
  *tracker.my_comm_time    += comm_time;
  int save=0;
  for (int j=0; j<cost_models.size(); j++){
    if (cost_models[j]=='1'){
      *(tracker.my_msg_count+save) += costs[j].first;
      *(tracker.my_wrd_count+save) += costs[j].second;
      save++;
    }
  }
  for (size_t i=0; i<comm_path_select_size; i++){
    *(tracker.critical_path_synch_time+i)   += tracker.synch_time;
    *(tracker.critical_path_comm_time+i)    += comm_time;
  }
  save=0;
  for (int j=0; j<cost_models.size(); j++){
    for (size_t i=0; i<comm_path_select_size; i++){
      if (cost_models[j]=='1'){
        *(tracker.critical_path_msg_count+save*comm_path_select_size+i) += costs[j].first;
        *(tracker.critical_path_wrd_count+save*comm_path_select_size+i) += costs[j].second;
      }
    }
    save++;
  }

  // Decompose measurements along multiple paths by symbol
  if (symbol_path_select_size>0 && symbol_stack.size()>0){
    for (auto i=0; i<symbol_path_select_size; i++){
      // update all communication-related measures for the top symbol in stack
      size_t save=0;
      for (int j=0; j<cost_models.size(); j++){
        if (cost_models[j]=='1'){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][save] += costs[j].second;
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][cost_models.size()+save] += costs[j].first;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][save] += costs[j].second;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][cost_models.size()+save] += costs[j].first;
        }
        save++;
      }
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-6] += tracker.barrier_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-5] += comm_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-4] += tracker.synch_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-1] += comm_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-6] += tracker.barrier_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-5] += comm_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-4] += tracker.synch_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-1] += comm_time;
    }
    // update all communication-related measures for the top symbol in stack
    size_t save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]=='1'){
        symbol_timers[symbol_stack.top()].pp_exclusive_measure[save] += costs[j].second;
        symbol_timers[symbol_stack.top()].pp_exclusive_measure[cost_models.size()+save] += costs[j].first;
        symbol_timers[symbol_stack.top()].pp_excl_measure[save] += costs[j].second;
        symbol_timers[symbol_stack.top()].pp_excl_measure[cost_models.size()+save] += costs[j].first;
      }
      save++;
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-6] += tracker.barrier_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-5] += comm_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-4] += tracker.synch_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += (comm_time+tracker.barrier_time);
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-6] += tracker.barrier_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-5] += comm_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-4] += tracker.synch_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += (comm_time+tracker.barrier_time);
  }

  // Update measurements that define the critical path for each metric.
  save=0;
  for (int j=0; j<cost_models.size(); j++){
    if (cost_models[j]=='1'){
      critical_path_costs[save]                 += costs[j].second;		// update critical path estimated communication cost
      critical_path_costs[cost_model_size+save] += costs[j].first;		// update critical path estimated synchronization cost
      volume_costs[save]                        += costs[j].second;		// update local estimated communication cost
      volume_costs[cost_model_size+save]        += costs[j].first;		// update local estimated synchronization cost
      save++;
    }
  }
  critical_path_costs[num_critical_path_measures-5] += comm_time;		// update critical path communication time (for what this process has seen thus far)
  critical_path_costs[num_critical_path_measures-4] += tracker.synch_time;	// update critical path synchronization time
  critical_path_costs[num_critical_path_measures-1] += comm_time;		// update critical path runtime

  volume_costs[num_volume_measures-6] += tracker.barrier_time;			// update local barrier/idle time
  volume_costs[num_volume_measures-5] += comm_time;				// update local communication time (not volume until after the completion of the program)
  volume_costs[num_volume_measures-4] += tracker.synch_time;			// update local synchronization time
  volume_costs[num_volume_measures-1] += (tracker.barrier_time+comm_time);	// update local runtime with idle time and comm time

  // Note that this block of code below is left in solely for blocking communication to avoid over-counting the idle time
  //   (which does not get subtracted by the min idle time any one process incurs due to efficiency complications with matching nonblocking+blocking p2p communications).
  //   Its handled correctly for blocking collectives.
  // If per-process execution-time gets larger than execution-time along the execution-time critical path, subtract out the difference from idle time.
  volume_costs[num_volume_measures-6] -= std::max((float)0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
  if (symbol_path_select_size>0 && symbol_stack.size()>0){
    // Special handling of excessively large idle time caused by suspected tool interference
    // Specifically, this interference is caused by not subtracting out the barrier time of the last process to enter the barrier (which ideally is 0).
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] -= std::max((float)0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1]     -= std::max((float)0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-6] -= std::max((float)0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-6]     -= std::max((float)0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
  }

  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
  volume_costs[num_volume_measures-5] = volume_costs[num_volume_measures-5] > critical_path_costs[num_critical_path_measures-5]
                                          ? critical_path_costs[num_critical_path_measures-5] : volume_costs[num_volume_measures-5];
  volume_costs[num_volume_measures-4] = volume_costs[num_volume_measures-4] > critical_path_costs[num_critical_path_measures-4]
                                          ? critical_path_costs[num_critical_path_measures-4] : volume_costs[num_volume_measures-4];
  volume_costs[num_volume_measures-3] = volume_costs[num_volume_measures-3] > critical_path_costs[num_critical_path_measures-3]
                                          ? critical_path_costs[num_critical_path_measures-3] : volume_costs[num_volume_measures-3];
  volume_costs[num_volume_measures-2] = volume_costs[num_volume_measures-2] > critical_path_costs[num_critical_path_measures-2]
                                          ? critical_path_costs[num_critical_path_measures-2] : volume_costs[num_volume_measures-2];
  volume_costs[num_volume_measures-1] = volume_costs[num_volume_measures-1] > critical_path_costs[num_critical_path_measures-1]
                                          ? critical_path_costs[num_critical_path_measures-1] : volume_costs[num_volume_measures-1];

  // Propogate critical paths for all processes in communicator based on what each process has seen up until now (not including this communication)
  propagate(tracker);

  void* temp_buf; int temp_size;
  MPI_Buffer_detach(&temp_buf,&temp_size);

  // Prepare to leave interception and re-enter user code by restarting computation timers.
  bsp_counter++;
  tracker.start_time = MPI_Wtime();
  computation_timer = tracker.start_time;
  if (symbol_path_select_size>0 && symbol_stack.size()>0){ symbol_timers[symbol_stack.top()].start_timer.top() = tracker.start_time; }
}

// Called by both nonblocking p2p and nonblocking collectives
bool path::initiate_comm(nonblocking& tracker, volatile float curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm, bool is_sender, int partner){

  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  tracker.comp_time = curtime - computation_timer;
  // Note that use of MPI_ANY_SOURCE may in fact be ok, but it has not been tested.
  assert(partner != MPI_ANY_SOURCE);

  critical_path_costs[num_critical_path_measures-3] += tracker.comp_time;		// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += tracker.comp_time;		// update critical path runtime
  volume_costs[num_volume_measures-3]        += tracker.comp_time;		// update local computation time
  volume_costs[num_volume_measures-1]        += tracker.comp_time;		// update local runtime
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-1-i-comm_path_select_size] += tracker.comp_time; }
  if (symbol_path_select_size>0 && symbol_stack.size()>0){
    assert(symbol_stack.size()>0);
    assert(symbol_timers[symbol_stack.top()].start_timer.size()>0);
    float save_time = curtime - symbol_timers[symbol_stack.top()].start_timer.top();
    for (auto i=0; i<symbol_path_select_size; i++){
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-1] += save_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-3] += save_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-1] += save_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-3] += save_time;
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += save_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-3] += save_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-1] += save_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-3] += save_time;
  }

  // Note: routine below will be called immediately afterward.
  return true;
}

// Called by both nonblocking p2p and nonblocking collectives
void path::initiate_comm(nonblocking& tracker, volatile float itime, int64_t nelem,
                            MPI_Datatype t, MPI_Comm comm, MPI_Request* request, bool is_sender, int partner){

  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  tracker.comp_time = itime;
  critical_path_costs[num_critical_path_measures-3] += tracker.comp_time;		// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += tracker.comp_time;		// update critical path runtime
  volume_costs[num_volume_measures-3]        += tracker.comp_time;		// update local computation time
  volume_costs[num_volume_measures-1]        += tracker.comp_time;		// update local runtime
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-1-i-comm_path_select_size] += tracker.comp_time; }
  if (symbol_path_select_size>0 && symbol_stack.size()>0){
    assert(symbol_stack.size()>0);
    assert(symbol_timers[symbol_stack.top()].start_timer.size()>0);
    float save_time = itime;
    for (auto i=0; i<symbol_path_select_size; i++){
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-1] += save_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-3] += save_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-1] += save_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-3] += save_time;
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += save_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-3] += save_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-1] += save_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-3] += save_time;
  }

  int el_size,p;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(comm, &p);
  bool eager = tracker.nbytes <= eager_limit;
  nonblocking_info msg_info(is_sender,partner,comm,(float)nbytes,(float)p,&tracker);
  nonblocking_internal_info[*request] = msg_info;

  MPI_Buffer_attach(&eager_pad[0],eager_pad.size());
  // Issue the barrier/max-barrier/synchronization calls, regardless of msg size
  if (partner!=-1){// Branch protects against nonblocking collectives
    int rank; MPI_Comm_rank(comm,&rank); 
    float max_barrier_time = 0;// counter-intuitively, a blocking partner should determine the idle time
    if (is_sender && rank != partner){
      PMPI_Bsend(&barrier_pad_send[0], 1, MPI_CHAR, partner, internal_tag3, comm);
    }
    else if (!is_sender && rank != partner){
      // The motivation behind use of PMPI_Recv here is that there is every chance that the Sender would have sent this already and its cheap to receive it right here and now.
      PMPI_Recv(&barrier_pad_recv[0], 1, MPI_CHAR, partner, internal_tag3, comm, MPI_STATUS_IGNORE);
    }
  }

  tracker.comm = comm;
  tracker.is_sender = is_sender;
  tracker.partner1 = partner;
  tracker.partner2 = -1;
  tracker.comm_size = p;
  tracker.nbytes = nbytes;
  propagate(tracker);

  // Forces buffered messages to send. Ideally we should wait till the next invocation of 'path::initiate(blocking&,...)' to call this,
  //   but to be safe and avoid stalls caused by MPI implementation not sending until this routine is called, we call it here.
  void* temp_buf; int temp_size;
  MPI_Buffer_detach(&temp_buf,&temp_size);
  tracker.start_time = MPI_Wtime();
  computation_timer = tracker.start_time;
  if (symbol_path_select_size>0 && symbol_stack.size()>0){ symbol_timers[symbol_stack.top()].start_timer.top() = tracker.start_time; }
}

void path::complete_comm(nonblocking& tracker, MPI_Request* request, float comp_time, float comm_time){
  auto info_it = nonblocking_internal_info.find(*request);
  assert(info_it != nonblocking_internal_info.end());

  tracker.is_sender = info_it->second.is_sender;
  tracker.comm = info_it->second.comm;
  tracker.partner1 = info_it->second.partner;
  tracker.partner2 = -1;
  tracker.nbytes = info_it->second.nbytes;
  tracker.comm_size = info_it->second.comm_size;
  tracker.synch_time=0;
  bool eager = tracker.nbytes <= eager_limit;

  int rank; MPI_Comm_rank(tracker.comm, &rank);
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

  // Both sender and receiver will now update its critical path with the data from the communication
  std::pair<float,float> cost_bsp  = tracker.cost_func_bsp(tracker.nbytes,tracker.comm_size);
  std::pair<float,float> cost_alphabeta = tracker.cost_func_alphabeta(tracker.nbytes,tracker.comm_size);
  if (is_first_request) cost_bsp.first=1.;	// this is usually zero, but we force it to be 1 in special circumstances (for nonblocking p2p with wait_id one)
  std::vector<std::pair<float,float>> costs = {cost_bsp,cost_alphabeta};

  // Update measurements that define the critical path for each metric.
  int save=0;
  for (int j=0; j<cost_models.size(); j++){
    if (cost_models[j]=='1'){
      critical_path_costs[save]                 += costs[j].second;	// update critical path estimated communication cost
      critical_path_costs[cost_model_size+save] += costs[j].first;	// update critical path estimated synchronization cost
      volume_costs[save]                 += costs[j].second;		// update local estimated communication cost
      volume_costs[cost_model_size+save] += costs[j].first;		// update local estimated synchronization cost
      save++;
    }
  }
  critical_path_costs[num_critical_path_measures-5] += comm_time;			// update critical path communication time (for what this process has seen thus far)
  critical_path_costs[num_critical_path_measures-4] += 0.;				// update critical path synchronization time
  critical_path_costs[num_critical_path_measures-3] += comp_time;			// update critical path runtime
  critical_path_costs[num_critical_path_measures-1] += comp_time+comm_time;		// update critical path runtime
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-1-i-comm_path_select_size] += comp_time; }

  volume_costs[num_volume_measures-5] += comm_time;				// update local communication time (not volume until after the completion of the program)
  volume_costs[num_volume_measures-4] += 0.;					// update local synchronization time
  volume_costs[num_volume_measures-3] += comp_time;				// update local runtime
  volume_costs[num_volume_measures-1] += comp_time+comm_time;			// update local runtime
  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
  volume_costs[num_volume_measures-5] = volume_costs[num_volume_measures-5] > critical_path_costs[num_critical_path_measures-5]
                                          ? critical_path_costs[num_critical_path_measures-5] : volume_costs[num_volume_measures-5];
  volume_costs[num_volume_measures-4] = volume_costs[num_volume_measures-4] > critical_path_costs[num_critical_path_measures-4]
                                          ? critical_path_costs[num_critical_path_measures-4] : volume_costs[num_volume_measures-4];
  volume_costs[num_volume_measures-3] = volume_costs[num_volume_measures-3] > critical_path_costs[num_critical_path_measures-3]
                                          ? critical_path_costs[num_critical_path_measures-3] : volume_costs[num_volume_measures-3];
  volume_costs[num_volume_measures-2] = volume_costs[num_volume_measures-2] > critical_path_costs[num_critical_path_measures-2]
                                          ? critical_path_costs[num_critical_path_measures-2] : volume_costs[num_volume_measures-2];
  volume_costs[num_volume_measures-1] = volume_costs[num_volume_measures-1] > critical_path_costs[num_critical_path_measures-1]
                                          ? critical_path_costs[num_critical_path_measures-1] : volume_costs[num_volume_measures-1];

  // Decompose measurements along multiple paths by MPI routine.
  // Accumuate MPI routine-local measurements. The "my_..." members will never modify the accumulations, while the "critical_path_..." will first accumulate before path propagation.
  *tracker.my_synch_time   += 0;			// Nonblocking routines will have no synchronization time component
  *tracker.my_comm_time    += comm_time;
  save=0;
  for (int j=0; j<cost_models.size(); j++){
    if (cost_models[j]=='1'){
      *(tracker.my_msg_count+save) += costs[j].first;
      *(tracker.my_wrd_count+save) += costs[j].second;
      save++;
    }
  }
  save=0;
  for (size_t i=0; i<comm_path_select_size; i++){
    *(tracker.critical_path_synch_time+i)   += tracker.synch_time;
    *(tracker.critical_path_comm_time+i)    += comm_time;
  }
  for (int j=0; j<cost_models.size(); j++){
    for (size_t i=0; i<comm_path_select_size; i++){
      if (cost_models[j]=='1'){
        *(tracker.critical_path_msg_count+save*comm_path_select_size+i) += costs[j].first;
        *(tracker.critical_path_wrd_count+save*comm_path_select_size+i) += costs[j].second;
      }
    }
    save++;
  }

  // Decompose measurements along multiple paths by symbol
  if (symbol_path_select_size>0 && symbol_stack.size()>0){
    for (auto i=0; i<symbol_path_select_size; i++){
      // update all communication-related measures for the top symbol in stack
      size_t save=0;
      for (int j=0; j<cost_models.size(); j++){
        if (cost_models[j]=='1'){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][save] += costs[j].second;
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][cost_models.size()+save] += costs[j].first;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][save] += costs[j].second;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][cost_models.size()+save] += costs[j].first;
        }
        save++;
      }
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-5] += comm_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-4] += 0.;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-3] += comp_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-1] += (comp_time+comm_time);
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-5] += comm_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-4] += 0.;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-3] += comp_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-1] += (comp_time+comm_time);
    }
    // update all communication-related measures for the top symbol in stack
    size_t save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]=='1'){
        symbol_timers[symbol_stack.top()].pp_exclusive_measure[save] += costs[j].second;
        symbol_timers[symbol_stack.top()].pp_exclusive_measure[cost_models.size()+save] += costs[j].first;
        symbol_timers[symbol_stack.top()].pp_excl_measure[save] += costs[j].second;
        symbol_timers[symbol_stack.top()].pp_excl_measure[cost_models.size()+save] += costs[j].first;
      }
      save++;
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-5] += comm_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-4] += 0.;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-3] += comp_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += (comp_time+comm_time);
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-5] += comm_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-4] += 0.;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-3] += comp_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += (comp_time+comm_time);
  }

  nonblocking_internal_info.erase(*request);
  tracker.start_time = MPI_Wtime();
}

int path::complete_comm(float curtime, MPI_Request* request, MPI_Status* status){
  float comp_time = curtime - computation_timer;
  int ret = MPI_SUCCESS;
  auto info_it = nonblocking_internal_info.find(*request);
  assert(info_it != nonblocking_internal_info.end());
  MPI_Request save_request = info_it->first;
  volatile float last_start_time = MPI_Wtime();
  ret = PMPI_Wait(request, status);
  float save_comm_time = MPI_Wtime() - last_start_time;
  complete_path_update();
  if (info_it->second.partner == MPI_ANY_SOURCE) { info_it->second.track->partner1 = status->MPI_SOURCE; }
  complete_comm(*info_it->second.track, &save_request, comp_time, save_comm_time);
  computation_timer = MPI_Wtime();
  if (symbol_path_select_size>0){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
  return ret;
}

int path::complete_comm(float curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){

  float waitany_comp_time = curtime - computation_timer;
  int ret = MPI_SUCCESS;
  // We must force 'track_p2p_idle' to be zero because we don't know which request the MPI implementation will choose before
  //   it chooses it. Note that this is a not a problem for MPI_Waitall because all requests are chosen.
  //   Thus, we cannot participate in any idle/synch time exchanges. This is not a big deal at all if MPI_Waitany
  //   is used for both sender and receiver side, as idle/synch time are zero anyways. This decision is limiting in the sense that
  //   if MPI_Waitany is used to close nonblocking requests on both the sender and receiver side, then we forfeit tracking of idle/synch
  //   time caused by blocking send/recv.
  // We must save the requests before the completition of a request by the MPI implementation because its tag is set to MPI_REQUEST_NULL and lost forever
  std::vector<MPI_Request> pt(count); for (int i=0;i<count;i++){pt[i]=(array_of_requests)[i];}
  volatile float last_start_time = MPI_Wtime();
  ret = PMPI_Waitany(count,array_of_requests,indx,status);
  float waitany_comm_time = MPI_Wtime() - last_start_time;
  complete_path_update();
  MPI_Request request = pt[*indx];
  auto info_it = nonblocking_internal_info.find(request);
  assert(info_it != nonblocking_internal_info.end());
  if (info_it->second.partner == MPI_ANY_SOURCE) { info_it->second.track->partner1 = status->MPI_SOURCE; }
  complete_comm(*info_it->second.track, &request, waitany_comp_time, waitany_comm_time);
  computation_timer = MPI_Wtime();
  if (symbol_path_select_size>0){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
  return ret;
}

int path::complete_comm(float curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[],
                        MPI_Status array_of_statuses[]){

  float waitsome_comp_time = curtime - computation_timer;
  int ret = MPI_SUCCESS;
  is_first_request=true;
  // We must save the requests before the completition of a request by the MPI implementation because its tag is set to MPI_REQUEST_NULL and lost forever
  std::vector<MPI_Request> pt(incount); for (int i=0;i<incount;i++){pt[i]=(array_of_requests)[i];}
  volatile float last_start_time = MPI_Wtime();
  ret = PMPI_Waitsome(incount,array_of_requests,outcount,array_of_indices,array_of_statuses);
  float waitsome_comm_time = MPI_Wtime() - last_start_time;
  complete_path_update();
  for (int i=0; i<*outcount; i++){
    MPI_Request request = pt[(array_of_indices)[i]];
    auto info_it = nonblocking_internal_info.find(request);
    assert(info_it != nonblocking_internal_info.end());
    if (info_it->second.partner == MPI_ANY_SOURCE) { info_it->second.track->partner1 = (array_of_statuses)[i].MPI_SOURCE; }
    complete_comm(*info_it->second.track, &request, waitsome_comp_time, waitsome_comm_time);
    waitsome_comp_time=0; waitsome_comm_time=0;
    if (i==0){is_first_request=false;}
  }
  computation_timer = MPI_Wtime();
  if (symbol_path_select_size>0){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
  return ret;
}

int path::complete_comm(float curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
  float waitall_comp_time = curtime - computation_timer;
  int ret = MPI_SUCCESS;
  is_first_request=true;
  std::vector<MPI_Request> internal_requests(3*count,MPI_REQUEST_NULL);
  // Issue all barrier/synch communications at once because request order is not guaranteed to be sequenced together on all processes.
  // Necessary to avoid corruption of idle time calculation that would occur if sending out in some sequence after each request is completed.
  // Presumably the sending communications will utilize the eager protocol, but as the Waitall is issued immediately following the loop, its irrelevant.
  // TODO: Staging nonblocking receives with nonblocking sends might be a problem because the Sends will likely use the eager protocol (since 1-byte messages),
  //         while the Recvs will block until the message buffer is available. This will force the matching blocking receives of the nonblocking senders to possibely start counting communication time while
  //         the sender is stuck in the Waitall loop with the other sends and receives, thus corrupting the measurement of communication time. I propose that the Sends issue their 3-messages together followed
  //         by a Waitall. At this point, the choice is whether to issue the user-communication sends one-by-one via Waitany loop, issue the nonblocking recv 3-messages together followed by waitall,
  //         or actually issue a Waitall for the user communication rather than a loop over Waits. I am now leaning towards supporting the user-communication Waitall instead, as that might be a source of overhead
  //         with CTF. Each process can utilize its timer and record the same communication time, and then issue the exchange of path information via nonblocking communications.
  // We must save the requests before the completition of a request by the MPI implementation because its tag is set to MPI_REQUEST_NULL and lost forever
  std::vector<MPI_Request> pt(count); for (int i=0;i<count;i++){pt[i]=(array_of_requests)[i];}
  volatile float last_start_time = MPI_Wtime();
  ret = PMPI_Waitall(count,array_of_requests,array_of_statuses);
  float waitall_comm_time = MPI_Wtime() - last_start_time;
  complete_path_update();
  for (int i=0; i<count; i++){
    MPI_Request request = pt[i];
    auto info_it = nonblocking_internal_info.find(request);
    assert(info_it != nonblocking_internal_info.end());
    if (info_it->second.partner == MPI_ANY_SOURCE) { info_it->second.track->partner1 = (array_of_statuses)[i].MPI_SOURCE; }
    complete_comm(*info_it->second.track, &request, waitall_comp_time, waitall_comm_time);
    // Although we have to exchange the path data for each request, we do not want to float-count the computation time nor the communicaion time
    waitall_comp_time=0; waitall_comm_time=0;
    if (i==0){is_first_request=false;}
  }
  computation_timer = MPI_Wtime();
  if (symbol_path_select_size>0){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
  return ret;
}

// *****************************************************************************************************************************************************
void path::propagate_symbols(nonblocking& tracker, int rank){
  if (tracker.is_sender){
    MPI_Request internal_request[4];
    int* send_envelope1 = nullptr; int* send_envelope2 = nullptr; float* send_envelope3 = nullptr; char* send_envelope5 = nullptr;
    int ftimer_size = symbol_timers.size();
    int data_len_size = symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*ftimer_size;
    int num_chars = 0;
    for (int i=0; i<ftimer_size; i++) { num_chars += symbol_order[i].size(); }
    send_envelope1 = (int*)malloc(sizeof(int)); *send_envelope1 = ftimer_size;
    send_envelope2 = (int*)malloc(sizeof(int)*(ftimer_size));
    send_envelope3 = (float*)malloc(sizeof(float)*data_len_size);
    send_envelope5 = (char*)malloc(sizeof(char)*num_chars);
    int symbol_offset = 0;
    for (auto i=0; i<ftimer_size; i++){
      send_envelope2[i] = symbol_order[i].size();
      for (auto j=0; j<symbol_order[i].size(); j++){
        send_envelope5[symbol_offset+j] = symbol_order[i][j];
      }
      symbol_offset += symbol_order[i].size();
    }
    std::memcpy(send_envelope3,&symbol_timer_pad_local_cp[0],sizeof(float)*data_len_size);
    PMPI_Isend(&send_envelope1[0],1,MPI_INT,tracker.partner1,internal_tag1,tracker.comm,&internal_request[0]);
    PMPI_Isend(&send_envelope2[0],ftimer_size,MPI_INT,tracker.partner1,internal_tag2,tracker.comm,&internal_request[1]);
    PMPI_Isend(&send_envelope3[0],data_len_size,MPI_FLOAT,tracker.partner1,internal_tag3,tracker.comm,&internal_request[2]);
    PMPI_Isend(&send_envelope5[0],symbol_offset,MPI_CHAR,tracker.partner1,internal_tag5,tracker.comm,&internal_request[3]);

    for (int i=0; i<4; i++) { internal_timer_prop_req.push_back(internal_request[i]); }
    internal_timer_prop_int.push_back(send_envelope1);
    internal_timer_prop_int.push_back(send_envelope2);
    internal_timer_prop_float.push_back(send_envelope3);
    internal_timer_prop_char.push_back(send_envelope5);
  } else{
    MPI_Request internal_request[4];
    int* recv_envelope1 = nullptr; int* recv_envelope2 = nullptr; float* recv_envelope3 = nullptr; char* recv_envelope5 = nullptr;
    recv_envelope1 = (int*)malloc(sizeof(int));
    recv_envelope2 = (int*)malloc(sizeof(int)*(max_num_symbols));
    recv_envelope3 = (float*)malloc(sizeof(float)*symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*max_num_symbols);
    recv_envelope5 = (char*)malloc(sizeof(char)*max_timer_name_length*max_num_symbols);
    PMPI_Irecv(recv_envelope1,1,MPI_INT,tracker.partner1,internal_tag1,tracker.comm,&internal_request[0]);
    PMPI_Irecv(recv_envelope2,max_num_symbols,MPI_INT,tracker.partner1,internal_tag2,tracker.comm,&internal_request[1]);
    PMPI_Irecv(recv_envelope3,symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*max_num_symbols,MPI_FLOAT,tracker.partner1,internal_tag3,tracker.comm,&internal_request[2]);
    PMPI_Irecv(recv_envelope5,max_timer_name_length*max_num_symbols,MPI_CHAR,tracker.partner1,internal_tag5,tracker.comm,&internal_request[3]);

    for (int i=0; i<4; i++) { internal_timer_prop_req.push_back(internal_request[i]); }
    internal_timer_prop_int.push_back(recv_envelope1);
    internal_timer_prop_int.push_back(recv_envelope2);
    internal_timer_prop_float.push_back(recv_envelope3);
    internal_timer_prop_char.push_back(recv_envelope5);
  }
}


/*
 Its important to note here that a blocking p2p call will already know whether its the cp root or not, regardless of whether its partner used a nonblocking p2p routine.
   But, because that potential nonblocking partner does not have this knowledge, and thus posted both sends and recvs, the blocking partner also has to do so as well, even if its partner (unknown to him) used a blocking p2p routine.
*/
void path::propagate_symbols(blocking& tracker, int rank){
  bool eager = ((tracker.partner1 != -1) && (tracker.tag!=13) && (tracker.tag!=14) && (tracker.tag!=15) && (tracker.nbytes <= eager_limit));
  std::vector<int> ftimer_size_cp(symbol_path_select_size,0);
  int ftimer_size_ncp1=0;
  int ftimer_size_ncp2=0;
  for (auto i=0; i<symbol_path_select_size; i++){
    ftimer_size_cp[i] = symbol_timers.size();
    // Specialized multi-root reductions for propagating user collectives require non-roots to set values to 0.
    // This is not the case for p2p, so we must separate the logic.
    if ((tracker.partner1 == -1) && (rank != info_receiver[symbol_path_select_index[i]].second)){
      ftimer_size_cp[i] = 0;
    }
  }

  if (tracker.partner1 == -1){
    PMPI_Allreduce(MPI_IN_PLACE,&ftimer_size_cp[0],symbol_path_select_size,MPI_INT,MPI_SUM,tracker.comm);
    memset(&symbol_len_pad_cp[0],0,sizeof(int)*symbol_len_pad_cp.size());// not as simple as 'ftimer_size_cp' for blocking collectives. Dependent on the entries in that array
    size_t symbol_offset_cp = 0;
    for (auto k=0; k<symbol_path_select_size; k++){
      // Only the roots determining each path will write the symbol length for its symbols.
      //   The rest must keep the zero set above in the memset, but they will still increment the symbol offset counter.
      if (rank==info_receiver[symbol_path_select_index[k]].second){
        for (auto i=0; i<ftimer_size_cp[k]; i++){
          symbol_len_pad_cp[symbol_offset_cp++] = symbol_order[i].size();
        }
      }
      else{
        symbol_offset_cp += ftimer_size_cp[k];
      }
    }
    PMPI_Allreduce(MPI_IN_PLACE,&symbol_len_pad_cp[0],symbol_offset_cp,MPI_INT,MPI_SUM,tracker.comm);
    symbol_offset_cp = 0;
    int char_count_cp = 0; int char_count_ncp1 = 0; int char_count_ncp2 = 0;
    size_t pad_global_offset = 0;
    for (auto k=0; k<symbol_path_select_size; k++){
      // Only the roots determining each path will write the symbol length for its symbols.
      //   The rest must keep the zero set above in the memset, but they will still increment the symbol offset counter.
      if (rank==info_receiver[symbol_path_select_index[k]].second){
        for (auto i=0; i<ftimer_size_cp[k]; i++){
          size_t pad_local_offset = (i*symbol_path_select_size+k)*(cp_symbol_class_count*num_per_process_measures+1);
          std::memcpy(&symbol_timer_pad_global_cp[pad_global_offset],
                      &symbol_timer_pad_local_cp[pad_local_offset],
                      sizeof(float)*(cp_symbol_class_count*num_per_process_measures+1));
          for (auto j=0; j<symbol_len_pad_cp[symbol_offset_cp]; j++){
            symbol_pad_cp[char_count_cp+j] = symbol_order[i][j];
          }
          pad_global_offset += (cp_symbol_class_count*num_per_process_measures+1);
          char_count_cp += symbol_len_pad_cp[symbol_offset_cp++];
        }
      }
      else{
        memset(&symbol_timer_pad_global_cp[pad_global_offset],0,sizeof(float)*ftimer_size_cp[k]*(cp_symbol_class_count*num_per_process_measures+1));
        pad_global_offset += ftimer_size_cp[k]*(cp_symbol_class_count*num_per_process_measures+1);
        for (auto i=0; i<ftimer_size_cp[k]; i++){
          memset(&symbol_pad_cp[char_count_cp],0,sizeof(char)*symbol_len_pad_cp[symbol_offset_cp]);
          char_count_cp += symbol_len_pad_cp[symbol_offset_cp++];
        }
      }
    }
    PMPI_Allreduce(MPI_IN_PLACE,&symbol_timer_pad_global_cp[0],pad_global_offset,MPI_FLOAT,MPI_SUM,tracker.comm);
    PMPI_Allreduce(MPI_IN_PLACE,&symbol_pad_cp[0],char_count_cp,MPI_CHAR,MPI_SUM,tracker.comm);
    pad_global_offset = 0;
    size_t symbol_pad_offset = 0;
    size_t symbol_len_pad_offset=0;
    for (auto k=0; k<symbol_path_select_size; k++){
      if (rank != info_receiver[symbol_path_select_index[k]].second){
        for (int i=0; i<ftimer_size_cp[k]; i++){
          auto reconstructed_symbol = std::string(symbol_pad_cp.begin()+symbol_pad_offset,symbol_pad_cp.begin()+symbol_pad_offset+symbol_len_pad_cp[symbol_len_pad_offset]);
          if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
            symbol_timers[reconstructed_symbol] = symbol_tracker(reconstructed_symbol);
            symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
          }
          std::memcpy(symbol_timers[reconstructed_symbol].cp_numcalls[k],
                      &symbol_timer_pad_global_cp[pad_global_offset],
                      sizeof(float)*(cp_symbol_class_count*num_per_process_measures+1));
          symbol_timers[reconstructed_symbol].has_been_processed = true;
          pad_global_offset += (cp_symbol_class_count*num_per_process_measures+1);
          symbol_pad_offset += symbol_len_pad_cp[symbol_len_pad_offset];
          symbol_len_pad_offset++;
        }
        // Now cycle through and find the symbols that were not processed and set their accumulated measures to 0
        for (auto& it : symbol_timers){
          if (it.second.has_been_processed){ it.second.has_been_processed = false; }
          else{
            it.second.cp_numcalls[k][0] = 0;
            for (int j=0; j<num_per_process_measures; j++){
              it.second.cp_incl_measure[k][j] = 0;
              it.second.cp_excl_measure[k][j] = 0;
            }
          }
        }
      }
      else{
        for (int i=0; i<ftimer_size_cp[k]; i++){
          symbol_pad_offset += symbol_len_pad_cp[symbol_len_pad_offset++];
        }
        pad_global_offset += ftimer_size_cp[k]*(cp_symbol_class_count*num_per_process_measures+1);   
      }
    }
  }
  else if (!eager){
    // This propagation for p2p user communication is agnostic (for now) to which process determines the root for a specific metric.
    //   Its simply an exchange. Note that this allows each process to only send the bare minimum. As an example, each process need only send
    //     a single integer representing its symbol size. There is no need to send symbol_path_select_size integers with the same value.
    MPI_Request symbol_exchance_reqs[8]; int exchange_count=0;
    PMPI_Isend(&ftimer_size_cp[0],1,MPI_INT,tracker.partner1,internal_tag1,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
    PMPI_Irecv(&ftimer_size_ncp1,1,MPI_INT,tracker.partner1,internal_tag1,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
    if (tracker.partner1 != tracker.partner2){ PMPI_Isend(&ftimer_size_cp[0],1,MPI_INT,tracker.partner2,internal_tag1,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
                                               PMPI_Irecv(&ftimer_size_ncp2,1,MPI_INT,tracker.partner2,internal_tag1,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
                                             }
    PMPI_Waitall(exchange_count,&symbol_exchance_reqs[0],MPI_STATUSES_IGNORE);
    memset(&symbol_len_pad_cp[0],0,sizeof(int)*symbol_len_pad_cp.size());// not as simple as 'ftimer_size_cp' for blocking collectives. Dependent on the entries in that array
    memset(&symbol_len_pad_ncp1[0],0,sizeof(int)*ftimer_size_ncp1);
    memset(&symbol_len_pad_ncp2[0],0,sizeof(int)*ftimer_size_ncp2);
    size_t symbol_offset_cp = 0; size_t symbol_offset_ncp1 = 0; size_t symbol_offset_ncp2 = 0;
    // Each process will determine the symbol length for each of its symbols first
    //   while incrementing simply the counters to prepare to receive.
    for (auto i=0; i<ftimer_size_cp[0]; i++){
      symbol_len_pad_cp[i] = symbol_order[i].size();
    }
    // This propagation for p2p user communication is agnostic (for now) to which process determines the root for a specific metric.
    //   Its simply an exchange.
    exchange_count=0;
    PMPI_Isend(&symbol_len_pad_cp[0],ftimer_size_cp[0],MPI_INT,tracker.partner1,internal_tag2,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
    PMPI_Irecv(&symbol_len_pad_ncp1[0],ftimer_size_ncp1,MPI_INT,tracker.partner1,internal_tag2,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
    if (tracker.partner1 != tracker.partner2){ PMPI_Isend(&symbol_len_pad_cp[0],ftimer_size_cp[0],MPI_INT,tracker.partner2,internal_tag2,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
                                               PMPI_Irecv(&symbol_len_pad_ncp2[0],ftimer_size_ncp2,MPI_INT,tracker.partner2,internal_tag2,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
                                             }
    PMPI_Waitall(exchange_count,&symbol_exchance_reqs[0],MPI_STATUSES_IGNORE);
    symbol_offset_cp = 0; symbol_offset_ncp1 = 0; symbol_offset_ncp2 = 0;
    int char_count_cp = 0; int char_count_ncp1 = 0; int char_count_ncp2 = 0;
    size_t pad_global_offset = 0;
    // Each process will determine the symbol length for each of its symbols first
    //   while incrementing simply the counters to prepare to receive.
    for (auto i=0; i<ftimer_size_cp[0]; i++){
      for (auto j=0; j<symbol_len_pad_cp[symbol_offset_cp]; j++){
        symbol_pad_cp[char_count_cp+j] = symbol_order[i][j];
      }
      char_count_cp += symbol_len_pad_cp[symbol_offset_cp++];
    }
    for (auto i=0; i<ftimer_size_ncp1; i++){
      char_count_ncp1 += symbol_len_pad_ncp1[symbol_offset_ncp1++];
    }
    if (tracker.partner1 != tracker.partner2){
      for (auto i=0; i<ftimer_size_ncp2; i++){
        char_count_ncp2 += symbol_len_pad_ncp2[symbol_offset_ncp2++];
      }
    }
    // This propagation for p2p user communication is agnostic (for now) to which process determines the root for a specific metric.
    //   Its simply an exchange. No special copying is needed as in collectives case.
    exchange_count=0;
    int data_len_cp = symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*ftimer_size_cp[0];
    int data_len_ncp1 = symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*ftimer_size_ncp1;
    int data_len_ncp2 = symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*ftimer_size_ncp2;
    PMPI_Isend(&symbol_timer_pad_local_cp[0],data_len_cp,MPI_FLOAT,tracker.partner1,internal_tag3,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
    PMPI_Isend(&symbol_pad_cp[0],char_count_cp,MPI_CHAR,tracker.partner1,internal_tag5,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
    PMPI_Irecv(&symbol_timer_pad_global_cp[0],data_len_ncp1,MPI_FLOAT,tracker.partner1,internal_tag3,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
    PMPI_Irecv(&symbol_pad_ncp1[0],char_count_ncp1,MPI_CHAR,tracker.partner1,internal_tag5,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
    if (tracker.partner1 != tracker.partner2){
      PMPI_Isend(&symbol_timer_pad_local_cp[0],data_len_cp,MPI_FLOAT,tracker.partner2,internal_tag3,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
      PMPI_Isend(&symbol_pad_cp[0],char_count_cp,MPI_CHAR,tracker.partner2,internal_tag5,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
      PMPI_Irecv(&symbol_timer_pad_global_cp2[0],data_len_ncp2,MPI_FLOAT,tracker.partner2,internal_tag3,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
      PMPI_Irecv(&symbol_pad_ncp2[0],char_count_ncp2,MPI_CHAR,tracker.partner2,internal_tag5,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
    }
    PMPI_Waitall(exchange_count,&symbol_exchance_reqs[0],MPI_STATUSES_IGNORE);
    for (auto k=0; k<symbol_path_select_size; k++){
      bool foreign_root = true;
      if (rank == info_receiver[symbol_path_select_index[k]].second){
        foreign_root=false;
      }
      else if (tracker.partner1 == info_receiver[symbol_path_select_index[k]].second){
        size_t symbol_pad_offset = 0;
        size_t symbol_len_pad_offset=0;
        std::string reconstructed_symbol;
        for (int i=0; i<ftimer_size_ncp1; i++){
          reconstructed_symbol = std::string(symbol_pad_ncp1.begin()+symbol_pad_offset,symbol_pad_ncp1.begin()+symbol_pad_offset+symbol_len_pad_ncp1[symbol_len_pad_offset]);
          if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
            symbol_timers[reconstructed_symbol] = symbol_tracker(reconstructed_symbol);
            symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
          }
          std::memcpy(symbol_timers[reconstructed_symbol].cp_numcalls[k],
                      &symbol_timer_pad_global_cp[(i*symbol_path_select_size+k)*(cp_symbol_class_count*num_per_process_measures+1)],
                      sizeof(float)*(cp_symbol_class_count*num_per_process_measures+1));
          symbol_timers[reconstructed_symbol].has_been_processed = true;
          symbol_pad_offset += symbol_len_pad_ncp1[symbol_len_pad_offset++];
        }
      }
      else{
        size_t symbol_pad_offset = 0;
        size_t symbol_len_pad_offset=0;
        std::string reconstructed_symbol;
        for (int i=0; i<ftimer_size_ncp2; i++){
          reconstructed_symbol = std::string(symbol_pad_ncp2.begin()+symbol_pad_offset,symbol_pad_ncp2.begin()+symbol_pad_offset+symbol_len_pad_ncp2[symbol_len_pad_offset]);
          if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
            symbol_timers[reconstructed_symbol] = symbol_tracker(reconstructed_symbol);
            symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
          }
          std::memcpy(&symbol_timers[reconstructed_symbol].cp_numcalls[k],
                      &symbol_timer_pad_global_cp2[(i*symbol_path_select_size+k)*(cp_symbol_class_count*num_per_process_measures+1)],
                      sizeof(float)*(cp_symbol_class_count*num_per_process_measures+1));
          symbol_timers[reconstructed_symbol].has_been_processed = true;
          symbol_pad_offset += symbol_len_pad_ncp2[symbol_len_pad_offset];
          symbol_len_pad_offset++;
        }
      }
      if (foreign_root){
        // Now cycle through and find the symbols that were not processed and set their accumulated measures to 0
        for (auto& it : symbol_timers){
          if (it.second.has_been_processed){ it.second.has_been_processed = false; }
          else{
            it.second.cp_numcalls[k][0] = 0;
            for (int j=0; j<num_per_process_measures; j++){
              it.second.cp_incl_measure[k][j] = 0;
              it.second.cp_excl_measure[k][j] = 0;
            }
          }
        }
      }
    }
  }
  else{
    // This propagation for p2p user communication is agnostic (for now) to which process determines the root for a specific metric.
    //   Its simply an exchange. Note that this allows each process to only send the bare minimum. As an example, each process need only send
    //     a single integer representing its symbol size. There is no need to send symbol_path_select_size integers with the same value.
    if (tracker.is_sender){
      PMPI_Bsend(&ftimer_size_cp[0],1,MPI_INT,tracker.partner1,internal_tag1,tracker.comm);
      memset(&symbol_len_pad_cp[0],0,sizeof(int)*symbol_len_pad_cp.size());// not as simple as 'ftimer_size_cp' for blocking collectives. Dependent on the entries in that array
      // Each process will determine the symbol length for each of its symbols first
      //   while incrementing simply the counters to prepare to receive.
      for (auto i=0; i<ftimer_size_cp[0]; i++){
        symbol_len_pad_cp[i] = symbol_order[i].size();
      }
      PMPI_Bsend(&symbol_len_pad_cp[0],ftimer_size_cp[0],MPI_INT,tracker.partner1,internal_tag2,tracker.comm);
      int char_count_cp = 0;
      size_t pad_global_offset = 0;
      // Each process will determine the symbol length for each of its symbols first
      //   while incrementing simply the counters to prepare to receive.
      for (auto i=0; i<ftimer_size_cp[0]; i++){
        for (auto j=0; j<symbol_len_pad_cp[i]; j++){
          symbol_pad_cp[char_count_cp+j] = symbol_order[i][j];
        }
        char_count_cp += symbol_len_pad_cp[i];
      }
      int data_len_cp = symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*ftimer_size_cp[0];
      PMPI_Bsend(&symbol_timer_pad_local_cp[0],data_len_cp,MPI_FLOAT,tracker.partner1,internal_tag3,tracker.comm);
      PMPI_Bsend(&symbol_pad_cp[0],char_count_cp,MPI_CHAR,tracker.partner1,internal_tag5,tracker.comm);
    } else{
      PMPI_Recv(&ftimer_size_cp[0],1,MPI_INT,tracker.partner1,internal_tag1,tracker.comm,MPI_STATUS_IGNORE);
      memset(&symbol_len_pad_cp[0],0,sizeof(int)*symbol_len_pad_cp.size());// not as simple as 'ftimer_size_cp' for blocking collectives. Dependent on the entries in that array
      PMPI_Recv(&symbol_len_pad_cp[0],ftimer_size_cp[0],MPI_INT,tracker.partner1,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      int char_count_cp = 0;
      size_t pad_global_offset = 0;
      for (auto i=0; i<ftimer_size_cp[0]; i++){
        char_count_cp += symbol_len_pad_cp[i];
      }
      int data_len_cp = symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*ftimer_size_cp[0];
      PMPI_Recv(&symbol_timer_pad_global_cp[0],data_len_cp,MPI_FLOAT,tracker.partner1,internal_tag3,tracker.comm,MPI_STATUS_IGNORE);
      PMPI_Recv(&symbol_pad_cp[0],char_count_cp,MPI_CHAR,tracker.partner1,internal_tag5,tracker.comm,MPI_STATUS_IGNORE);
    }
    for (auto k=0; k<symbol_path_select_size; k++){
      if (info_sender[symbol_path_select_index[k]].second < info_receiver[symbol_path_select_index[k]].second){
        size_t symbol_pad_offset = 0;
        size_t symbol_len_pad_offset=0;
        std::string reconstructed_symbol;
        for (int i=0; i<ftimer_size_cp[0]; i++){
          reconstructed_symbol = std::string(symbol_pad_cp.begin()+symbol_pad_offset,symbol_pad_cp.begin()+symbol_pad_offset+symbol_len_pad_cp[symbol_len_pad_offset]);
          if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
            symbol_timers[reconstructed_symbol] = symbol_tracker(reconstructed_symbol);
            symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
          }
          std::memcpy(symbol_timers[reconstructed_symbol].cp_numcalls[k],
                      &symbol_timer_pad_global_cp[(i*symbol_path_select_size+k)*(cp_symbol_class_count*num_per_process_measures+1)],
                      sizeof(float)*(cp_symbol_class_count*num_per_process_measures+1));
          symbol_timers[reconstructed_symbol].has_been_processed = true;
          symbol_pad_offset += symbol_len_pad_cp[symbol_len_pad_offset++];
        }
      }
      // Now cycle through and find the symbols that were not processed and set their accumulated measures to 0
      for (auto& it : symbol_timers){
        if (it.second.has_been_processed){ it.second.has_been_processed = false; }
        else{
          it.second.cp_numcalls[k][0] = 0;
          for (int j=0; j<num_per_process_measures; j++){
            it.second.cp_incl_measure[k][j] = 0;
            it.second.cp_excl_measure[k][j] = 0;
          }
        }
      }
    }
  }
}

void path::propagate(blocking& tracker){
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  if ((rank == tracker.partner1) && (rank == tracker.partner2)) { return; } 
  bool eager = ((tracker.partner1 != -1) && (tracker.tag!=13) && (tracker.tag!=14) && (tracker.tag!=15) && (tracker.nbytes <= eager_limit));
  if (symbol_path_select_size>0){
    for (int i=0; i<num_critical_path_measures; i++){
      info_sender[i].first = critical_path_costs[i];
      info_sender[i].second = rank;
    }
    if (tracker.partner1 == -1){
      PMPI_Allreduce(&info_sender[0].first, &info_receiver[0].first, num_critical_path_measures, MPI_FLOAT_INT, MPI_MAXLOC, tracker.comm);
    }
    else{
      if (!eager && send_dependency){
        PMPI_Sendrecv(&info_sender[0].first, num_critical_path_measures, MPI_FLOAT_INT, tracker.partner1, internal_tag,
                      &info_receiver[0].first, num_critical_path_measures, MPI_FLOAT_INT, tracker.partner2, internal_tag, tracker.comm, MPI_STATUS_IGNORE);
        for (int i=0; i<num_critical_path_measures; i++){
          if (info_sender[i].first>info_receiver[i].first){info_receiver[i].second = rank;}
          else if (info_sender[i].first==info_receiver[i].first){ info_receiver[i].second = std::min(rank,tracker.partner1); }
          info_receiver[i].first = std::max(info_sender[i].first, info_receiver[i].first);
        }
        if (tracker.partner2 != tracker.partner1){
          // Assuming the sender is always dependent on the receiver (not necessarily true or eager protocol, but we make this assumption), this condition signifies a 3-process exchange.
          for (int i=0; i<num_critical_path_measures; i++){
            info_sender[i].first = info_receiver[i].first;
            info_sender[i].second = info_receiver[i].second;
          }
          PMPI_Sendrecv(&info_sender[0].first, num_critical_path_measures, MPI_FLOAT_INT, tracker.partner2, internal_tag,
                        &info_receiver[0].first, num_critical_path_measures, MPI_FLOAT_INT, tracker.partner1, internal_tag, tracker.comm, MPI_STATUS_IGNORE);
          for (int i=0; i<num_critical_path_measures; i++){
            if (info_sender[i].first>info_receiver[i].first){info_receiver[i].second = rank;}
            else if (info_sender[i].first==info_receiver[i].first){ info_receiver[i].second = std::min(rank,tracker.partner1); }
            info_receiver[i].first = std::max(info_sender[i].first, info_receiver[i].first);
          }
        }
      }
      else{
        if (tracker.is_sender){
          PMPI_Bsend(&info_sender[0].first, num_critical_path_measures, MPI_FLOAT_INT, tracker.partner1, internal_tag, tracker.comm);
        } else{
          PMPI_Recv(&info_receiver[0].first, num_critical_path_measures, MPI_FLOAT_INT, tracker.partner1, internal_tag, tracker.comm, MPI_STATUS_IGNORE);
        }
      }
    }
  }
  // Exchange the tracked routine critical path data
  if (tracker.partner1 == -1){
    MPI_Op op; MPI_Op_create((MPI_User_function*) propagate_critical_path_op,0,&op);
    PMPI_Allreduce(MPI_IN_PLACE, &critical_path_costs[0], critical_path_costs.size(), MPI_FLOAT, op, tracker.comm);
    MPI_Op_free(&op);
  }
  else{
    // Note that a blocking sendrecv allows exchanges even when the other party issued a request via nonblocking communication, as the process with the nonblocking request posts both sends and receives.
    if (!eager && send_dependency){
      PMPI_Sendrecv(&critical_path_costs[0], critical_path_costs.size(), MPI_FLOAT, tracker.partner1, internal_tag2, &new_cs[0], critical_path_costs.size(),
                    MPI_FLOAT, tracker.partner2, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
      update_critical_path(&new_cs[0],&critical_path_costs[0],critical_path_costs_size);
      if (tracker.partner2 != tracker.partner1){
        PMPI_Sendrecv(&critical_path_costs[0], critical_path_costs.size(), MPI_FLOAT, tracker.partner2, internal_tag2, &new_cs[0], critical_path_costs.size(), MPI_FLOAT, tracker.partner1, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
        update_critical_path(&new_cs[0],&critical_path_costs[0],critical_path_costs_size);
      }
    }
    else{
      if (tracker.is_sender){
        PMPI_Bsend(&critical_path_costs[0], critical_path_costs.size(), MPI_FLOAT, tracker.partner1, internal_tag2, tracker.comm);
      } else{
        PMPI_Recv(&new_cs[0], critical_path_costs.size(), MPI_FLOAT, tracker.partner1, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
        update_critical_path(&new_cs[0],&critical_path_costs[0],critical_path_costs_size);
      }
    }
  }
  if (symbol_path_select_size>0) { propagate_symbols(tracker,rank); }
}

void path::propagate(nonblocking& tracker){
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  if (rank == tracker.partner1) { return; } 
  if (symbol_path_select_size>0){
    if (tracker.partner1 == -1){ assert(0); }
    for (int i=0; i<num_critical_path_measures; i++){
      info_sender[i].first = critical_path_costs[i];
      info_sender[i].second = rank;
    }
    if (tracker.is_sender){
      MPI_Request req1;
      float_int* send_pathdata = (float_int*)malloc(num_critical_path_measures*sizeof(float_int));
      memcpy(&send_pathdata[0].first, &info_sender[0].first, num_critical_path_measures*sizeof(float_int));
      PMPI_Isend(&send_pathdata[0].first, num_critical_path_measures, MPI_FLOAT_INT, tracker.partner1, internal_tag, tracker.comm, &req1);
      internal_timer_prop_req.push_back(req1);
      internal_timer_prop_float_int.push_back(send_pathdata);
    } else{
      MPI_Request req1;
      float_int* recv_pathdata = (float_int*)malloc(num_critical_path_measures*sizeof(float_int));
      PMPI_Irecv(&recv_pathdata[0].first, num_critical_path_measures, MPI_FLOAT_INT, tracker.partner1, internal_tag, tracker.comm, &req1);
      internal_timer_prop_req.push_back(req1);
      internal_timer_prop_float_int.push_back(recv_pathdata);
    }
  }
  // Exchange the tracked routine critical path data
  if (tracker.partner1 == -1){
    MPI_Op op; MPI_Op_create((MPI_User_function*) propagate_critical_path_op,0,&op);
    MPI_Request req1;
    float* local_path_data = (float*)malloc(critical_path_costs.size()*sizeof(float));
    std::memcpy(local_path_data, &critical_path_costs[0], critical_path_costs.size()*sizeof(float));
    PMPI_Iallreduce(MPI_IN_PLACE,local_path_data,critical_path_costs.size(),MPI_FLOAT,op,tracker.comm,&req1);
    //MPI_Op_free(&op);
    internal_comm_prop.push_back(std::make_pair(local_path_data,true));
    internal_comm_prop_req.push_back(req1);
  }
  else{
    MPI_Request req1;
    if (tracker.is_sender){
      float* local_path_data = (float*)malloc(critical_path_costs.size()*sizeof(float));
      std::memcpy(local_path_data, &critical_path_costs[0], critical_path_costs.size()*sizeof(float));
      PMPI_Isend(local_path_data, critical_path_costs.size(), MPI_FLOAT, tracker.partner1, internal_tag2, tracker.comm, &req1);
      internal_comm_prop.push_back(std::make_pair(local_path_data,true));
      internal_comm_prop_req.push_back(req1);
    }
    else{
      float* remote_path_data = (float*)malloc(critical_path_costs.size()*sizeof(float));
      PMPI_Irecv(remote_path_data, critical_path_costs.size(), MPI_FLOAT, tracker.partner1, internal_tag2, tracker.comm, &req1);
      internal_comm_prop.push_back(std::make_pair(remote_path_data,false));
      internal_comm_prop_req.push_back(req1);
    }
  }
  if (symbol_path_select_size>0) { propagate_symbols(tracker,rank); }
}

}
}
}
