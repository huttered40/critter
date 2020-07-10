#include "decomposition.h"
#include "../../container/symbol_tracker.h"
#include "../../util.h"

namespace critter{
namespace internal{

void decomposition::allocate(MPI_Comm comm){

  cp_symbol_class_count = 4;
  pp_symbol_class_count = 4;
  vol_symbol_class_count = 4;// should truly be 2, but set to 4 to conform to pp_symbol_class_count
  // The '2*comm_path_select_size' used below are used to track the computation time and idle time along each of the 'comm_path_select_size' paths.
  critical_path_costs_size            	= num_critical_path_measures+num_tracker_critical_path_measures*comm_path_select_size*list_size+2*comm_path_select_size;
  per_process_costs_size              	= num_per_process_measures+num_tracker_per_process_measures*comm_path_select_size*list_size+2*comm_path_select_size;
  volume_costs_size                   	= num_volume_measures+num_tracker_volume_measures*list_size;

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

static void add_critical_path_data_op(int_int_double* in, int_int_double* inout, int* len, MPI_Datatype* dtype){
  int_int_double* invec = in;
  int_int_double* inoutvec = inout;
  for (int i=0; i<*len; i++){
    inoutvec[i].first = std::max(inoutvec[i].first,invec[i].first);
    inoutvec[i].second = std::max(inoutvec[i].second,invec[i].second);
    inoutvec[i].third = std::max(inoutvec[i].third,invec[i].third);
  }
}

static void update_critical_path(double* in, double* inout, size_t len){
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

static void propagate_critical_path_op(double* in, double* inout, int* len, MPI_Datatype* dtype){
  update_critical_path(in,inout,static_cast<size_t>(*len));
}

static void complete_timers(double* remote_path_data, size_t msg_id){
  if (eager_p2p==0){
    int* envelope_int[2] = { internal_timer_prop_int[4*msg_id+2], internal_timer_prop_int[4*msg_id+3] };
    double* envelope_double[2] = { remote_path_data, internal_timer_prop_double[2*msg_id+1] };
    char* envelope_char = internal_timer_prop_char[2*msg_id+1];
    for (auto k=0; k<symbol_path_select_size; k++){
      // Up until this very point, we had no idea whether we, or our partner rank, determined the path for a specific metric.
      if (envelope_double[0][symbol_path_select_index[k]] > critical_path_costs[symbol_path_select_index[k]]){
        int ftimer_size = *envelope_int[0];
        int symbol_offset = 0;
        for (int i=0; i<ftimer_size; i++){
          auto reconstructed_symbol = std::string(envelope_char+symbol_offset,envelope_char+symbol_offset+envelope_int[1][i]);
          if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
            symbol_timers[reconstructed_symbol] = symbol_tracker(reconstructed_symbol);
            symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
          }
          std::memcpy(symbol_timers[reconstructed_symbol].cp_numcalls[k],
                      &envelope_double[1][(i*symbol_path_select_size+k)*(cp_symbol_class_count*num_per_process_measures+1)],
                      sizeof(double)*(cp_symbol_class_count*num_per_process_measures+1));
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
  } else{
    // Note that only receives will enter this function anyways, so no need to branch
    int* envelope_int[2] = { internal_timer_prop_int[2*msg_id], internal_timer_prop_int[2*msg_id+1] };
    double* envelope_double[2] = { remote_path_data, internal_timer_prop_double[msg_id] };
    char* envelope_char = internal_timer_prop_char[msg_id];
    for (auto k=0; k<symbol_path_select_size; k++){
      // Up until this very point, we had no idea whether we, or our partner rank, determined the path for a specific metric.
      if (envelope_double[0][symbol_path_select_index[k]] > critical_path_costs[symbol_path_select_index[k]]){
        int ftimer_size = *envelope_int[0];
        int symbol_offset = 0;
        for (int i=0; i<ftimer_size; i++){
          auto reconstructed_symbol = std::string(envelope_char+symbol_offset,envelope_char+symbol_offset+envelope_int[1][i]);
          if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
            symbol_timers[reconstructed_symbol] = symbol_tracker(reconstructed_symbol);
            symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
          }
          std::memcpy(symbol_timers[reconstructed_symbol].cp_numcalls[k],
                      &envelope_double[1][(i*symbol_path_select_size+k)*(cp_symbol_class_count*num_per_process_measures+1)],
                      sizeof(double)*(cp_symbol_class_count*num_per_process_measures+1));
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
  for (auto& it : internal_timer_prop_double){ free(it); }
  for (auto& it : internal_timer_prop_double_int){ free(it); }
  for (auto& it : internal_timer_prop_char){ free(it); }
  internal_timer_prop_int.clear(); internal_timer_prop_double.clear(); internal_timer_prop_double_int.clear();
  internal_timer_prop_char.clear(); internal_timer_prop_req.clear();
}


void decomposition::initiate(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                            bool is_sender, int partner1, int partner2){
  // Save and accumulate the computation time between last communication routine as both execution-time and computation time
  //   into both the execution-time critical path data structures and the per-process data structures.
  tracker.comp_time = curtime - computation_timer;

  assert(comm != 0);
  int rank; MPI_Comm_rank(comm, &rank);
  // We consider usage of Sendrecv variants to forfeit usage of eager internal communication.
  // Note that the reason we can't force user Bsends to be 'true_eager_p2p' is because the corresponding Receives would be expecting internal communications
  bool true_eager_p2p = ((eager_p2p == 1) && (tracker.tag!=13) && (tracker.tag!=14));
  if (true_eager_p2p){
    MPI_Buffer_attach(&eager_pad[0],eager_pad.size());
  }

  tracker.barrier_time=0.;// might get updated below
  if ((partner1==-1) || (track_p2p_idle==1)){// if blocking collective, or if p2p and idle time is requested to be tracked
    assert(partner1 != MPI_ANY_SOURCE);
    if ((tracker.tag == 13) || (tracker.tag == 14)){ assert(partner2 != MPI_ANY_SOURCE); }

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

    volatile double init_time = MPI_Wtime();
    if (partner1 == -1){ PMPI_Barrier(comm); }
    else {
      MPI_Request barrier_reqs[3]; int barrier_count=0;
      char sbuf='H'; char rbuf='H';
      if ((is_sender) && (rank != partner1)){
        if (true_eager_p2p) { PMPI_Bsend(&sbuf, 1, MPI_CHAR, partner1, internal_tag3, comm); }
        else                { PMPI_Issend(&sbuf, 1, MPI_CHAR, partner1, internal_tag3, comm, &barrier_reqs[barrier_count]); barrier_count++; }
      }
      if ((!is_sender) && (rank != partner1)){
        PMPI_Irecv(&rbuf, 1, MPI_CHAR, partner1, internal_tag3, comm, &barrier_reqs[barrier_count]); barrier_count++;
      }
      if ((partner2 != -1) && (rank != partner2)){
        PMPI_Irecv(&rbuf, 1, MPI_CHAR, partner2, internal_tag3, comm, &barrier_reqs[barrier_count]); barrier_count++;
      }
      PMPI_Waitall(barrier_count,&barrier_reqs[0],MPI_STATUSES_IGNORE);
    }
    tracker.barrier_time = MPI_Wtime() - init_time;

    // If eager protocol is enabled, its assumed that any message latency the sender incurs is negligable, and thus the receiver incurs its true idle time above
    // Again, the gray-area is with Sendrecv variants, and we assume they are treated without eager protocol
    if (!true_eager_p2p){
      // We need to subtract out the idle time of the path-root along the execution-time cp so that it appears as this path has no idle time.
      // Ideally we would do this for the last process to enter this barrier (which would always determine the execution-time cp anyway, but would apply for a path defined by any metric in its distribution).
      // This is more invasive for p2p, because it requires more handling on the nonblocking side in case of a nonblocking+blocking communication pattern.
      // TODO: This might function within a new 2-stage mechanism, which is where I want to go next. I think its more efficient as well. See notes.
      double min_idle_time=tracker.barrier_time;
      double recv_idle_time1=std::numeric_limits<double>::max();
      double recv_idle_time2=std::numeric_limits<double>::max();
      if (partner1 == -1){ PMPI_Allreduce(MPI_IN_PLACE, &min_idle_time, 1, MPI_DOUBLE, MPI_MIN, comm); }
      else {
        MPI_Request barrier_reqs[3]; int barrier_count=0;
        if ((is_sender) && (rank != partner1)){
          PMPI_Issend(&min_idle_time, 1, MPI_DOUBLE, partner1, internal_tag4, comm, &barrier_reqs[barrier_count]); barrier_count++;
        }
        if ((!is_sender) && (rank != partner1)){
          PMPI_Irecv(&recv_idle_time1, 1, MPI_DOUBLE, partner1, internal_tag4, comm, &barrier_reqs[barrier_count]); barrier_count++;
        }
        if ((partner2 != -1) && (rank != partner2)){
          PMPI_Irecv(&recv_idle_time2, 1, MPI_DOUBLE, partner2, internal_tag4, comm, &barrier_reqs[barrier_count]); barrier_count++;
        }
        PMPI_Waitall(barrier_count,&barrier_reqs[0],MPI_STATUSES_IGNORE);
        min_idle_time = std::min(min_idle_time,std::min(recv_idle_time1,recv_idle_time2));
      }
      tracker.barrier_time -= min_idle_time;
    }
    for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-1-i-comm_path_select_size] += tracker.barrier_time; }
  }

  critical_path_costs[num_critical_path_measures-2] += tracker.comp_time;	// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += tracker.comp_time;	// update critical path runtime
  volume_costs[num_volume_measures-2]        += tracker.comp_time;		// update local computation time
  volume_costs[num_volume_measures-1]        += tracker.comp_time;		// update local runtime
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-1-i] += tracker.comp_time; }// update each metric's critical path's computation time
  if (symbol_path_select_size>0 && symbol_stack.size()>0){
    // Get the current symbol's execution-time since last communication routine or its inception.
    // Accumulate as both execution-time and computation time into both the execution-time critical path data structures and the per-process data structures.
    auto last_symbol_time = curtime - symbol_timers[symbol_stack.top()].start_timer.top();
    for (auto i=0; i<symbol_path_select_size; i++){
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-1] += last_symbol_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-2] += last_symbol_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-1] += last_symbol_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-2] += last_symbol_time;
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-2] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-2] += last_symbol_time;
  }

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

  if ((partner1==-1) || (track_p2p_idle==1)){// if blocking collective, or if p2p and idle time is requested to be tracked
    assert(partner1 != MPI_ANY_SOURCE);
    if ((tracker.tag == 13) || (tracker.tag == 14)){ assert(partner2 != MPI_ANY_SOURCE); }

    // Use the user communication routine to measre synchronization time.
    // Note the following consequences of using a tiny 1-byte message (note that 0-byte is trivially handled by most MPI implementations) on measuring synchronization time:
    // 	1) The collective communication algorithm is likely different for small messages than large messages.
    // 	2) The eager sending protocol will be utilized, which would incur a potentially significant difference in synchronization time than if rendezvous protocol was invoked.
    // 		On second thought. I will force usage of Ssend. TODO: Check whether this breaks any correctness semantics.

    // Special arrays for use in the collective -v routines as well as Reduce_scatter.
    std::vector<int> counts(np,1); std::vector<int> disp(np,0);// TODO: could be moved as global variables
    if (tracker.tag>=9 && tracker.tag<=12) for (int i=1; i<np; i++) disp[i]=disp[i-1]+1;
    // start synchronization timer for communication routine
    tracker.start_time = MPI_Wtime();
    switch (tracker.tag){
      case 0:
        PMPI_Barrier(comm);
        break;
      case 1:
        PMPI_Bcast(&synch_pad_send[0], 1, MPI_CHAR, 0, comm);// arbitrary root 0
        break;
      case 2:
        PMPI_Reduce(&synch_pad_send[0], &synch_pad_recv[0], 1, MPI_CHAR, MPI_MAX, 0, comm);// arbitrary root 0
        break;
      case 3:
        PMPI_Allreduce(MPI_IN_PLACE, &synch_pad_send[0], 1, MPI_CHAR, MPI_MAX, comm);
        break;
      case 4:
        PMPI_Gather(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, 0, comm);// arbitrary root 0
        break;
      case 5:
        PMPI_Allgather(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, comm);
        break;
      case 6:
        PMPI_Scatter(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, 0, comm);// arbitrary root 0
        break;
      case 7:
        PMPI_Reduce_scatter(&synch_pad_send[0], &synch_pad_recv[0], &counts[0], MPI_CHAR, MPI_MAX, comm);
        break;
      case 8:
        PMPI_Alltoall(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, comm);
        break;
      case 9:
        PMPI_Gatherv(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], &counts[0], &disp[0], MPI_CHAR, 0, comm);// arbitrary root 0
        break;
      case 10:
        PMPI_Allgatherv(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], &counts[0], &disp[0], MPI_CHAR, comm);
        break;
      case 11:
        PMPI_Scatterv(&synch_pad_send[0], &counts[0], &disp[0], MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, 0, comm);// arbitrary root 0
        break;
      case 12:
        PMPI_Alltoallv(&synch_pad_send[0], &counts[0], &disp[0], MPI_CHAR, &synch_pad_recv[0], &counts[0], &disp[0], MPI_CHAR, comm);
        break;
      case 13:
        PMPI_Sendrecv(&synch_pad_send[0], 1, MPI_CHAR, partner1, internal_tag, &synch_pad_recv[0], 1, MPI_CHAR, partner2, internal_tag, comm, MPI_STATUS_IGNORE);
        break;
      case 14:
        PMPI_Sendrecv_replace(&synch_pad_send[0], 1, MPI_CHAR, partner1, internal_tag, partner2, internal_tag, comm, MPI_STATUS_IGNORE);
        break;
      case 15:
        PMPI_Ssend(&synch_pad_send[0], 1, MPI_CHAR, partner1, internal_tag, comm);
        break;
      case 16:
        if (true_eager_p2p) { PMPI_Bsend(&synch_pad_send[0], 1, MPI_CHAR, partner1, internal_tag, comm); }
        else            { PMPI_Ssend(&synch_pad_send[0], 1, MPI_CHAR, partner1, internal_tag, comm); }// forced usage of synchronous send to avoid eager sends for large messages.
        break;
      case 17:
        PMPI_Recv(&synch_pad_recv[0], 1, MPI_CHAR, partner1, internal_tag, comm, MPI_STATUS_IGNORE);
        break;
      case 32:
        PMPI_Bsend(&synch_pad_send[0], 1, MPI_CHAR, partner1, internal_tag, comm);
        break;
    }
    tracker.synch_time = MPI_Wtime()-tracker.start_time;
  }

  // start communication timer for communication routine
  tracker.start_time = MPI_Wtime();
}

// Used only for p2p communication. All blocking collectives use sychronous protocol
void decomposition::complete(blocking& tracker, int recv_source){
  // We handle wildcard sources (for MPI_Recv variants) only after the user communication.
  if (recv_source != -1){
    if ((tracker.tag == 13) || (tracker.tag == 14)){
      tracker.partner2=recv_source;
    }
    else{
      assert(tracker.tag==17);
      tracker.partner1=recv_source;
    }
  }
  volatile double comm_time = MPI_Wtime() - tracker.start_time;	// complete communication time
  std::pair<double,double> cost_bsp    = tracker.cost_func_bsp(tracker.nbytes, tracker.comm_size);
  std::pair<double,double> cost_alphabeta = tracker.cost_func_alphabeta(tracker.nbytes, tracker.comm_size);
  std::vector<std::pair<double,double>> costs = {cost_bsp,cost_alphabeta};

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
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-5] += tracker.barrier_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-4] += comm_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-3] += tracker.synch_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-1] += comm_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-5] += tracker.barrier_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-4] += comm_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-3] += tracker.synch_time;
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
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-5] += tracker.barrier_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-4] += comm_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-3] += tracker.synch_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += (comm_time+tracker.barrier_time);
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-5] += tracker.barrier_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-4] += comm_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-3] += tracker.synch_time;
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
  critical_path_costs[num_critical_path_measures-4] += comm_time;		// update critical path communication time (for what this process has seen thus far)
  critical_path_costs[num_critical_path_measures-3] += tracker.synch_time;	// update critical path synchronization time
  critical_path_costs[num_critical_path_measures-1] += comm_time;		// update critical path runtime

  volume_costs[num_volume_measures-5] += tracker.barrier_time;			// update local barrier/idle time
  volume_costs[num_volume_measures-4] += comm_time;				// update local communication time (not volume until after the completion of the program)
  volume_costs[num_volume_measures-3] += tracker.synch_time;			// update local synchronization time
  volume_costs[num_volume_measures-1] += (tracker.barrier_time+comm_time);	// update local runtime with idle time and comm time

  // Note that this block of code below is left in solely for blocking communication to avoid over-counting the idle time
  //   (which does not get subtracted by the min idle time any one process incurs due to efficiency complications with matching nonblocking+blocking p2p communications).
  //   Its handled correctly for blocking collectives.
  // If per-process execution-time gets larger than execution-time along the execution-time critical path, subtract out the difference from idle time.
  volume_costs[num_volume_measures-5] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
  if (symbol_path_select_size>0 && symbol_stack.size()>0){
    // Special handling of excessively large idle time caused by suspected tool interference
    // Specifically, this interference is caused by not subtracting out the barrier time of the last process to enter the barrier (which ideally is 0).
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1]     -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-5] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-5]     -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
  }

  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
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
  // Prepare to leave interception and re-enter user code by restarting computation timers.
  tracker.start_time = MPI_Wtime();
  computation_timer = tracker.start_time;
  if (symbol_path_select_size>0 && symbol_stack.size()>0){ symbol_timers[symbol_stack.top()].start_timer.top() = tracker.start_time; }
}

// Called by both nonblocking p2p and nonblocking collectives
void decomposition::initiate(nonblocking& tracker, volatile double curtime, volatile double itime, int64_t nelem,
                            MPI_Datatype t, MPI_Comm comm, MPI_Request* request, bool is_sender, int partner){

  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  tracker.comp_time = curtime - computation_timer + itime;
  critical_path_costs[num_critical_path_measures-2] += tracker.comp_time;		// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += tracker.comp_time;		// update critical path runtime
  volume_costs[num_volume_measures-2]        += tracker.comp_time;		// update local computation time
  volume_costs[num_volume_measures-1]        += tracker.comp_time;		// update local runtime
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-1-i] += tracker.comp_time; }
  if (symbol_path_select_size>0 && symbol_stack.size()>0){
    assert(symbol_stack.size()>0);
    assert(symbol_timers[symbol_stack.top()].start_timer.size()>0);
    double save_time = curtime - symbol_timers[symbol_stack.top()].start_timer.top()+itime;
    for (auto i=0; i<symbol_path_select_size; i++){
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-1] += save_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-2] += save_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-1] += save_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-2] += save_time;
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += save_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-2] += save_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-1] += save_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-2] += save_time;
  }

  int el_size,p;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(comm, &p);

  internal_comm_info[*request] = is_sender;
  internal_comm_comm[*request] = std::make_pair(comm,partner);// Note 'partner' might be MPI_ANY_SOURCE
  internal_comm_data[*request] = std::make_pair((double)nbytes,(double)p);
  internal_comm_track[*request] = &tracker;

  if (eager_p2p==1){
    tracker.comm = comm;
    tracker.is_sender = is_sender;
    tracker.partner1 = partner;
    tracker.partner2 = -1;
    tracker.comm_size = p;
    tracker.nbytes = nbytes;
    propagate(tracker);
  }

  tracker.start_time = MPI_Wtime();
  computation_timer = tracker.start_time;
  if (symbol_path_select_size>0 && symbol_stack.size()>0){ symbol_timers[symbol_stack.top()].start_timer.top() = tracker.start_time; }
}

void decomposition::complete(nonblocking& tracker, MPI_Request* request, double comp_time, double comm_time){
  auto comm_info_it = internal_comm_info.find(*request);
  auto comm_comm_it = internal_comm_comm.find(*request);
  auto comm_data_it = internal_comm_data.find(*request);
  auto comm_track_it = internal_comm_track.find(*request);
  assert(comm_info_it != internal_comm_info.end());
  assert(comm_comm_it != internal_comm_comm.end());
  assert(comm_data_it != internal_comm_data.end());
  assert(comm_track_it != internal_comm_track.end());

  tracker.is_sender = comm_info_it->second;
  tracker.comm = comm_comm_it->second.first;
  tracker.partner1 = comm_comm_it->second.second;
  tracker.partner2 = -1;
  tracker.nbytes = comm_data_it->second.first;
  tracker.comm_size = comm_data_it->second.second;
  tracker.synch_time=0;
  
  // Both sender and receiver will now update its critical path with the data from the communication
  std::pair<double,double> cost_bsp  = tracker.cost_func_bsp(tracker.nbytes,tracker.comm_size);
  std::pair<double,double> cost_alphabeta = tracker.cost_func_alphabeta(tracker.nbytes,tracker.comm_size);
  if ((tracker.tag<20) && (wait_id)) cost_bsp.first=1.;	// this is usually zero, but we force it to be 1 in special circumstances (for nonblocking p2p with wait_id one)
  std::vector<std::pair<double,double>> costs = {cost_bsp,cost_alphabeta};

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
  critical_path_costs[num_critical_path_measures-4] += comm_time;			// update critical path communication time (for what this process has seen thus far)
  critical_path_costs[num_critical_path_measures-3] += 0.;				// update critical path synchronization time
  critical_path_costs[num_critical_path_measures-2] += comp_time;			// update critical path runtime
  critical_path_costs[num_critical_path_measures-1] += comp_time+comm_time;		// update critical path runtime
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-1-i] += comp_time; }

  volume_costs[num_volume_measures-4] += comm_time;				// update local communication time (not volume until after the completion of the program)
  volume_costs[num_volume_measures-3] += 0.;					// update local synchronization time
  volume_costs[num_volume_measures-2] += comp_time;				// update local runtime
  volume_costs[num_volume_measures-1] += comp_time+comm_time;			// update local runtime
  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
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
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-4] += comm_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-3] += 0.;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-2] += comp_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-1] += (comp_time+comm_time);
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-4] += comm_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-3] += 0.;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-2] += comp_time;
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
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-4] += comm_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-3] += 0.;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-3] += comp_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += (comp_time+comm_time);
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-4] += comm_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-3] += 0.;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-2] += comp_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += (comp_time+comm_time);
  }

  if (eager_p2p==0) { propagate(tracker); }

  internal_comm_info.erase(*request);
  internal_comm_comm.erase(*request);
  internal_comm_data.erase(*request);
  internal_comm_track.erase(*request);

  tracker.start_time = MPI_Wtime();
}

void decomposition::complete(double curtime, MPI_Request* request, MPI_Status* status){
  double comp_time = curtime - computation_timer;
  auto comm_track_it = internal_comm_track.find(*request);
  assert(comm_track_it != internal_comm_track.end());
  auto comm_info_it = internal_comm_info.find(*request);
  auto comm_comm_it = internal_comm_comm.find(*request);
  MPI_Request save_request = comm_info_it->first;
  if ((comm_comm_it->second.second!=-1) && (track_p2p_idle==1)){// if p2p and idle time is requested to be tracked (first case prevents nonblocking collectives
    assert(comm_comm_it->second.first != 0);
    int comm_rank; MPI_Comm_rank(comm_comm_it->second.first,&comm_rank); 
    double max_barrier_time = 0;// counter-intuitively, a blocking partner should determine the idle time
    if (comm_info_it->second && comm_rank != comm_comm_it->second.second){
      PMPI_Send(&barrier_pad_send[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag3, comm_comm_it->second.first);
      PMPI_Send(&max_barrier_time, 1, MPI_DOUBLE, comm_comm_it->second.second, internal_tag4, comm_comm_it->second.first);
      PMPI_Send(&synch_pad_send[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag, comm_comm_it->second.first);
    }
    else if (!comm_info_it->second && comm_rank != comm_comm_it->second.second){
      PMPI_Recv(&barrier_pad_recv[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag3, comm_comm_it->second.first, MPI_STATUS_IGNORE);
      PMPI_Recv(&max_barrier_time, 1, MPI_DOUBLE, comm_comm_it->second.second, internal_tag4, comm_comm_it->second.first, MPI_STATUS_IGNORE);
      PMPI_Recv(&synch_pad_recv[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag, comm_comm_it->second.first, MPI_STATUS_IGNORE);
    }
  }
  volatile double last_start_time = MPI_Wtime();
  PMPI_Wait(request, status);
  double save_comm_time = MPI_Wtime() - last_start_time;
  if (eager_p2p==1) { complete_path_update(); }
  if (comm_comm_it->second.second == MPI_ANY_SOURCE) { comm_track_it->second->partner1 = status->MPI_SOURCE; }
  complete(*comm_track_it->second, &save_request, comp_time, save_comm_time);
  if (eager_p2p==0) { complete_path_update(); }
  computation_timer = MPI_Wtime();
  if (symbol_path_select_size>0){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
}

void decomposition::complete(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){

  double waitany_comp_time = curtime - computation_timer;
  // We must force 'track_p2p_idle' to be zero because we don't know which request the MPI implementation will choose before
  //   it chooses it. Note that this is a not a problem for MPI_Waitall because all requests are chosen.
  //   Thus, we cannot participate in any idle/synch time exchanges. This is not a big deal at all if MPI_Waitany
  //   is used for both sender and receiver side, as idle/synch time are zero anyways. This decision is limiting in the sense that
  //   if MPI_Waitany is used to close nonblocking requests on both the sender and receiver side, then we forfeit tracking of idle/synch
  //   time caused by blocking send/recv.
  // TODO: For that reason alone, it may make sense to remove the assert and inform the user about the possibility of a hang if disobeying our rules.
  assert(track_p2p_idle==0);
  // We must save the requests before the completition of a request by the MPI implementation because its tag is set to MPI_REQUEST_NULL and lost forever
  std::vector<MPI_Request> pt(count); for (int i=0;i<count;i++){pt[i]=(array_of_requests)[i];}
  volatile double last_start_time = MPI_Wtime();
  PMPI_Waitany(count,array_of_requests,indx,status);
  double waitany_comm_time = MPI_Wtime() - last_start_time;
  if (eager_p2p==1) { complete_path_update(); }
  MPI_Request request = pt[*indx];
  auto comm_track_it = internal_comm_track.find(request);
  auto comm_comm_it = internal_comm_comm.find(request);
  assert(comm_track_it != internal_comm_track.end());
  if (comm_comm_it->second.second == MPI_ANY_SOURCE) { comm_track_it->second->partner1 = status->MPI_SOURCE; }
  complete(*comm_track_it->second, &request, waitany_comp_time, waitany_comm_time);
  if (eager_p2p==0) { complete_path_update(); }
  computation_timer = MPI_Wtime();
  if (symbol_path_select_size>0){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
}

void decomposition::complete(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[],
                        MPI_Status array_of_statuses[]){

  double waitsome_comp_time = curtime - computation_timer;
  wait_id=true;
  // Read comment in function above. Same ideas apply for Waitsome.
  assert(track_p2p_idle==0);
  // We must save the requests before the completition of a request by the MPI implementation because its tag is set to MPI_REQUEST_NULL and lost forever
  std::vector<MPI_Request> pt(incount); for (int i=0;i<incount;i++){pt[i]=(array_of_requests)[i];}
  volatile double last_start_time = MPI_Wtime();
  PMPI_Waitsome(incount,array_of_requests,outcount,array_of_indices,array_of_statuses);
  double waitsome_comm_time = MPI_Wtime() - last_start_time;
  if (eager_p2p==1) { complete_path_update(); }
  for (int i=0; i<*outcount; i++){
    MPI_Request request = pt[(array_of_indices)[i]];
    auto comm_track_it = internal_comm_track.find(request);
    auto comm_comm_it = internal_comm_comm.find(request);
    assert(comm_track_it != internal_comm_track.end());
    if (comm_comm_it->second.second == MPI_ANY_SOURCE) { comm_track_it->second->partner1 = (array_of_statuses)[i].MPI_SOURCE; }
    complete(*comm_track_it->second, &request, waitsome_comp_time, waitsome_comm_time);
    waitsome_comp_time=0;
    waitsome_comm_time=0;
    if (i==0){wait_id=false;}
  }
  if (eager_p2p==0) { complete_path_update(); }
  computation_timer = MPI_Wtime();
  if (symbol_path_select_size>0){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
}

void decomposition::complete(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
  double waitall_comp_time = curtime - computation_timer;
  wait_id=true;
  if (track_p2p_idle==1){// nonblocking collectives won't pass the if statements below anyway.
    std::vector<MPI_Request> internal_requests(3*count,MPI_REQUEST_NULL);
    if (count > barrier_pad_send.size()){
      barrier_pad_send.resize(count);
      barrier_pad_recv.resize(count);
      synch_pad_send.resize(count);
      synch_pad_recv.resize(count);
    }
    // Issue all barrier/synch communications at once because request order is not guaranteed to be sequenced together on all processes.
    // Necessary to avoid corruption of idle time calculation that would occur if sending out in some sequence after each request is completed.
    // Presumably the sending communications will utilize the eager protocol, but as the Waitall is issued immediately following the loop, its irrelevant.
    // TODO: Staging nonblocking receives with nonblocking sends might be a problem because the Sends will likely use the eager protocol (since 1-byte messages),
    //         while the Recvs will block until the message buffer is available. This will force the matching blocking receives of the nonblocking senders to possibely start counting communication time while
    //         the sender is stuck in the Waitall loop with the other sends and receives, thus corrupting the measurement of communication time. I propose that the Sends issue their 3-messages together followed
    //         by a Waitall. At this point, the choice is whether to issue the user-communication sends one-by-one via Waitany loop, issue the nonblocking recv 3-messages together followed by waitall,
    //         or actually issue a Waitall for the user communication rather than a loop over Waits. I am now leaning towards supporting the user-communication Waitall instead, as that might be a source of overhead
    //         with CTF. Each process can utilize its timer and record the same communication time, and then issue the exchange of path information via nonblocking communications.
    double max_barrier_time = 0;// counter-intuitively, a blocking partner should determine the idle time
    for (int i=0; i<count; i++){
      auto comm_info_it = internal_comm_info.find(*(array_of_requests+i));
      assert(comm_info_it != internal_comm_info.end());
      auto comm_comm_it = internal_comm_comm.find(*(array_of_requests+i));
      assert(comm_comm_it != internal_comm_comm.end());
      assert(comm_comm_it->second.second != MPI_ANY_SOURCE);
      if (comm_info_it->second && comm_comm_it->second.second != -1){
        PMPI_Isend(&barrier_pad_send[i], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag3,
          comm_comm_it->second.first, &internal_requests[3*i]);
        if (eager_p2p==0) { PMPI_Isend(&max_barrier_time, 1, MPI_DOUBLE, comm_comm_it->second.second, internal_tag4,
          comm_comm_it->second.first, &internal_requests[3*i+1]); }
        PMPI_Isend(&synch_pad_send[i], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag,
          comm_comm_it->second.first, &internal_requests[3*i+2]);
      }
      else if (!comm_info_it->second && comm_comm_it->second.second != -1){
        PMPI_Irecv(&barrier_pad_recv[i], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag3,
          comm_comm_it->second.first, &internal_requests[3*i]);
        if (eager_p2p==0) { PMPI_Irecv(&max_barrier_time, 1, MPI_DOUBLE, comm_comm_it->second.second, internal_tag4,
          comm_comm_it->second.first, &internal_requests[3*i+1]); }
        PMPI_Irecv(&synch_pad_recv[i], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag,
          comm_comm_it->second.first, &internal_requests[3*i+2]);
      }
    }
    PMPI_Waitall(internal_requests.size(), &internal_requests[0], MPI_STATUSES_IGNORE);
  }
  // We must save the requests before the completition of a request by the MPI implementation because its tag is set to MPI_REQUEST_NULL and lost forever
  std::vector<MPI_Request> pt(count); for (int i=0;i<count;i++){pt[i]=(array_of_requests)[i];}
  volatile double last_start_time = MPI_Wtime();
  PMPI_Waitall(count,array_of_requests,array_of_statuses);
  double waitall_comm_time = MPI_Wtime() - last_start_time;
  if (eager_p2p==1) { complete_path_update(); }
  for (int i=0; i<count; i++){
    MPI_Request request = pt[i];
    auto comm_track_it = internal_comm_track.find(request);
    auto comm_comm_it = internal_comm_comm.find(request);
    assert(comm_track_it != internal_comm_track.end());
    if (comm_comm_it->second.second == MPI_ANY_SOURCE) { comm_track_it->second->partner1 = (array_of_statuses)[i].MPI_SOURCE; }
    complete(*comm_track_it->second, &request, waitall_comp_time, waitall_comm_time);
    // Although we have to exchange the path data for each request, we do not want to double-count the computation time nor the communicaion time
    waitall_comp_time=0;
    waitall_comm_time=0;
    if (i==0){wait_id=false;}
  }
  wait_id=true;
  if (eager_p2p==0) { complete_path_update(); }
  computation_timer = MPI_Wtime();
  if (symbol_path_select_size>0){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
}

void decomposition::propagate_symbols(nonblocking& tracker, int rank){
  if (eager_p2p==0){
    MPI_Request internal_request[8];
    int* send_envelope1 = nullptr; int* send_envelope2 = nullptr; double* send_envelope3 = nullptr; char* send_envelope5 = nullptr;
    int* recv_envelope1 = nullptr; int* recv_envelope2 = nullptr; double* recv_envelope3 = nullptr; char* recv_envelope5 = nullptr;
    int ftimer_size = symbol_timers.size();
    int data_len_size = symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*ftimer_size;
    int num_chars = 0;
    for (int i=0; i<ftimer_size; i++) { num_chars += symbol_order[i].size(); }
    send_envelope1 = (int*)malloc(sizeof(int)); *send_envelope1 = ftimer_size;
    send_envelope2 = (int*)malloc(sizeof(int)*(ftimer_size));
    send_envelope3 = (double*)malloc(sizeof(double)*data_len_size);
    send_envelope5 = (char*)malloc(sizeof(char)*num_chars);
    int symbol_offset = 0;
    for (auto i=0; i<ftimer_size; i++){
      send_envelope2[i] = symbol_order[i].size();
      for (auto j=0; j<symbol_order[i].size(); j++){
        send_envelope5[symbol_offset+j] = symbol_order[i][j];
      }
      symbol_offset += symbol_order[i].size();
    }
    std::memcpy(send_envelope3,&symbol_timer_pad_local_cp[0],sizeof(double)*data_len_size);
    PMPI_Isend(&send_envelope1[0],1,MPI_INT,tracker.partner1,internal_tag1,tracker.comm,&internal_request[0]);
    PMPI_Isend(&send_envelope2[0],ftimer_size,MPI_INT,tracker.partner1,internal_tag2,tracker.comm,&internal_request[1]);
    PMPI_Isend(&send_envelope3[0],data_len_size,MPI_DOUBLE,tracker.partner1,internal_tag3,tracker.comm,&internal_request[2]);
    PMPI_Isend(&send_envelope5[0],symbol_offset,MPI_CHAR,tracker.partner1,internal_tag5,tracker.comm,&internal_request[3]);

    recv_envelope1 = (int*)malloc(sizeof(int));
    recv_envelope2 = (int*)malloc(sizeof(int)*(max_num_symbols));
    recv_envelope3 = (double*)malloc(sizeof(double)*symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*max_num_symbols);
    recv_envelope5 = (char*)malloc(sizeof(char)*max_timer_name_length*max_num_symbols);
    PMPI_Irecv(recv_envelope1,1,MPI_INT,tracker.partner1,internal_tag1,tracker.comm,&internal_request[4]);
    PMPI_Irecv(recv_envelope2,max_num_symbols,MPI_INT,tracker.partner1,internal_tag2,tracker.comm,&internal_request[5]);
    PMPI_Irecv(recv_envelope3,symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*max_num_symbols,MPI_DOUBLE,tracker.partner1,internal_tag3,tracker.comm,&internal_request[6]);
    PMPI_Irecv(recv_envelope5,max_timer_name_length*max_num_symbols,MPI_CHAR,tracker.partner1,internal_tag5,tracker.comm,&internal_request[7]);

    for (int i=0; i<8; i++) { internal_timer_prop_req.push_back(internal_request[i]); }
    internal_timer_prop_int.push_back(send_envelope1); internal_timer_prop_int.push_back(send_envelope2);
    internal_timer_prop_int.push_back(recv_envelope1); internal_timer_prop_int.push_back(recv_envelope2);
    internal_timer_prop_double.push_back(send_envelope3);
    internal_timer_prop_double.push_back(recv_envelope3);
    internal_timer_prop_char.push_back(send_envelope5); internal_timer_prop_char.push_back(recv_envelope5);
  } else{
    if (tracker.is_sender){
      MPI_Request internal_request[4];
      int* send_envelope1 = nullptr; int* send_envelope2 = nullptr; double* send_envelope3 = nullptr; char* send_envelope5 = nullptr;
      int ftimer_size = symbol_timers.size();
      int data_len_size = symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*ftimer_size;
      int num_chars = 0;
      for (int i=0; i<ftimer_size; i++) { num_chars += symbol_order[i].size(); }
      send_envelope1 = (int*)malloc(sizeof(int)); *send_envelope1 = ftimer_size;
      send_envelope2 = (int*)malloc(sizeof(int)*(ftimer_size));
      send_envelope3 = (double*)malloc(sizeof(double)*data_len_size);
      send_envelope5 = (char*)malloc(sizeof(char)*num_chars);
      int symbol_offset = 0;
      for (auto i=0; i<ftimer_size; i++){
        send_envelope2[i] = symbol_order[i].size();
        for (auto j=0; j<symbol_order[i].size(); j++){
          send_envelope5[symbol_offset+j] = symbol_order[i][j];
        }
        symbol_offset += symbol_order[i].size();
      }
      std::memcpy(send_envelope3,&symbol_timer_pad_local_cp[0],sizeof(double)*data_len_size);
      PMPI_Isend(&send_envelope1[0],1,MPI_INT,tracker.partner1,internal_tag1,tracker.comm,&internal_request[0]);
      PMPI_Isend(&send_envelope2[0],ftimer_size,MPI_INT,tracker.partner1,internal_tag2,tracker.comm,&internal_request[1]);
      PMPI_Isend(&send_envelope3[0],data_len_size,MPI_DOUBLE,tracker.partner1,internal_tag3,tracker.comm,&internal_request[2]);
      PMPI_Isend(&send_envelope5[0],symbol_offset,MPI_CHAR,tracker.partner1,internal_tag5,tracker.comm,&internal_request[3]);

      for (int i=0; i<4; i++) { internal_timer_prop_req.push_back(internal_request[i]); }
      internal_timer_prop_int.push_back(send_envelope1);
      internal_timer_prop_int.push_back(send_envelope2);
      internal_timer_prop_double.push_back(send_envelope3);
      internal_timer_prop_char.push_back(send_envelope5);
    } else{
      MPI_Request internal_request[4];
      int* recv_envelope1 = nullptr; int* recv_envelope2 = nullptr; double* recv_envelope3 = nullptr; char* recv_envelope5 = nullptr;
      recv_envelope1 = (int*)malloc(sizeof(int));
      recv_envelope2 = (int*)malloc(sizeof(int)*(max_num_symbols));
      recv_envelope3 = (double*)malloc(sizeof(double)*symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*max_num_symbols);
      recv_envelope5 = (char*)malloc(sizeof(char)*max_timer_name_length*max_num_symbols);
      PMPI_Irecv(recv_envelope1,1,MPI_INT,tracker.partner1,internal_tag1,tracker.comm,&internal_request[0]);
      PMPI_Irecv(recv_envelope2,max_num_symbols,MPI_INT,tracker.partner1,internal_tag2,tracker.comm,&internal_request[1]);
      PMPI_Irecv(recv_envelope3,symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1)*max_num_symbols,MPI_DOUBLE,tracker.partner1,internal_tag3,tracker.comm,&internal_request[2]);
      PMPI_Irecv(recv_envelope5,max_timer_name_length*max_num_symbols,MPI_CHAR,tracker.partner1,internal_tag5,tracker.comm,&internal_request[3]);

      for (int i=0; i<4; i++) { internal_timer_prop_req.push_back(internal_request[i]); }
      internal_timer_prop_int.push_back(recv_envelope1);
      internal_timer_prop_int.push_back(recv_envelope2);
      internal_timer_prop_double.push_back(recv_envelope3);
      internal_timer_prop_char.push_back(recv_envelope5);
    }
  }
}

/*
 Its important to note here that a blocking p2p call will already know whether its the cp root or not, regardless of whether its partner used a nonblocking p2p routine.
   But, because that potential nonblocking partner does not have this knowledge, and thus posted both sends and recvs, the blocking partner also has to do so as well, even if its partner (unknown to him) used a blocking p2p routine.
*/
void decomposition::propagate_symbols(blocking& tracker, int rank){
  bool true_eager_p2p = ((eager_p2p == 1) && (tracker.tag!=13) && (tracker.tag!=14));
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
                      sizeof(double)*(cp_symbol_class_count*num_per_process_measures+1));
          for (auto j=0; j<symbol_len_pad_cp[symbol_offset_cp]; j++){
            symbol_pad_cp[char_count_cp+j] = symbol_order[i][j];
          }
          pad_global_offset += (cp_symbol_class_count*num_per_process_measures+1);
          char_count_cp += symbol_len_pad_cp[symbol_offset_cp++];
        }
      }
      else{
        memset(&symbol_timer_pad_global_cp[pad_global_offset],0,sizeof(double)*ftimer_size_cp[k]*(cp_symbol_class_count*num_per_process_measures+1));
        pad_global_offset += ftimer_size_cp[k]*(cp_symbol_class_count*num_per_process_measures+1);
        for (auto i=0; i<ftimer_size_cp[k]; i++){
          memset(&symbol_pad_cp[char_count_cp],0,sizeof(char)*symbol_len_pad_cp[symbol_offset_cp]);
          char_count_cp += symbol_len_pad_cp[symbol_offset_cp++];
        }
      }
    }
    PMPI_Allreduce(MPI_IN_PLACE,&symbol_timer_pad_global_cp[0],pad_global_offset,MPI_DOUBLE,MPI_SUM,tracker.comm);
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
                      sizeof(double)*(cp_symbol_class_count*num_per_process_measures+1));
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
  else if (!true_eager_p2p){
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
    PMPI_Isend(&symbol_timer_pad_local_cp[0],data_len_cp,MPI_DOUBLE,tracker.partner1,internal_tag3,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
    PMPI_Isend(&symbol_pad_cp[0],char_count_cp,MPI_CHAR,tracker.partner1,internal_tag5,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
    PMPI_Irecv(&symbol_timer_pad_global_cp[0],data_len_ncp1,MPI_DOUBLE,tracker.partner1,internal_tag3,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
    PMPI_Irecv(&symbol_pad_ncp1[0],char_count_ncp1,MPI_CHAR,tracker.partner1,internal_tag5,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
    if (tracker.partner1 != tracker.partner2){
      PMPI_Isend(&symbol_timer_pad_local_cp[0],data_len_cp,MPI_DOUBLE,tracker.partner2,internal_tag3,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
      PMPI_Isend(&symbol_pad_cp[0],char_count_cp,MPI_CHAR,tracker.partner2,internal_tag5,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
      PMPI_Irecv(&symbol_timer_pad_global_cp2[0],data_len_ncp2,MPI_DOUBLE,tracker.partner2,internal_tag3,tracker.comm,&symbol_exchance_reqs[exchange_count]); exchange_count++;
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
                      sizeof(double)*(cp_symbol_class_count*num_per_process_measures+1));
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
                      sizeof(double)*(cp_symbol_class_count*num_per_process_measures+1));
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
      PMPI_Bsend(&symbol_timer_pad_local_cp[0],data_len_cp,MPI_DOUBLE,tracker.partner1,internal_tag3,tracker.comm);
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
      PMPI_Recv(&symbol_timer_pad_global_cp[0],data_len_cp,MPI_DOUBLE,tracker.partner1,internal_tag3,tracker.comm,MPI_STATUS_IGNORE);
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
                      sizeof(double)*(cp_symbol_class_count*num_per_process_measures+1));
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

void decomposition::propagate(blocking& tracker){
  assert(tracker.comm != 0);
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  if ((rank == tracker.partner1) && (rank == tracker.partner2)) { return; } 
  bool true_eager_p2p = ((eager_p2p == 1) && (tracker.tag!=13) && (tracker.tag!=14));
  if (symbol_path_select_size>0){
    //TODO: Idea for 2-stage reduction: move this out of the mode>=2 if statement, and then after this, scan the critical_path_costs and zero out what is not defining a critical path and then post a MPI_Allreduce (via multi-root hack)
    for (int i=0; i<num_critical_path_measures; i++){
      info_sender[i].first = critical_path_costs[i];
      info_sender[i].second = rank;
    }
    if (tracker.partner1 == -1){
      PMPI_Allreduce(&info_sender[0].first, &info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, MPI_MAXLOC, tracker.comm);
    }
    else{
      if (!true_eager_p2p){
        PMPI_Sendrecv(&info_sender[0].first, num_critical_path_measures, MPI_DOUBLE_INT, tracker.partner1, internal_tag,
                      &info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, tracker.partner2, internal_tag, tracker.comm, MPI_STATUS_IGNORE);
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
          PMPI_Sendrecv(&info_sender[0].first, num_critical_path_measures, MPI_DOUBLE_INT, tracker.partner2, internal_tag,
                        &info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, tracker.partner1, internal_tag, tracker.comm, MPI_STATUS_IGNORE);
          for (int i=0; i<num_critical_path_measures; i++){
            if (info_sender[i].first>info_receiver[i].first){info_receiver[i].second = rank;}
            else if (info_sender[i].first==info_receiver[i].first){ info_receiver[i].second = std::min(rank,tracker.partner1); }
            info_receiver[i].first = std::max(info_sender[i].first, info_receiver[i].first);
          }
        }
      }
      else{
        if (tracker.is_sender){
          PMPI_Bsend(&info_sender[0].first, num_critical_path_measures, MPI_DOUBLE_INT, tracker.partner1, internal_tag, tracker.comm);
        } else{
          PMPI_Recv(&info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, tracker.partner1, internal_tag, tracker.comm, MPI_STATUS_IGNORE);
        }
      }
    }
  }
  // Exchange the tracked routine critical path data
  if (tracker.partner1 == -1){
    MPI_Op op; MPI_Op_create((MPI_User_function*) propagate_critical_path_op,0,&op);
    PMPI_Allreduce(MPI_IN_PLACE, &critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, op, tracker.comm);
    MPI_Op_free(&op);
  }
  else{
    // Note that a blocking sendrecv allows exchanges even when the other party issued a request via nonblocking communication, as the process with the nonblocking request posts both sends and receives.
    if (true_eager_p2p){ PMPI_Bsend(&critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, tracker.partner1, internal_tag2, tracker.comm); }
    else { PMPI_Sendrecv(&critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, tracker.partner1, internal_tag2, &new_cs[0], critical_path_costs.size(),
                         MPI_DOUBLE, tracker.partner2, internal_tag2, tracker.comm, MPI_STATUS_IGNORE); }
    update_critical_path(&new_cs[0],&critical_path_costs[0],critical_path_costs_size);
    if (tracker.partner2 != tracker.partner1){
      // This if-statement will never be breached if 'true_eager_p2p'=true anyways.
      PMPI_Sendrecv(&critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, tracker.partner2, internal_tag2, &new_cs[0], critical_path_costs.size(), MPI_DOUBLE, tracker.partner1, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
      update_critical_path(&new_cs[0],&critical_path_costs[0],critical_path_costs_size);
    }
  }
  if (symbol_path_select_size>0) { propagate_symbols(tracker,rank); }
  if (true_eager_p2p){
    void* temp_buf; int temp_size;
    // Forces buffered messages to send. Ideally we should wait till the next invocation of 'decomposition::initiate(blocking&,...)' to call this,
    //   but to be safe and avoid stalls caused by MPI implementation not sending until this routine is called, we call it here.
    MPI_Buffer_detach(&temp_buf,&temp_size);
  }
}

void decomposition::propagate(nonblocking& tracker){
  assert(tracker.comm != 0);
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  if (rank == tracker.partner1) { return; } 
  if (symbol_path_select_size>0){
    for (int i=0; i<num_critical_path_measures; i++){
      info_sender[i].first = critical_path_costs[i];
      info_sender[i].second = rank;
    }
    if (tracker.partner1 == -1){ assert(0); }
    else if (eager_p2p==0){
      MPI_Request req1,req2;
      double_int* send_pathdata = (double_int*)malloc(num_critical_path_measures*sizeof(double_int));
      double_int* recv_pathdata = (double_int*)malloc(num_critical_path_measures*sizeof(double_int));
      memcpy(&send_pathdata[0].first, &info_sender[0].first, num_critical_path_measures*sizeof(double_int));
      PMPI_Isend(&send_pathdata[0].first, num_critical_path_measures, MPI_DOUBLE_INT, tracker.partner1, internal_tag, tracker.comm, &req1);
      PMPI_Irecv(&recv_pathdata[0].first, num_critical_path_measures, MPI_DOUBLE_INT, tracker.partner1, internal_tag, tracker.comm, &req2);
      internal_timer_prop_req.push_back(req1); internal_timer_prop_req.push_back(req2);
      internal_timer_prop_double_int.push_back(send_pathdata); internal_timer_prop_double_int.push_back(recv_pathdata);
    }
    else{
      if (tracker.is_sender){
        MPI_Request req1;
        double_int* send_pathdata = (double_int*)malloc(num_critical_path_measures*sizeof(double_int));
        memcpy(&send_pathdata[0].first, &info_sender[0].first, num_critical_path_measures*sizeof(double_int));
        PMPI_Isend(&send_pathdata[0].first, num_critical_path_measures, MPI_DOUBLE_INT, tracker.partner1, internal_tag, tracker.comm, &req1);
        internal_timer_prop_req.push_back(req1);
        internal_timer_prop_double_int.push_back(send_pathdata);
      } else{
        MPI_Request req1;
        double_int* recv_pathdata = (double_int*)malloc(num_critical_path_measures*sizeof(double_int));
        PMPI_Irecv(&recv_pathdata[0].first, num_critical_path_measures, MPI_DOUBLE_INT, tracker.partner1, internal_tag, tracker.comm, &req1);
        internal_timer_prop_req.push_back(req1);
        internal_timer_prop_double_int.push_back(recv_pathdata);
      }
    }
  }
  // Exchange the tracked routine critical path data
  if (tracker.partner1 == -1){
    MPI_Op op; MPI_Op_create((MPI_User_function*) propagate_critical_path_op,0,&op);
    MPI_Request req1;
    double* local_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
    std::memcpy(local_path_data, &critical_path_costs[0], critical_path_costs.size()*sizeof(double));
    PMPI_Iallreduce(MPI_IN_PLACE,local_path_data,critical_path_costs.size(),MPI_DOUBLE,op,tracker.comm,&req1);
    //MPI_Op_free(&op);
    internal_comm_prop.push_back(std::make_pair(local_path_data,true));
    internal_comm_prop_req.push_back(req1);
  }
  else if (eager_p2p==0){
    MPI_Request req1,req2;
    double* local_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
    std::memcpy(local_path_data, &critical_path_costs[0], critical_path_costs.size()*sizeof(double));
    double* remote_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
    PMPI_Isend(local_path_data, critical_path_costs.size(), MPI_DOUBLE, tracker.partner1, internal_tag2, tracker.comm, &req1);
    PMPI_Irecv(remote_path_data, critical_path_costs.size(), MPI_DOUBLE, tracker.partner1, internal_tag2, tracker.comm, &req2);
    internal_comm_prop.push_back(std::make_pair(local_path_data,true));
    internal_comm_prop_req.push_back(req1);
    internal_comm_prop.push_back(std::make_pair(remote_path_data,false));
    internal_comm_prop_req.push_back(req2);
  }
  else{
    MPI_Request req1;
    if (tracker.is_sender){
      double* local_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
      std::memcpy(local_path_data, &critical_path_costs[0], critical_path_costs.size()*sizeof(double));
      PMPI_Isend(local_path_data, critical_path_costs.size(), MPI_DOUBLE, tracker.partner1, internal_tag2, tracker.comm, &req1);
      internal_comm_prop.push_back(std::make_pair(local_path_data,true));
      internal_comm_prop_req.push_back(req1);
    }
    else{
      double* remote_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
      PMPI_Irecv(remote_path_data, critical_path_costs.size(), MPI_DOUBLE, tracker.partner1, internal_tag2, tracker.comm, &req1);
      internal_comm_prop.push_back(std::make_pair(remote_path_data,false));
      internal_comm_prop_req.push_back(req1);
    }
  }
  if (symbol_path_select_size>0) { propagate_symbols(tracker,rank); }
}

void decomposition::final_accumulate(double last_time){
  critical_path_costs[num_critical_path_measures-2]+=(last_time-computation_timer);	// update critical path computation time
  critical_path_costs[num_critical_path_measures-1]+=(last_time-computation_timer);	// update critical path runtime
  volume_costs[num_volume_measures-2]+=(last_time-computation_timer);			// update computation time volume
  volume_costs[num_volume_measures-1]+=(last_time-computation_timer);			// update runtime volume
  // update the computation time (i.e. time between last MPI synchronization point and this function invocation) along all paths decomposed by MPI communication routine
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-1-i] += (last_time-computation_timer); }
}

}
}
