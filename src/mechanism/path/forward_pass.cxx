#include "forward_pass.h"
#include "../../container/symbol_tracker.h"
#include "../../util.h"

namespace critter{
namespace internal{

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
  if (breakdown_size > 0){
    size_t breakdown_idx=0;
    for (int i=0; i<num_critical_path_measures; i++){
      if (breakdown[i]=='1'){ decisions[breakdown_idx++] = inout[i] > in[i]; }
    }
    for (int i=0; i<num_critical_path_measures; i++){
      inout[i] = std::max(inout[i],in[i]);
    }
    for (int i=num_critical_path_measures; i<critical_path_costs_size; i++){
      int idx = (i-num_critical_path_measures)%breakdown_size;
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
  int* envelope_int[2] = { internal_timer_prop_int[4*msg_id+2], internal_timer_prop_int[4*msg_id+3] };
  double* envelope_double[3] = { remote_path_data, internal_timer_prop_double[4*msg_id+2], internal_timer_prop_double[4*msg_id+3] };
  char* envelope_char = internal_timer_prop_char[2*msg_id+1];
  if (envelope_double[0][num_critical_path_measures-1] > critical_path_costs[num_critical_path_measures-1]){
    int ftimer_size = *envelope_int[0];
    int symbol_offset = 0;
    for (int i=0; i<ftimer_size; i++){
      auto reconstructed_symbol = std::string(envelope_char+symbol_offset,envelope_char+symbol_offset+envelope_int[1][i]);
      if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
        symbol_timers[reconstructed_symbol] = symbol_tracker(reconstructed_symbol);
        symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
      }
      *symbol_timers[reconstructed_symbol].cp_numcalls = envelope_double[1][(num_ftimer_measures*num_critical_path_measures+1)*i];
      for (int j=0; j<num_critical_path_measures; j++){
        *symbol_timers[reconstructed_symbol].cp_incl_measure[j] = envelope_double[1][(num_ftimer_measures*num_critical_path_measures+1)*i+2*j+1];
        *symbol_timers[reconstructed_symbol].cp_excl_measure[j] = envelope_double[1][(num_ftimer_measures*num_critical_path_measures+1)*i+2*(j+1)];
      }
      for (int j=0; j<symbol_timers[reconstructed_symbol].cp_exclusive_measure.size(); j++){ symbol_timers[reconstructed_symbol].cp_exclusive_measure[j]=0.; }
      for (int j=0; j<num_critical_path_measures; j++){
        symbol_timers[reconstructed_symbol].cp_exclusive_contributions[j] = envelope_double[2][i*num_critical_path_measures+j];
      }
      for (int k=0; k<num_critical_path_measures; k++){
        symbol_timers[reconstructed_symbol].cp_exclusive_measure[k] = envelope_double[2][ftimer_size*num_critical_path_measures+i*num_critical_path_measures+k];
      }
      symbol_timers[reconstructed_symbol].has_been_processed = true;
      symbol_offset += envelope_int[1][i];
    }
  }
}

static void complete_path_update(){
  PMPI_Waitall(internal_comm_prop_req.size(), &internal_comm_prop_req[0], MPI_STATUSES_IGNORE);
  if (mode>=2) { PMPI_Waitall(internal_timer_prop_req.size(), &internal_timer_prop_req[0], MPI_STATUSES_IGNORE); }
  size_t msg_id=0;
  for (auto& it : internal_comm_prop){
    if (!it.second){
      if (mode>=2) complete_timers(it.first,msg_id++);
      update_critical_path(it.first,&critical_path_costs[0],critical_path_costs_size);
    }
    free(it.first);
  }
  internal_comm_prop.clear(); internal_comm_prop_req.clear();
  for (auto& it : internal_timer_prop_int){ free(it); }
  for (auto& it : internal_timer_prop_double){ free(it); }
  for (auto& it : internal_timer_prop_double_int){ free(it); }
  for (auto& it : internal_timer_prop_char){ free(it); }
  internal_timer_prop_int.clear(); internal_timer_prop_double.clear(); internal_timer_prop_double_int.clear(); internal_timer_prop_char.clear(); internal_timer_prop_req.clear();
}


void forward_pass::initiate(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
                            bool is_sender, int partner1, int partner2){
  if (mode>=3){
    bool is_p2p = tracker.tag>=13;
    if (comm_pattern_table1.find(std::make_pair(tracker.tag,cm)) == comm_pattern_table1.end()){
      comm_pattern_table1.insert(std::make_pair(tracker.tag,cm));
    }
    if (is_p2p){
      if (p2p_table.find(std::make_pair(cm,partner1)) == p2p_table.end()){
        p2p_table.insert(std::make_pair(cm,partner1));
      }
    }
    comm_pattern_seq.push_back(std::make_pair(std::make_pair(tracker.tag,cm),partner1));
  }

  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  tracker.save_comp_time    = curtime - computation_timer;
  critical_path_costs[num_critical_path_measures-2] += tracker.save_comp_time;	// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += tracker.save_comp_time;	// update critical path runtime
  volume_costs[num_volume_measures-2]        += tracker.save_comp_time;		// update local computation time
  volume_costs[num_volume_measures-1]        += tracker.save_comp_time;		// update local runtime
  for (size_t i=0; i<breakdown_size; i++){ critical_path_costs[critical_path_costs_size-1-i] += tracker.save_comp_time; }
  if (mode>=2 && symbol_stack.size()>0){
    auto last_symbol_time = curtime - symbol_timers[symbol_stack.top()].start_timer.top();
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-2] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-2] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-1] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-2] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-2] += last_symbol_time;
  }

  int el_size,np,rank;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(cm, &np);
  MPI_Comm_rank(cm, &rank);
  tracker.last_nbytes = nbytes;
  tracker.last_cm = cm;
  tracker.last_p = np;
  tracker.last_is_sender = is_sender;
  tracker.last_partner1 = partner1;
  tracker.last_partner2 = partner2;

  volatile double init_time = MPI_Wtime();
  if (tracker.last_partner1 == -1){ PMPI_Barrier(cm); }
  else {
    char sbuf='e'; char rbuf='d';
    if ((tracker.tag < 14) || (tracker.tag > 15)){
      if ((is_sender) && (rank != tracker.last_partner1)){
        PMPI_Ssend(&sbuf, 1, MPI_CHAR, tracker.last_partner1, internal_tag3, cm);
      } else if ((!is_sender) && (rank != tracker.last_partner1)){
        PMPI_Recv(&rbuf, 1, MPI_CHAR, tracker.last_partner1, internal_tag3, cm, MPI_STATUS_IGNORE);
      }
    } else{
        PMPI_Sendrecv(&sbuf, 1, MPI_CHAR, tracker.last_partner1, internal_tag3, &rbuf, 1, MPI_CHAR, tracker.last_partner1, internal_tag3, cm, MPI_STATUS_IGNORE);
    }
/*    if ((is_sender) && (tracker.last_partner2 != -1)){
      PMPI_Recv(&rbuf, 1, MPI_CHAR, tracker.last_partner2, internal_tag3, cm, MPI_STATUS_IGNORE);
    }
*/
  }
  tracker.last_barrier_time = MPI_Wtime() - init_time;
  //TODO: We want to sbtract out the idle time of the path-root along the execution-time cp so that it appears as this path has no idle time.
  //      Ideally we would do this for the last process to enter this barrier (which would always determine the execution-time cp
  //      This would require use of the DOUBLE_INT 2-stage propagation mechanism, which is where I want to go next. I think its more efficient as well. See notes.
  for (size_t i=0; i<breakdown_size; i++){ critical_path_costs[critical_path_costs_size-1-i-breakdown_size] += tracker.last_barrier_time; }
  std::vector<int> counts(np,1); std::vector<int> disp(np,0);// could be moved as global variables
  if (tracker.tag>=9 && tracker.tag<=12) for (int i=1; i<np; i++) disp[i]=disp[i-1]+1;
  tracker.last_start_time = MPI_Wtime();

  switch (tracker.tag){
    case 0:
      PMPI_Barrier(cm);
      break;
    case 1:
      PMPI_Bcast(&synch_pad_send[0], 1, MPI_CHAR, 0, cm);// arbitrary root 0
      break;
    case 2:
      PMPI_Reduce(&synch_pad_send[0], &synch_pad_recv[0], 1, MPI_CHAR, MPI_MAX, 0, cm);// arbitrary root 0
      break;
    case 3:
      PMPI_Allreduce(MPI_IN_PLACE, &synch_pad_send[0], 1, MPI_CHAR, MPI_MAX, cm);
      break;
    case 4:
      PMPI_Gather(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, 0, cm);// arbitrary root 0
      break;
    case 5:
      PMPI_Allgather(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, cm);
      break;
    case 6:
      PMPI_Scatter(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, 0, cm);// arbitrary root 0
      break;
    case 7:
      PMPI_Reduce_scatter(&synch_pad_send[0], &synch_pad_recv[0], &counts[0], MPI_CHAR, MPI_MAX, cm);
      break;
    case 8:
      PMPI_Alltoall(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, cm);
      break;
    case 9:
      PMPI_Gatherv(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], &counts[0], &disp[0], MPI_CHAR, 0, cm);// arbitrary root 0
      break;
    case 10:
      PMPI_Allgatherv(&synch_pad_send[0], 1, MPI_CHAR, &synch_pad_recv[0], &counts[0], &disp[0], MPI_CHAR, cm);
      break;
    case 11:
      PMPI_Scatterv(&synch_pad_send[0], &counts[0], &disp[0], MPI_CHAR, &synch_pad_recv[0], 1, MPI_CHAR, 0, cm);// arbitrary root 0
      break;
    case 12:
      PMPI_Alltoallv(&synch_pad_send[0], &counts[0], &disp[0], MPI_CHAR, &synch_pad_recv[0], &counts[0], &disp[0], MPI_CHAR, cm);
      break;
    case 13:
      PMPI_Sendrecv(&synch_pad_send[0], 1, MPI_CHAR, partner1, internal_tag, &synch_pad_recv[0], 1, MPI_CHAR, partner2, internal_tag, cm, MPI_STATUS_IGNORE);
      break;
    case 14:
      PMPI_Sendrecv_replace(&synch_pad_send[0], 1, MPI_CHAR, partner1, internal_tag, partner2, internal_tag, cm, MPI_STATUS_IGNORE);
      break;
    case 15:
      PMPI_Ssend(&synch_pad_send[0], 1, MPI_CHAR, partner1, internal_tag, cm);
      break;
    case 16:
      PMPI_Send(&synch_pad_send[0], 1, MPI_CHAR, partner1, internal_tag, cm);
      break;
    case 17:
      PMPI_Recv(&synch_pad_recv[0], 1, MPI_CHAR, partner1, internal_tag, cm, MPI_STATUS_IGNORE);
      break;
  }

  // start synchronization timer for communication routine
  volatile double synchTime = MPI_Wtime();
  tracker.last_synch_time = synchTime-tracker.last_start_time;
  // start communication timer for communication routine
  tracker.last_start_time = MPI_Wtime();
}

// Used only for p2p communication. All blocking collectives use sychronous protocol
void forward_pass::complete(blocking& tracker){
  volatile double new_time = MPI_Wtime();
  volatile double dt = new_time - tracker.last_start_time;	// complete communication time
  double datamvt_time = std::max(0.,(dt-tracker.last_synch_time));
  std::pair<double,double> dcost_bsp    = tracker.cost_func_bsp(tracker.last_nbytes, tracker.last_p);
  std::pair<double,double> dcost_alphabeta_butterfly = tracker.cost_func_alphabeta_butterfly(tracker.last_nbytes, tracker.last_p);
  std::vector<std::pair<double,double>> dcosts = {dcost_bsp,dcost_alphabeta_butterfly};

  if (mode>=1){
    *tracker.my_synch_time   += tracker.last_synch_time;
    *tracker.my_datamvt_time += datamvt_time;
    *tracker.my_comm_time    += dt;
    int save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]=='1'){
        *(tracker.my_msg_count+save) += dcosts[j].first;
        *(tracker.my_wrd_count+save) += dcosts[j].second;
        save++;
      }
    }
    for (size_t i=0; i<breakdown_size; i++){
      *(tracker.critical_path_synch_time+i)   += tracker.last_synch_time;
      *(tracker.critical_path_datamvt_time+i) += datamvt_time;
      *(tracker.critical_path_comm_time+i)    += dt;
    }
    save=0;
    for (int j=0; j<cost_models.size(); j++){
      for (size_t i=0; i<breakdown_size; i++){
        if (cost_models[j]=='1'){
          *(tracker.critical_path_msg_count+save*breakdown_size+i) += dcosts[j].first;
          *(tracker.critical_path_wrd_count+save*breakdown_size+i) += dcosts[j].second;
        }
      }
      save++;
    }
  }
  if (mode>=2 && symbol_stack.size()>0){
    // update all communication-related measures for the top symbol in stack
    size_t save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]=='1'){
        symbol_timers[symbol_stack.top()].cp_exclusive_measure[save] += dcosts[j].second;
        symbol_timers[symbol_stack.top()].cp_exclusive_measure[cost_models.size()+save] += dcosts[j].first;
        symbol_timers[symbol_stack.top()].pp_exclusive_measure[save] += dcosts[j].second;
        symbol_timers[symbol_stack.top()].pp_exclusive_measure[cost_models.size()+save] += dcosts[j].first;
        *symbol_timers[symbol_stack.top()].cp_excl_measure[save] += dcosts[j].second;
        *symbol_timers[symbol_stack.top()].cp_excl_measure[cost_models.size()+save] += dcosts[j].first;
        *symbol_timers[symbol_stack.top()].pp_excl_measure[save] += dcosts[j].second;
        *symbol_timers[symbol_stack.top()].pp_excl_measure[cost_models.size()+save] += dcosts[j].first;
      }
      save++;
    }
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-5] += dt;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-4] += tracker.last_synch_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-3] += datamvt_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-1] += dt;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-6] += tracker.last_barrier_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-6] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-5] += dt;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-4] += tracker.last_synch_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-3] += datamvt_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += (dt+tracker.last_barrier_time);
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-5] += dt;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-4] += tracker.last_synch_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-3] += datamvt_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-1] += dt;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-6] += tracker.last_barrier_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-6] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-5] += dt;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-4] += tracker.last_synch_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-3] += datamvt_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += (dt+tracker.last_barrier_time);
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
  }

  int save=0;
  for (int j=0; j<cost_models.size(); j++){
    if (cost_models[j]=='1'){
      critical_path_costs[save]                 += dcosts[j].second;		// update critical path estimated communication cost
      critical_path_costs[cost_model_size+save] += dcosts[j].first;		// update critical path estimated synchronization cost
      save++;
    }
  }
  critical_path_costs[num_critical_path_measures-5] += dt;				// update critical path communication time (for what this process has seen thus far)
  critical_path_costs[num_critical_path_measures-4] += tracker.last_synch_time;	// update critical path synchronization time
  critical_path_costs[num_critical_path_measures-3] += datamvt_time;		// update critical path data mvt time
  critical_path_costs[num_critical_path_measures-1] += dt;				// update critical path runtime

  save=0;
  for (int j=0; j<cost_models.size(); j++){
    if (cost_models[j]=='1'){
      volume_costs[save]                 += dcosts[j].second;		// update local estimated communication cost
      volume_costs[cost_model_size+save] += dcosts[j].first;		// update local estimated synchronization cost
      save++;
    }
  }
  volume_costs[num_volume_measures-6] += tracker.last_barrier_time;	// update local barrier/idle time
  volume_costs[num_volume_measures-5] += dt;				// update local communication time (not volume until after the completion of the program)
  volume_costs[num_volume_measures-4] += tracker.last_synch_time;		// update local synchronization time
  volume_costs[num_volume_measures-3] += datamvt_time;			// update local data mvt time
  volume_costs[num_volume_measures-1] += (tracker.last_barrier_time+dt);	// update local runtime with idle time and comm time
  volume_costs[num_volume_measures-6] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);

  if (mode>=2 && symbol_stack.size()>0){
    // Special handling of excessively large idle time caused by suspected tool interference
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-6] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-6] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
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
  propagate(tracker.last_cm, tracker.tag, tracker.last_is_sender, tracker.last_partner1, tracker.last_partner2);
  // Prepare to leave interception and re-enter user code
  tracker.last_start_time = MPI_Wtime();
  computation_timer = tracker.last_start_time;
  if (mode>=2 && symbol_stack.size()>0){ symbol_timers[symbol_stack.top()].start_timer.top() = tracker.last_start_time; }
}

// Called by both nonblocking p2p and nonblocking collectives
void forward_pass::initiate(nonblocking& tracker, volatile double curtime, volatile double itime, int64_t nelem,
                            MPI_Datatype t, MPI_Comm cm, MPI_Request* request, bool is_sender, int partner){

  if (mode>=3){
    bool is_p2p = tracker.tag<20;
    if (comm_pattern_table1.find(std::make_pair(tracker.tag,cm)) == comm_pattern_table1.end()){
      comm_pattern_table1.insert(std::make_pair(tracker.tag,cm));
    }
    if (is_p2p){
      if (p2p_table.find(std::make_pair(cm,partner)) == p2p_table.end()){
        p2p_table.insert(std::make_pair(cm,partner));
      }
    }
    comm_pattern_seq.push_back(std::make_pair(std::make_pair(tracker.tag,cm),partner));
  }

  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  tracker.save_comp_time = curtime - computation_timer + itime;
  critical_path_costs[num_critical_path_measures-2] += tracker.save_comp_time;		// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += tracker.save_comp_time;		// update critical path runtime
  volume_costs[num_volume_measures-2]        += tracker.save_comp_time;		// update local computation time
  volume_costs[num_volume_measures-1]        += tracker.save_comp_time;		// update local runtime
  for (size_t i=0; i<breakdown_size; i++){ critical_path_costs[critical_path_costs_size-1-i] += tracker.save_comp_time; }
  if (mode>=2 && symbol_stack.size()>0){
    assert(symbol_stack.size()>0);
    assert(symbol_timers[symbol_stack.top()].start_timer.size()>0);
    tracker.save_time = curtime - symbol_timers[symbol_stack.top()].start_timer.top()+itime;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-1] += tracker.save_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-2] += tracker.save_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += tracker.save_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-2] += tracker.save_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-1] += tracker.save_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-2] += tracker.save_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-1] += tracker.save_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-2] += tracker.save_time;
  }

  int el_size,p;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(cm, &p);
  int rank; MPI_Comm_rank(cm, &rank);

  internal_comm_info[*request] = is_sender;
  internal_comm_comm[*request] = std::make_pair(cm,partner);
  internal_comm_data[*request] = std::make_pair((double)nbytes,(double)p);
  internal_comm_track[*request] = &tracker;

  tracker.last_start_time = MPI_Wtime();
  computation_timer = tracker.last_start_time;
  if (mode>=2 && symbol_stack.size()>0){ symbol_timers[symbol_stack.top()].start_timer.top() = tracker.last_start_time; }
}

void forward_pass::complete(nonblocking& tracker, MPI_Request* request, double comp_time, double comm_time){
  auto comm_info_it = internal_comm_info.find(*request);
  auto comm_comm_it = internal_comm_comm.find(*request);
  auto comm_data_it = internal_comm_data.find(*request);
  auto comm_track_it = internal_comm_track.find(*request);
  assert(comm_info_it != internal_comm_info.end());
  assert(comm_comm_it != internal_comm_comm.end());
  assert(comm_data_it != internal_comm_data.end());
  assert(comm_track_it != internal_comm_track.end());

  tracker.last_is_sender = comm_info_it->second;
  tracker.last_cm = comm_comm_it->second.first;
  tracker.last_partner1 = comm_comm_it->second.second;
  tracker.last_partner2 = -1;
  tracker.last_nbytes = comm_data_it->second.first;
  tracker.last_p = comm_data_it->second.second;
  tracker.last_synch_time=0;
  
  // Both sender and receiver will now update its critical path with the data from the communication
  std::pair<double,double> dcost_bsp  = tracker.cost_func_bsp(tracker.last_nbytes,tracker.last_p);
  std::pair<double,double> dcost_alphabeta_butterfly = tracker.cost_func_alphabeta_butterfly(tracker.last_nbytes,tracker.last_p);
  if ((tracker.tag<20) && (wait_id)) dcost_bsp.first=1.;	// this is usually zero, but we force it to be 1 in special circumstances (for nonblocking p2p with wait_id one)
  std::vector<std::pair<double,double>> dcosts = {dcost_bsp,dcost_alphabeta_butterfly};

  if (mode >= 1){
    *tracker.my_synch_time   += 0;			// Nonblocking routines will have no synchronization time component
    *tracker.my_datamvt_time += comm_time;
    *tracker.my_comm_time    += comm_time;
    int save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]=='1'){
        *(tracker.my_msg_count+save) += dcosts[j].first;
        *(tracker.my_wrd_count+save) += dcosts[j].second;
        save++;
      }
    }
    save=0;
    for (size_t i=0; i<breakdown_size; i++){
      *(tracker.critical_path_synch_time+i)   += tracker.last_synch_time;
      *(tracker.critical_path_datamvt_time+i) += comm_time;
      *(tracker.critical_path_comm_time+i)    += comm_time;
    }
    for (int j=0; j<cost_models.size(); j++){
      for (size_t i=0; i<breakdown_size; i++){
        if (cost_models[j]=='1'){
          *(tracker.critical_path_msg_count+save*breakdown_size+i) += dcosts[j].first;
          *(tracker.critical_path_wrd_count+save*breakdown_size+i) += dcosts[j].second;
        }
      }
      save++;
    }
  }
  if (mode>=2 && symbol_stack.size()>0){
    size_t save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]=='1'){
        symbol_timers[symbol_stack.top()].cp_exclusive_measure[save] += dcosts[j].second;
        symbol_timers[symbol_stack.top()].cp_exclusive_measure[cost_models.size()+save] += dcosts[j].first;
        symbol_timers[symbol_stack.top()].pp_exclusive_measure[save] += dcosts[j].second;
        symbol_timers[symbol_stack.top()].pp_exclusive_measure[cost_models.size()+save] += dcosts[j].first;
        *symbol_timers[symbol_stack.top()].cp_excl_measure[save] += dcosts[j].second;
        *symbol_timers[symbol_stack.top()].cp_excl_measure[cost_models.size()+save] += dcosts[j].first;
        *symbol_timers[symbol_stack.top()].pp_excl_measure[save] += dcosts[j].second;
        *symbol_timers[symbol_stack.top()].pp_excl_measure[cost_models.size()+save] += dcosts[j].first;
      }
      save++;
    }
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-5] += comm_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-4] += 0;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-3] += comm_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-2] += comp_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-1] += (comp_time+comm_time);
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-5] += comm_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-4] += 0;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-3] += comm_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-2] += comp_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += (comp_time+comm_time);
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-5] += comm_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-4] += 0;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-3] += comm_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-2] += comp_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-1] += (comp_time+comm_time);
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-5] += comm_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-4] += 0;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-3] += comm_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-2] += comp_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += (comp_time+comm_time);
  }

  int save=0;
  for (int j=0; j<cost_models.size(); j++){
    if (cost_models[j]=='1'){
      critical_path_costs[save]                 += dcosts[j].second;		// update critical path estimated communication cost
      critical_path_costs[cost_model_size+save] += dcosts[j].first;		// update critical path estimated synchronization cost
      save++;
    }
  }
  critical_path_costs[num_critical_path_measures-5] += comm_time;				// update critical path communication time (for what this process has seen thus far)
  critical_path_costs[num_critical_path_measures-4] += 0.;				// update critical path synchronization time
  critical_path_costs[num_critical_path_measures-3] += comm_time;				// update critical path runtime
  critical_path_costs[num_critical_path_measures-2] += comp_time;				// update critical path runtime
  critical_path_costs[num_critical_path_measures-1] += comp_time+comm_time;				// update critical path runtime
  for (size_t i=0; i<breakdown_size; i++){
    critical_path_costs[critical_path_costs_size-1-i] += comp_time;
  }

  save=0;
  for (int j=0; j<cost_models.size(); j++){
    if (cost_models[j]=='1'){
      volume_costs[save]                 += dcosts[j].second;		// update local estimated communication cost
      volume_costs[cost_model_size+save] += dcosts[j].first;		// update local estimated synchronization cost
      save++;
    }
  }
  volume_costs[num_volume_measures-5] += comm_time;				// update local communication time (not volume until after the completion of the program)
  volume_costs[num_volume_measures-4] += 0.;		// update local synchronization time
  volume_costs[num_volume_measures-3] += comm_time;			// update local data mvt time
  volume_costs[num_volume_measures-2] += comp_time;				// update local runtime
  volume_costs[num_volume_measures-1] += comp_time+comm_time;				// update local runtime

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

  if (mode>=1){
    propagate(tracker.last_cm,tracker.tag,tracker.last_is_sender,tracker.last_partner1,tracker.last_partner2);
  }
  internal_comm_info.erase(*request);
  internal_comm_comm.erase(*request);
  internal_comm_data.erase(*request);
  internal_comm_track.erase(*request);

  tracker.last_start_time = MPI_Wtime();
}

void forward_pass::wait(double curtime, MPI_Request* request, MPI_Status* status){
  double save_comp_time = curtime - computation_timer;
  auto comm_track_it = internal_comm_track.find(*request);
  assert(comm_track_it != internal_comm_track.end());
  auto comm_info_it = internal_comm_info.find(*request);
  auto comm_comm_it = internal_comm_comm.find(*request);
  MPI_Request save_request = comm_info_it->first;
  int comm_rank; MPI_Comm_rank(comm_comm_it->second.first,&comm_rank);
  if (comm_info_it->second && comm_comm_it->second.second != -1 && comm_rank != comm_comm_it->second.second){
    PMPI_Ssend(&barrier_pad_send[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag3, comm_comm_it->second.first);
    PMPI_Ssend(&synch_pad_send[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag, comm_comm_it->second.first);
  }
  else if (!comm_info_it->second && comm_comm_it->second.second != -1 && comm_rank != comm_comm_it->second.second){
    PMPI_Recv(&barrier_pad_recv[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag3, comm_comm_it->second.first, MPI_STATUS_IGNORE);
    PMPI_Recv(&synch_pad_recv[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag, comm_comm_it->second.first, MPI_STATUS_IGNORE);
  }
  volatile double last_start_time = MPI_Wtime();
  PMPI_Wait(request, status);
  curtime = MPI_Wtime(); double save_comm_time = curtime - last_start_time;
  complete(*comm_track_it->second, &save_request, save_comp_time, save_comm_time);
  complete_path_update();
  computation_timer = MPI_Wtime();
  if (mode>=2){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
}

void forward_pass::wait(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){
  bool success=false;
  for (int i=0; i<count; i++){
    if (*(array_of_requests+i) != MPI_REQUEST_NULL){
      MPI_Request* req = (MPI_Request*)array_of_requests+i;
      MPI_Status* stat=(MPI_Status*)status+i;
      double save_comp_time = (i==0 ? curtime - computation_timer : 0);// avoid over-counting.
      auto comm_track_it = internal_comm_track.find(*req);
      assert(comm_track_it != internal_comm_track.end());
      auto comm_info_it = internal_comm_info.find(*req);
      auto comm_comm_it = internal_comm_comm.find(*req);
      MPI_Request save_request = comm_info_it->first;
      int comm_rank; MPI_Comm_rank(comm_comm_it->second.first,&comm_rank);
      if (comm_info_it->second && comm_comm_it->second.second != -1 && comm_rank != comm_comm_it->second.second){
        PMPI_Ssend(&barrier_pad_send[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag3, comm_comm_it->second.first);
        PMPI_Ssend(&synch_pad_send[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag, comm_comm_it->second.first);
      }
      else if (!comm_info_it->second && comm_comm_it->second.second != -1 && comm_rank != comm_comm_it->second.second){
        PMPI_Recv(&barrier_pad_recv[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag3, comm_comm_it->second.first, MPI_STATUS_IGNORE);
        PMPI_Recv(&synch_pad_recv[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag, comm_comm_it->second.first, MPI_STATUS_IGNORE);
      }
      volatile double last_start_time = MPI_Wtime();
      PMPI_Wait(req, stat);
      curtime = MPI_Wtime(); double save_comm_time = curtime - last_start_time;
      complete(*comm_track_it->second, &save_request, save_comp_time, save_comm_time);
      complete_path_update();
      computation_timer = MPI_Wtime();
      if (mode>=2){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
      *indx=i;
      success=true;
      break;
    }
  }
  if (!success) { *indx=MPI_UNDEFINED; }
}

void forward_pass::wait(int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){
  std::vector<MPI_Request> pt(count); for (int i=0;i<count;i++){pt[i]=(array_of_requests)[i];}
  volatile double last_start_time = MPI_Wtime();
  PMPI_Waitany(count, array_of_requests, indx, status);
  volatile double curtime = MPI_Wtime(); double save_comm_time = curtime - last_start_time;
  MPI_Request request = pt[*indx];
  auto comm_track_it = internal_comm_track.find(request);
  assert(comm_track_it != internal_comm_track.end());
  complete(*comm_track_it->second, &request, waitall_comp_time, save_comm_time);
  waitall_comp_time=0;
  if (!waitall_id){ complete_path_update(); }
}

void forward_pass::wait(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[],
                        MPI_Status array_of_statuses[]){
  for (int i=0; i<incount; i++){
    if (*(array_of_requests+i) != MPI_REQUEST_NULL){
      MPI_Request* req = array_of_requests+i; MPI_Status* stat=array_of_statuses+i;
      double save_comp_time = i==0 ? curtime - computation_timer : 0;
      auto comm_track_it = internal_comm_track.find(*req);
      assert(comm_track_it != internal_comm_track.end());
      auto comm_info_it = internal_comm_info.find(*req);
      auto comm_comm_it = internal_comm_comm.find(*req);
      MPI_Request save_request = comm_info_it->first;
      int comm_rank; MPI_Comm_rank(comm_comm_it->second.first,&comm_rank);
      if (comm_info_it->second && comm_comm_it->second.second != -1 && comm_rank != comm_comm_it->second.second){
        PMPI_Ssend(&barrier_pad_send[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag3, comm_comm_it->second.first);
        PMPI_Ssend(&synch_pad_send[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag, comm_comm_it->second.first);
      }
      else if (!comm_info_it->second && comm_comm_it->second.second != -1 && comm_rank != comm_comm_it->second.second){
        PMPI_Recv(&barrier_pad_recv[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag3, comm_comm_it->second.first, MPI_STATUS_IGNORE);
        PMPI_Recv(&synch_pad_recv[0], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag, comm_comm_it->second.first, MPI_STATUS_IGNORE);
      }
      volatile double last_start_time = MPI_Wtime();
      PMPI_Wait(req, stat);
      curtime = MPI_Wtime(); double save_comm_time = curtime - last_start_time;
      complete(*comm_track_it->second, &save_request, save_comp_time, save_comm_time);
      complete_path_update();
      computation_timer = MPI_Wtime();
      if (mode>=2){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
      array_of_indices[0]=i;
      *outcount=1;
    }
  }
}

void forward_pass::wait(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
  waitall_comp_time = curtime - computation_timer;
  wait_id=true;
  waitall_id=true;
  std::vector<MPI_Request> internal_requests(2*int(count),MPI_REQUEST_NULL);
  if (count > barrier_pad_send.size()){
    barrier_pad_send.resize(count);
    barrier_pad_recv.resize(count);
    synch_pad_send.resize(count);
    synch_pad_recv.resize(count);
  }
  for (int i=0; i<count; i++){
    auto comm_info_it = internal_comm_info.find(*(array_of_requests+i));
    assert(comm_info_it != internal_comm_info.end());
    auto comm_comm_it = internal_comm_comm.find(*(array_of_requests+i));
    assert(comm_comm_it != internal_comm_comm.end());
    if (comm_info_it->second && comm_comm_it->second.second != -1){
      PMPI_Isend(&barrier_pad_send[i], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag3,
        comm_comm_it->second.first, &internal_requests[2*i]);
      PMPI_Isend(&synch_pad_send[i], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag,
        comm_comm_it->second.first, &internal_requests[2*i+1]);
    }
    else if (!comm_info_it->second && comm_comm_it->second.second != -1){
      PMPI_Irecv(&barrier_pad_recv[i], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag3,
        comm_comm_it->second.first, &internal_requests[2*i]);
      PMPI_Irecv(&synch_pad_recv[i], 1, MPI_CHAR, comm_comm_it->second.second, internal_tag,
        comm_comm_it->second.first, &internal_requests[2*i+1]);
    }
  }
  PMPI_Waitall(internal_requests.size(), &internal_requests[0], MPI_STATUSES_IGNORE);
  int indx; MPI_Status stat;
  for (int i=0; i<count; i++){
    wait(count, array_of_requests, &indx, &stat);
    if (i==0){wait_id=false;}
    if ((MPI_Status*)array_of_statuses != (MPI_Status*)MPI_STATUSES_IGNORE) ((MPI_Status*)array_of_statuses)[indx] = stat;
  }
  wait_id=true;
  complete_path_update();
  waitall_id=false;
  computation_timer = MPI_Wtime();
  if (mode>=2){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
}

void forward_pass::initiate_timers(int rank, int tag, MPI_Comm cm, int partner1, int partner2){
  assert(mode>=2);
  MPI_Request internal_request[10];
  int* send_envelope1 = nullptr; int* send_envelope2 = nullptr; double* send_envelope3 = nullptr; double* send_envelope4 = nullptr; char* send_envelope5 = nullptr;
  int* recv_envelope1 = nullptr; int* recv_envelope2 = nullptr; double* recv_envelope3 = nullptr; double* recv_envelope4 = nullptr; char* recv_envelope5 = nullptr;
  int ftimer_size = symbol_timers.size();
  int num_chars = 0;
  for (int i=0; i<ftimer_size; i++) { num_chars += symbol_order[i].size(); }
  send_envelope1 = (int*)malloc(sizeof(int)); *send_envelope1 = ftimer_size;
  send_envelope2 = (int*)malloc(sizeof(int)*(ftimer_size));
  send_envelope3 = (double*)malloc(sizeof(double)*(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size);
  send_envelope4 = (double*)malloc(sizeof(double)*2*ftimer_size*num_critical_path_measures);
  send_envelope5 = (char*)malloc(sizeof(char)*num_chars);
  int symbol_offset = 0;
  for (auto i=0; i<ftimer_size; i++){
    send_envelope2[i] = symbol_order[i].size();
    for (auto j=0; j<symbol_order[i].size(); j++){
      send_envelope5[symbol_offset+j] = symbol_order[i][j];
    }
    for (auto j=0; j<num_critical_path_measures; j++){
      send_envelope4[i*num_critical_path_measures+j] = symbol_timers[symbol_order[i]].cp_exclusive_contributions[j];
    }
    for (auto k=0; k<num_critical_path_measures; k++){
      send_envelope4[ftimer_size*num_critical_path_measures + i*num_critical_path_measures+k] = symbol_timers[symbol_order[i]].cp_exclusive_measure[k];
    }
    symbol_offset += symbol_order[i].size();
  }
  for (int i=0; i<(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size; i++){ send_envelope3[i] = symbol_timer_pad_local_cp[i]; }
  PMPI_Isend(&send_envelope1[0],1,MPI_INT,partner1,internal_tag1,cm,&internal_request[0]);
  PMPI_Isend(&send_envelope2[0],ftimer_size,MPI_INT,partner1,internal_tag2,cm,&internal_request[1]);
  PMPI_Isend(&send_envelope3[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,partner1,internal_tag3,cm,&internal_request[2]);
  PMPI_Isend(&send_envelope4[0],2*ftimer_size*num_critical_path_measures,MPI_DOUBLE,partner1,internal_tag4,cm,&internal_request[3]);
  PMPI_Isend(&send_envelope5[0],symbol_offset,MPI_CHAR,partner1,internal_tag5,cm,&internal_request[4]);

  recv_envelope1 = (int*)malloc(sizeof(int));
  recv_envelope2 = (int*)malloc(sizeof(int)*(max_num_symbols));
  recv_envelope3 = (double*)malloc(sizeof(double)*(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols);
  recv_envelope4 = (double*)malloc(sizeof(double)*2*max_num_symbols*num_critical_path_measures);
  recv_envelope5 = (char*)malloc(sizeof(char)*max_timer_name_length*max_num_symbols);
  PMPI_Irecv(recv_envelope1,1,MPI_INT,partner1,internal_tag1,cm,&internal_request[5]);
  PMPI_Irecv(recv_envelope2,max_num_symbols,MPI_INT,partner1,internal_tag2,cm,&internal_request[6]);
  PMPI_Irecv(recv_envelope3,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols,MPI_DOUBLE,partner1,internal_tag3,cm,&internal_request[7]);
  PMPI_Irecv(recv_envelope4,2*max_num_symbols*num_critical_path_measures,MPI_DOUBLE,partner1,internal_tag4,cm,&internal_request[8]);
  PMPI_Irecv(recv_envelope5,max_timer_name_length*max_num_symbols,MPI_CHAR,partner1,internal_tag5,cm,&internal_request[9]);

  for (int i=0; i<10; i++) { internal_timer_prop_req.push_back(internal_request[i]); }
  internal_timer_prop_int.push_back(send_envelope1); internal_timer_prop_int.push_back(send_envelope2);
  internal_timer_prop_int.push_back(recv_envelope1); internal_timer_prop_int.push_back(recv_envelope2);
  internal_timer_prop_double.push_back(send_envelope3); internal_timer_prop_double.push_back(send_envelope4);
  internal_timer_prop_double.push_back(recv_envelope3); internal_timer_prop_double.push_back(recv_envelope4);
  internal_timer_prop_char.push_back(send_envelope5); internal_timer_prop_char.push_back(recv_envelope5);
}

/*
 Its important to note here that a blocking p2p call will already know whether its the cp root or not, regardless of whether its partner used a nonblocking p2p routine.
   But, because that potential nonblocking partner posted both sends and recvs, the blocking partner also has to do so, even if its partner (unknown to him) used a blocking p2p routine
*/
void forward_pass::propagate_timers(int rank, int tag, MPI_Comm cm, int partner1, int partner2){
  assert(mode>=2);
  int critical_path_runtime_root_rank = timer_info_receiver[num_critical_path_measures-1].second;
  int ftimer_size_cp = 0; int ftimer_size_ncp = 0;
  if (rank==critical_path_runtime_root_rank){ ftimer_size_cp = symbol_timers.size(); }
  else                                      { ftimer_size_ncp = symbol_timers.size(); }
  if (partner1 == -1){ PMPI_Allreduce(MPI_IN_PLACE,&ftimer_size_cp,1,MPI_INT,MPI_SUM,cm); }
  else{
    if (rank==critical_path_runtime_root_rank){ PMPI_Send(&ftimer_size_cp,1,MPI_INT,partner1,internal_tag1,cm);
                                                PMPI_Recv(&ftimer_size_ncp,1,MPI_INT,partner1,internal_tag1,cm,MPI_STATUS_IGNORE);
                                              }
    else { PMPI_Recv(&ftimer_size_cp,1,MPI_INT,critical_path_runtime_root_rank,internal_tag1,cm,MPI_STATUS_IGNORE);
           PMPI_Send(&ftimer_size_ncp,1,MPI_INT,critical_path_runtime_root_rank,internal_tag1,cm);
         }
  }
  for (auto i=0; i<symbol_len_pad_cp.size(); i++){ symbol_len_pad_cp[i] = 0.; symbol_len_pad_ncp[i] = 0.; }
  std::vector<double> cp_data(2*ftimer_size_cp*num_critical_path_measures,0);
  std::vector<double> ncp_data(2*ftimer_size_ncp*num_critical_path_measures,0);
  if (rank==critical_path_runtime_root_rank){
    int symbol_offset = 0;
    for (auto i=0; i<symbol_timers.size(); i++){
      symbol_len_pad_cp[i] = symbol_order[i].size();
      for (auto j=0; j<symbol_len_pad_cp[i]; j++){
        symbol_pad_cp[symbol_offset+j] = symbol_order[i][j];
      }
      for (auto j=0; j<num_critical_path_measures; j++){
        cp_data[i*num_critical_path_measures+j] = symbol_timers[symbol_order[i]].cp_exclusive_contributions[j];
      }
      for (auto k=0; k<num_critical_path_measures; k++){
        cp_data[ftimer_size_cp*num_critical_path_measures + i*num_critical_path_measures+k] = symbol_timers[symbol_order[i]].cp_exclusive_measure[k];
      }
      symbol_offset += symbol_len_pad_cp[i];
    }
  }
  else{
    int symbol_offset = 0;
    for (auto i=0; i<symbol_timers.size(); i++){
      symbol_len_pad_ncp[i] = symbol_order[i].size();
      for (auto j=0; j<symbol_len_pad_ncp[i]; j++){
        symbol_pad_ncp[symbol_offset+j] = symbol_order[i][j];
      }
      for (auto j=0; j<num_critical_path_measures; j++){
        ncp_data[i*num_critical_path_measures+j] = symbol_timers[symbol_order[i]].cp_exclusive_contributions[j];
      }
      for (auto k=0; k<num_critical_path_measures; k++){
        ncp_data[ftimer_size_ncp*num_critical_path_measures + i*num_critical_path_measures+k] = symbol_timers[symbol_order[i]].cp_exclusive_measure[k];
      }
      symbol_offset += symbol_len_pad_ncp[i];
    }
  }
  if (partner1 == -1){ PMPI_Allreduce(MPI_IN_PLACE,&symbol_len_pad_cp[0],ftimer_size_cp,MPI_INT,MPI_SUM,cm); }
  else{
    if (rank==critical_path_runtime_root_rank){ PMPI_Send(&symbol_len_pad_cp[0],ftimer_size_cp,MPI_INT,partner1,internal_tag2,cm);
                                                PMPI_Recv(&symbol_len_pad_ncp[0],ftimer_size_ncp,MPI_INT,partner1,internal_tag2,cm,MPI_STATUS_IGNORE); }
    else{ PMPI_Recv(&symbol_len_pad_cp[0],ftimer_size_cp,MPI_INT,critical_path_runtime_root_rank,internal_tag2,cm,MPI_STATUS_IGNORE);
          PMPI_Send(&symbol_len_pad_ncp[0],ftimer_size_ncp,MPI_INT,critical_path_runtime_root_rank,internal_tag2,cm);
        }
  }
  int num_chars_cp = 0; int num_chars_ncp = 0;
  for (auto i=0; i<ftimer_size_cp; i++){ num_chars_cp += symbol_len_pad_cp[i]; }
  for (auto i=0; i<ftimer_size_ncp; i++){ num_chars_ncp += symbol_len_pad_ncp[i]; }
  if (rank == critical_path_runtime_root_rank){
    if (partner1 == -1){
      PMPI_Bcast(&symbol_timer_pad_local_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size_cp,MPI_DOUBLE,rank,cm);
      PMPI_Bcast(&cp_data[0],2*ftimer_size_cp*num_critical_path_measures,MPI_DOUBLE,rank,cm);
      PMPI_Bcast(&symbol_pad_cp[0],num_chars_cp,MPI_CHAR,rank,cm);
    }
    else{
      PMPI_Send(&symbol_timer_pad_local_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size_cp,MPI_DOUBLE,partner1,internal_tag3,cm);
      PMPI_Send(&cp_data[0],2*ftimer_size_cp*num_critical_path_measures,MPI_DOUBLE,partner1,internal_tag4,cm);
      PMPI_Send(&symbol_pad_cp[0],num_chars_cp,MPI_CHAR,partner1,internal_tag5,cm);
      PMPI_Recv(&symbol_timer_pad_global_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size_ncp,MPI_DOUBLE,partner1,internal_tag3,cm,MPI_STATUS_IGNORE);
      PMPI_Recv(&ncp_data[0],2*ftimer_size_ncp*num_critical_path_measures,MPI_DOUBLE,partner1,internal_tag4,cm,MPI_STATUS_IGNORE);
      PMPI_Recv(&symbol_pad_ncp[0],num_chars_ncp,MPI_CHAR,partner1,internal_tag5,cm,MPI_STATUS_IGNORE);
    }
  }
  else{
    if (partner1 == -1){
      PMPI_Bcast(&symbol_timer_pad_global_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size_cp,MPI_DOUBLE,critical_path_runtime_root_rank,cm);
      PMPI_Bcast(&cp_data[0],2*ftimer_size_cp*num_critical_path_measures,MPI_DOUBLE,critical_path_runtime_root_rank,cm);
      PMPI_Bcast(&symbol_pad_cp[0],num_chars_cp,MPI_CHAR,critical_path_runtime_root_rank,cm);
    }
    else{
      PMPI_Recv(&symbol_timer_pad_global_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size_cp,MPI_DOUBLE,critical_path_runtime_root_rank,internal_tag3,cm,MPI_STATUS_IGNORE);
      PMPI_Recv(&cp_data[0],2*ftimer_size_cp*num_critical_path_measures,MPI_DOUBLE,critical_path_runtime_root_rank,internal_tag4,cm,MPI_STATUS_IGNORE);
      PMPI_Recv(&symbol_pad_cp[0],num_chars_cp,MPI_CHAR,critical_path_runtime_root_rank,internal_tag5,cm,MPI_STATUS_IGNORE);
      PMPI_Send(&symbol_timer_pad_local_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size_ncp,MPI_DOUBLE,critical_path_runtime_root_rank,internal_tag3,cm);
      PMPI_Send(&cp_data[0],2*ftimer_size_ncp*num_critical_path_measures,MPI_DOUBLE,critical_path_runtime_root_rank,internal_tag4,cm);
      PMPI_Send(&symbol_pad_ncp[0],num_chars_ncp,MPI_CHAR,critical_path_runtime_root_rank,internal_tag5,cm);
    }
    if (rank != critical_path_runtime_root_rank){
      int symbol_offset = 0;
      for (int i=0; i<ftimer_size_cp; i++){
        auto reconstructed_symbol = std::string(symbol_pad_cp.begin()+symbol_offset,symbol_pad_cp.begin()+symbol_offset+symbol_len_pad_cp[i]);
        if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
          symbol_timers[reconstructed_symbol] = symbol_tracker(reconstructed_symbol);
          symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
        }
        *symbol_timers[reconstructed_symbol].cp_numcalls = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i];
        for (int j=0; j<num_critical_path_measures; j++){
          *symbol_timers[reconstructed_symbol].cp_incl_measure[j] = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i+2*j+1];
          *symbol_timers[reconstructed_symbol].cp_excl_measure[j] = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i+2*(j+1)];
        }
        for (int j=0; j<symbol_timers[reconstructed_symbol].cp_exclusive_measure.size(); j++){ symbol_timers[reconstructed_symbol].cp_exclusive_measure[j]=0.; }
        for (int j=0; j<num_critical_path_measures; j++){
          symbol_timers[reconstructed_symbol].cp_exclusive_contributions[j] = cp_data[i*num_critical_path_measures+j];
        }
        for (int k=0; k<num_critical_path_measures; k++){
          symbol_timers[reconstructed_symbol].cp_exclusive_measure[k] = cp_data[ftimer_size_cp*num_critical_path_measures+i*num_critical_path_measures+k];
        }
        symbol_timers[reconstructed_symbol].has_been_processed = true;
        symbol_offset += symbol_len_pad_cp[i];
      }
      // Now cycle through and find the symbols that were not processed and set their accumulated measures to 0
      for (auto& it : symbol_timers){
        if (it.second.has_been_processed){ it.second.has_been_processed = false; }
        else{
          *it.second.cp_numcalls = 0;
          for (int j=0; j<num_critical_path_measures; j++){
            *it.second.cp_incl_measure[j] = 0;
            *it.second.cp_excl_measure[j] = 0;
          }
        }
      }
    }
  }
}

void forward_pass::propagate(MPI_Comm cm, int tag, bool is_sender, int partner1, int partner2){
  int rank; MPI_Comm_rank(cm,&rank);
  if ((rank == partner1) || (rank == partner2)) { return; } 
  if (mode>=2){
    for (int i=0; i<num_critical_path_measures; i++){
      timer_info_sender[i].first = critical_path_costs[i];
      timer_info_sender[i].second = rank;
    }
    if (partner1 == -1){
      if (tag < 13){
        PMPI_Allreduce(&timer_info_sender[0].first, &timer_info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, MPI_MAXLOC, cm);
      } else { assert(0); }
    } else{
      if (tag < 18){
        PMPI_Sendrecv(&timer_info_sender[0].first, num_critical_path_measures, MPI_DOUBLE_INT, partner1, internal_tag,
                      &timer_info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, partner1, internal_tag, cm, MPI_STATUS_IGNORE);
      }
      else{
        MPI_Request req1,req2;
        double_int* send_pathdata = (double_int*)malloc(num_critical_path_measures*sizeof(double_int));
        double_int* recv_pathdata = (double_int*)malloc(num_critical_path_measures*sizeof(double_int));
        memcpy(&send_pathdata[0].first, &timer_info_sender[0].first, num_critical_path_measures*sizeof(double_int));
        PMPI_Isend(&send_pathdata[0].first, num_critical_path_measures, MPI_DOUBLE_INT, partner1, internal_tag, cm, &req1);
        PMPI_Irecv(&recv_pathdata[0].first, num_critical_path_measures, MPI_DOUBLE_INT, partner1, internal_tag, cm, &req2);
        internal_timer_prop_req.push_back(req1); internal_timer_prop_req.push_back(req2);
        internal_timer_prop_double_int.push_back(send_pathdata); internal_timer_prop_double_int.push_back(recv_pathdata);
      }
      if (tag < 18){
        for (int i=0; i<num_critical_path_measures; i++){
          if (timer_info_sender[i].first>timer_info_receiver[i].first){timer_info_receiver[i].second = rank;}
          else if (timer_info_sender[i].first==timer_info_receiver[i].first){
            if (timer_info_sender[i].second < timer_info_receiver[i].second){ timer_info_receiver[i].second = rank; }
            else { timer_info_receiver[i].second = partner1; }
          }
          timer_info_receiver[i].first = std::max(timer_info_sender[i].first, timer_info_receiver[i].first);
        }
      }
    }
  }
  // Exchange the tracked routine critical path data
  if (partner1 == -1){
    MPI_Op op; MPI_Op_create((MPI_User_function*) propagate_critical_path_op,0,&op);
    if (tag < 20){
      PMPI_Allreduce(MPI_IN_PLACE, &critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, op, cm);
      MPI_Op_free(&op);
      if (mode >= 2) { propagate_timers(rank,tag,cm,partner1,partner2==-1 ? partner1 : partner2); }
    } else{
      MPI_Request req1;
      double* local_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
      std::memcpy(local_path_data, &critical_path_costs[0], critical_path_costs.size()*sizeof(double));
      PMPI_Iallreduce(MPI_IN_PLACE,local_path_data,critical_path_costs.size(),MPI_DOUBLE,op,cm,&req1);
      //MPI_Op_free(&op);
      internal_comm_prop.push_back(std::make_pair(local_path_data,true));
      internal_comm_prop_req.push_back(req1);
      if (mode>=2) { assert(0); initiate_timers(rank,tag,cm,partner1,partner2==-1 ? partner1 : partner2); }
    }
  }
  else{
   if (tag < 18){
      //TODO: does this work for mode=1 with isend+recv? It must have been...
      PMPI_Sendrecv(&critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, partner1, internal_tag2, &new_cs[0], critical_path_costs.size(), MPI_DOUBLE, partner1, internal_tag2, cm, MPI_STATUS_IGNORE);
      update_critical_path(&new_cs[0],&critical_path_costs[0],critical_path_costs_size);
      if (mode>=2) { propagate_timers(rank,tag,cm,partner1,partner2==-1 ? partner1 : partner2); }
    }
    else{// only possible for Isend and Irecv
      MPI_Request req1,req2;
      if (is_sender){
        double* local_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
        std::memcpy(local_path_data, &critical_path_costs[0], critical_path_costs.size()*sizeof(double));
        //TODO: Can I keep sending out `critical_path_costs` or must I make copies and send that out?
        PMPI_Isend(local_path_data, critical_path_costs.size(), MPI_DOUBLE, partner1, internal_tag2, cm, &req1);
        double* remote_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
        PMPI_Irecv(remote_path_data, critical_path_costs.size(), MPI_DOUBLE, partner1, internal_tag2, cm, &req2);
        internal_comm_prop.push_back(std::make_pair(local_path_data,true));
        internal_comm_prop_req.push_back(req1);
        internal_comm_prop.push_back(std::make_pair(remote_path_data,false));
        internal_comm_prop_req.push_back(req2);
      }
      else{
        double* local_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
        std::memcpy(local_path_data, &critical_path_costs[0], critical_path_costs.size()*sizeof(double));
        double* remote_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
        PMPI_Irecv(remote_path_data, critical_path_costs.size(), MPI_DOUBLE, partner1, internal_tag2, cm, &req2);
        PMPI_Isend(local_path_data, critical_path_costs.size(), MPI_DOUBLE, partner1, internal_tag2, cm, &req1);
        internal_comm_prop.push_back(std::make_pair(local_path_data,true));
        internal_comm_prop_req.push_back(req1);
        internal_comm_prop.push_back(std::make_pair(remote_path_data,false));
        internal_comm_prop_req.push_back(req2);
      }
      if (mode>=2) { initiate_timers(rank,tag,cm,partner1,partner2==-1 ? partner1 : partner2); }
    }
  }
}

}
}
