#include "path.h"
#include "../container/symbol_tracker.h"
#include "../../replay/path/path.h"
#include "../../discretization/util/util.h"
#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

void path::exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm){
  int new_comm_size; MPI_Comm_size(newcomm,&new_comm_size);
  int old_comm_rank; MPI_Comm_rank(oldcomm,&old_comm_rank);
  std::vector<int> gathered_ranks(new_comm_size,0);
  PMPI_Allgather(&old_comm_rank,1,MPI_INT,&gathered_ranks[0],1,MPI_INT,newcomm);
  // Now we detect the "color" (really the stride) via iteration
  int color = new_comm_size>1 ? gathered_ranks[1]-gathered_ranks[0] : 0;
  communicator_map[newcomm] = std::make_pair(new_comm_size,color);
  computation_timer = MPI_Wtime();
}

bool path::initiate_comp(size_t id, volatile double curtime, double flop_count, int param1, int param2, int param3, int param4, int param5){
  // accumulate computation time
  double save_comp_time = curtime - computation_timer;

  // Special exit if no kernels are to be scheduled -- the whole point is to get a reading on the total time it takes, which is
  //   to be attained with timers outside of critter..
  if (schedule_kernels==0){ return false; }
  volatile double overhead_start_time = MPI_Wtime();

  //critical_path_costs[num_critical_path_measures-2] += save_comp_time;	// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += save_comp_time;	// update critical path runtime
  //volume_costs[num_volume_measures-2]        += save_comp_time;		// update local computation time
  volume_costs[num_volume_measures-1]        += save_comp_time;		// update local runtime

  bool schedule_decision = schedule_kernels==1 ? true : false;
  if (autotuning_mode>0){
    comp_pattern_key p_id_1(-1,id,flop_count,param1,param2,param3,param4,param5);// '-1' argument is arbitrary, does not influence overloaded operators
      if (!(comp_pattern_map.find(p_id_1) == comp_pattern_map.end())){
        schedule_decision = should_schedule(comp_pattern_map[p_id_1])==1;
      }
  }

  comp_intercept_overhead += MPI_Wtime() - overhead_start_time;
  // start compunication timer for compunication routine
  comp_start_time = MPI_Wtime();
  return schedule_decision;
}

void path::complete_comp(size_t id, double flop_count, int param1, int param2, int param3, int param4, int param5){
  volatile double comp_time = MPI_Wtime() - comp_start_time;	// complete computation time

  // Special exit if no kernels are to be scheduled -- the whole point is to get a reading on the total time it takes, which is
  //   to be attained with timers outside of critter..
  if (schedule_kernels==0){ return; }
  volatile double overhead_start_time = MPI_Wtime();

  if (autotuning_mode>0){
    comp_pattern_key p_id_1(active_patterns.size(),id,flop_count,param1,param2,param3,param4,param5);// 'active_patterns.size()' argument is arbitrary, does not influence overloaded operators
      if (comp_pattern_map.find(p_id_1) == comp_pattern_map.end()){
        active_comp_pattern_keys.emplace_back(p_id_1);
        active_patterns.emplace_back(pattern());
        comp_pattern_map[p_id_1] = pattern_key_id(true,active_comp_pattern_keys.size()-1,active_patterns.size()-1,false);
      }
      if (should_schedule(comp_pattern_map[p_id_1]) == 0){
        comp_time = get_estimate(comp_pattern_map[p_id_1],comp_pattern_param,flop_count);
      }
      update(comp_pattern_map[p_id_1],comp_pattern_param,comp_time,flop_count);
  }
  //critical_path_costs[num_critical_path_measures-5] += flop_count;
  //critical_path_costs[num_critical_path_measures-2] += comp_time;
  critical_path_costs[num_critical_path_measures-1] += comp_time;

  //volume_costs[num_volume_measures-6] += flop_count;
  //volume_costs[num_volume_measures-2] += comp_time;
  volume_costs[num_volume_measures-1] += comp_time;

  comp_intercept_overhead += MPI_Wtime() - overhead_start_time;
  computation_timer = MPI_Wtime();
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
  for (int i=0; i<num_critical_path_measures; i++){
    inout[i] = std::max(inout[i],in[i]);
  }
}

static void propagate_critical_path_op(double* in, double* inout, int* len, MPI_Datatype* dtype){
  update_critical_path(in,inout,static_cast<size_t>(*len));
}

bool path::initiate_comm(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                            bool is_sender, int partner1, int partner2){
  // Save and accumulate the computation time between last communication routine as both execution-time and computation time
  //   into both the execution-time critical path data structures and the per-process data structures.
  tracker.comp_time = curtime - computation_timer;
  int rank; MPI_Comm_rank(comm, &rank);

  // Special exit if no kernels are to be scheduled -- the whole point is to get a reading on the total time it takes, which is
  //   to be attained with timers outside of critter..
  if (schedule_kernels==0){ return false; }
  volatile double overhead_start_time = MPI_Wtime();

  // We consider usage of Sendrecv variants to forfeit usage of eager internal communication.
  // Note that the reason we can't force user Bsends to be 'true_eager_p2p' is because the corresponding Receives would be expecting internal communications
  bool true_eager_p2p = ((eager_p2p == 1) && (tracker.tag!=13) && (tracker.tag!=14));
  if (true_eager_p2p){
    MPI_Buffer_attach(&eager_pad[0],eager_pad.size());
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
  tracker.barrier_time=0.;// might get updated below

  bool schedule_decision = schedule_kernels==1 ? true : false;
  int schedule_decision_int; int schedule_decision_foreign_int;
  if (autotuning_mode>0){
    comm_pattern_key p_id_1(-1,tracker.tag,communicator_map[tracker.comm].first,communicator_map[tracker.comm].second,tracker.nbytes,(tracker.partner1 == -1 ? -1 : rank - tracker.partner1));
      if (!(comm_pattern_map.find(p_id_1) == comm_pattern_map.end())){
        schedule_decision = should_schedule_global(comm_pattern_map[p_id_1])==1;
      }
  }

  if (autotuning_mode==0 || schedule_decision==true){
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
    }
  }

  comm_intercept_overhead_stage1 += MPI_Wtime() - overhead_start_time;
  overhead_start_time = MPI_Wtime();

  //critical_path_costs[num_critical_path_measures-2] += tracker.comp_time;	// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += tracker.comp_time;	// update critical path runtime
  //volume_costs[num_volume_measures-2]        += tracker.comp_time;		// update local computation time
  volume_costs[num_volume_measures-1]        += tracker.comp_time;		// update local runtime

  if (autotuning_mode>0 && schedule_decision == true){
    comm_pattern_key p_id_1(-1,tracker.tag,communicator_map[tracker.comm].first,communicator_map[tracker.comm].second,tracker.nbytes,(tracker.partner1 == -1 ? -1 : rank - tracker.partner1));
    if (!(comm_pattern_map.find(p_id_1) == comm_pattern_map.end())){
      schedule_decision = should_schedule(comm_pattern_map[p_id_1])==1;
    }
    schedule_decision_int = (int)schedule_decision;
    if (tracker.partner1 == -1){
      PMPI_Allreduce(MPI_IN_PLACE, &schedule_decision_int, 1, MPI_INT, MPI_SUM, tracker.comm);
      schedule_decision = schedule_decision_int > 0;
    }
    else{
      if (true_eager_p2p){
        if (tracker.is_sender){
          PMPI_Bsend(&schedule_decision_int, 1, MPI_INT, tracker.partner1, internal_tag2, tracker.comm);
        }
        else{
          PMPI_Recv(&schedule_decision_foreign_int, 1, MPI_INT, tracker.partner1, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
          schedule_decision = (schedule_decision_foreign_int > 0 ? true : false);
        }
      }
      else{
        PMPI_Sendrecv(&schedule_decision_int, 1, MPI_INT, tracker.partner1, internal_tag2, &schedule_decision_foreign_int, 1,
                      MPI_INT, tracker.partner2, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
        schedule_decision = ((schedule_decision_int > 0) || (schedule_decision_foreign_int>0) ? true : false);
      }
      if (tracker.partner2 != tracker.partner1){
        // This if-statement will never be breached if 'true_eager_p2p'=true anyways.
        PMPI_Sendrecv(&schedule_decision_int, 1, MPI_INT, tracker.partner2, internal_tag2, &schedule_decision_foreign_int, 1,
                           MPI_INT, tracker.partner2, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
        schedule_decision = ((schedule_decision_int > 0) || (schedule_decision_foreign_int>0) ? true : false);
      }
    }
    // If local steady_state has been reached, and we find out the other processes have reached the same, then we can set global_steady_state=1
    //   and never have to use another internal collective to check again.
    if (!(comm_pattern_map.find(p_id_1) == comm_pattern_map.end())){
      set_schedule(comm_pattern_map[p_id_1],schedule_decision);
      if (autotuning_mode == 3){// prevents critical path analysis from changing global steady_state
        active_patterns[comm_pattern_map[p_id_1].val_index].global_steady_state=0;
      }
    }
  } else if (schedule_decision == false){
    // This call is merely to increment the num_non_propagated member
    comm_pattern_key p_id_1(-1,tracker.tag,communicator_map[tracker.comm].first,communicator_map[tracker.comm].second,tracker.nbytes,(tracker.partner1 == -1 ? -1 : rank - tracker.partner1));
    set_schedule(comm_pattern_map[p_id_1],schedule_decision);
  }

  comm_intercept_overhead_stage2 += MPI_Wtime() - overhead_start_time;
  // start communication timer for communication routine
  tracker.start_time = MPI_Wtime();
  return schedule_decision;
}

// Used only for p2p communication. All blocking collectives use sychronous protocol
void path::complete_comm(blocking& tracker, int recv_source){
  volatile double comm_time = MPI_Wtime() - tracker.start_time;	// complete communication time

  // Special exit if no kernels are to be scheduled -- the whole point is to get a reading on the total time it takes, which is
  //   to be attained with timers outside of critter..
  if (schedule_kernels==0){ return; }
  volatile double overhead_start_time = MPI_Wtime();

  int rank; MPI_Comm_rank(tracker.comm,&rank);
  bool true_eager_p2p = ((eager_p2p == 1) && (tracker.tag!=13) && (tracker.tag!=14));
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

  bool autotuning_special_bool = false;
  if (autotuning_mode>0){
    assert(communicator_map.find(tracker.comm) != communicator_map.end());
    comm_pattern_key p_id_1(active_patterns.size(),tracker.tag,communicator_map[tracker.comm].first,communicator_map[tracker.comm].second,tracker.nbytes,(tracker.partner1 == -1 ? -1 : rank - tracker.partner1));
      if (comm_pattern_map.find(p_id_1) == comm_pattern_map.end()){
        active_comm_pattern_keys.emplace_back(p_id_1);
        active_patterns.emplace_back(pattern());
        comm_pattern_map[p_id_1] = pattern_key_id(true,active_comm_pattern_keys.size()-1,active_patterns.size()-1,false);
      }
      if (should_schedule(comm_pattern_map[p_id_1])==0){
        comm_time = get_estimate(comm_pattern_map[p_id_1],comm_pattern_param,tracker.nbytes);
      }
      if (should_schedule_global(comm_pattern_map[p_id_1])==0){
        autotuning_special_bool = true;
      } else{
      }
      update(comm_pattern_map[p_id_1],comm_pattern_param,comm_time,tracker.nbytes);
  }

  // Update measurements that define the critical path for each metric.
  //critical_path_costs[num_critical_path_measures-4] += comm_time;		// update critical path communication time (for what this process has seen thus far)
  //critical_path_costs[num_critical_path_measures-3] += tracker.synch_time;	// update critical path synchronization time
  critical_path_costs[num_critical_path_measures-1] += comm_time;		// update critical path runtime

  //TODO: For autotuning_mode>0 and a skipped schedule, we want to add the barrier time as well.
  if (autotuning_mode>0 && autotuning_special_bool){
    //critical_path_costs[num_critical_path_measures-1] += tracker.barrier_time;
  }

  //volume_costs[num_volume_measures-5] += tracker.barrier_time;			// update local barrier/idle time
  //volume_costs[num_volume_measures-4] += comm_time;				// update local communication time (not volume until after the completion of the program)
  //volume_costs[num_volume_measures-3] += tracker.synch_time;			// update local synchronization time
  volume_costs[num_volume_measures-1] += (tracker.barrier_time+comm_time);	// update local runtime with idle time and comm time

  // Note that this block of code below is left in solely for blocking communication to avoid over-counting the idle time
  //   (which does not get subtracted by the min idle time any one process incurs due to efficiency complications with matching nonblocking+blocking p2p communications).
  //   Its handled correctly for blocking collectives.
  // If per-process execution-time gets larger than execution-time along the execution-time critical path, subtract out the difference from idle time.
  //volume_costs[num_volume_measures-5] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);

  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
  //volume_costs[num_volume_measures-4] = volume_costs[num_volume_measures-4] > critical_path_costs[num_critical_path_measures-4]
  //                                        ? critical_path_costs[num_critical_path_measures-4] : volume_costs[num_volume_measures-4];
  //volume_costs[num_volume_measures-3] = volume_costs[num_volume_measures-3] > critical_path_costs[num_critical_path_measures-3]
  //                                        ? critical_path_costs[num_critical_path_measures-3] : volume_costs[num_volume_measures-3];
  //volume_costs[num_volume_measures-2] = volume_costs[num_volume_measures-2] > critical_path_costs[num_critical_path_measures-2]
  //                                        ? critical_path_costs[num_critical_path_measures-2] : volume_costs[num_volume_measures-2];
  volume_costs[num_volume_measures-1] = volume_costs[num_volume_measures-1] > critical_path_costs[num_critical_path_measures-1]
                                          ? critical_path_costs[num_critical_path_measures-1] : volume_costs[num_volume_measures-1];

  comm_intercept_overhead_stage3 += MPI_Wtime() - overhead_start_time;
  overhead_start_time = MPI_Wtime();

  // Propogate critical paths for all processes in communicator based on what each process has seen up until now (not including this communication)
  if ((autotuning_mode==0) || (autotuning_special_bool==false) || (tracker.comm==MPI_COMM_WORLD && tracker.tag==0)){
    if ((rank == tracker.partner1) && (rank == tracker.partner2)) { ; }
    else{
      propagate(tracker);
      if (autotuning_mode>2 && autotuning_propagate>0){
        propagate_patterns(tracker,rank);
        // check for world communication, in which case we can flush the steady-state kernels out of the active buffers for more efficient propagation
        if ((tracker.comm == MPI_COMM_WORLD) && (tracker.partner1 == -1)){
          flush_patterns(tracker);
        }
      }
    }
    if (true_eager_p2p){
      void* temp_buf; int temp_size;
      // Forces buffered messages to send. Ideally we should wait till the next invocation of 'path::initiate(blocking&,...)' to call this,
      //   but to be safe and avoid stalls caused by MPI implementation not sending until this routine is called, we call it here.
      MPI_Buffer_detach(&temp_buf,&temp_size);
    }
  }

  comm_intercept_overhead_stage4 += MPI_Wtime() - overhead_start_time;
  // Prepare to leave interception and re-enter user code by restarting computation timers.
  tracker.start_time = MPI_Wtime();
  computation_timer = tracker.start_time;
}

// Called by both nonblocking p2p and nonblocking collectives
bool path::initiate_comm(nonblocking& tracker, volatile double curtime, volatile double itime, int64_t nelem,
                            MPI_Datatype t, MPI_Comm comm, MPI_Request* request, bool is_sender, int partner){
  return true;
}

void path::complete_comm(nonblocking& tracker, MPI_Request* request, double comp_time, double comm_time){}

void path::complete_comm(double curtime, MPI_Request* request, MPI_Status* status){}

void path::complete_comm(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){}

void path::complete_comm(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[],
                        MPI_Status array_of_statuses[]){}

void path::complete_comm(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){}

void path::flush_patterns(blocking& tracker){
  // Iterate over all computation and communication kernel pattern and
  //   flush steady-state patterns currently residing in active buffers into steady-state buffers to avoid propagation cost,
  //     as these patterns are no longer being scheduled, and thus their arithmetic mean is fixed for the rest of the program.
  // As I iterate over the entries, I will mark "is_updated=true" for those that were flushed.
  // A 3rd and final loop over active_patterns will collapse the array and update the key and value indices in the corresponding map values
  std::vector<comm_pattern_key> active_comm_pattern_keys_mirror;
  std::vector<comp_pattern_key> active_comp_pattern_keys_mirror;
  std::vector<pattern> active_patterns_mirror;
  for (auto it : comm_pattern_map){
    // We only care about those patterns that are stil active. If its not, the profile data will belong to the steady state buffers
    if (it.second.is_active){
      // If the corresponding kernel reached steady state in the last phase, flush it.
      if (active_patterns[it.second.val_index].steady_state == 1){// Reason for not using "global_steady_state" is because comp patterns do not use them.
        // Flush to steady_state_patterns
        steady_state_patterns.push_back(active_patterns[it.second.val_index]);
        steady_state_comm_pattern_keys.push_back(active_comm_pattern_keys[it.second.key_index]);
        // Update active_patterns and active_comm_pattern_keys to remove that "hole" in each array
        // Update the val/key indices to reflect the change in buffers
        comm_pattern_map[it.first].val_index = steady_state_patterns.size()-1;
        comm_pattern_map[it.first].key_index = steady_state_comm_pattern_keys.size()-1;
        comm_pattern_map[it.first].is_active = false;	// final update to make sure its known that this pattern resides in the steady-state buffers
        steady_state_patterns[steady_state_patterns.size()-1].global_steady_state=1;// prevents any more schedules in subsequent phases
        steady_state_comm_pattern_keys[steady_state_comm_pattern_keys.size()-1].pattern_index = steady_state_patterns.size()-1;
      }
      else{
        // Update active_patterns and active_comm_pattern_keys to remove that "hole" in each array
        active_comm_pattern_keys_mirror.push_back(active_comm_pattern_keys[it.second.key_index]);
        active_patterns_mirror.push_back(active_patterns[it.second.val_index]);
        comm_pattern_map[it.first].val_index = active_patterns_mirror.size()-1;
        comm_pattern_map[it.first].key_index = active_comm_pattern_keys_mirror.size()-1;
        active_comm_pattern_keys_mirror[active_comm_pattern_keys_mirror.size()-1].pattern_index = active_patterns_mirror.size()-1;
      }
    }
  }
  for (auto it : comp_pattern_map){
    // We only care about those patterns that are stil active
    if (it.second.is_active){
      if (active_patterns[it.second.val_index].steady_state == 1){// Reason for not using "global_steady_state" is because comp patterns do not use them.
        // Flush to steady_state_patterns
        steady_state_patterns.push_back(active_patterns[it.second.val_index]);
        steady_state_comp_pattern_keys.push_back(active_comp_pattern_keys[it.second.key_index]);
        // Update active_patterns and active_comp_pattern_keys to remove that "hole" in each array
        // Update the val/key indices to reflect the change in buffers
        comp_pattern_map[it.first].val_index = steady_state_patterns.size()-1;
        comp_pattern_map[it.first].key_index = steady_state_comp_pattern_keys.size()-1;
        comp_pattern_map[it.first].is_active = false;	// final update to make sure its known that this pattern resides in the steady-state buffers
        steady_state_patterns[steady_state_patterns.size()-1].global_steady_state=1;// prevents any more schedules in subsequent phases
        steady_state_comp_pattern_keys[steady_state_comp_pattern_keys.size()-1].pattern_index = steady_state_patterns.size()-1;
      }
      else{
        // Update active_patterns and active_comp_pattern_keys to remove that "hole" in each array
        active_comp_pattern_keys_mirror.push_back(active_comp_pattern_keys[it.second.key_index]);
        active_patterns_mirror.push_back(active_patterns[it.second.val_index]);
        comp_pattern_map[it.first].val_index = active_patterns_mirror.size()-1;
        comp_pattern_map[it.first].key_index = active_comp_pattern_keys_mirror.size()-1;
        active_comp_pattern_keys_mirror[active_comp_pattern_keys_mirror.size()-1].pattern_index = active_patterns_mirror.size()-1;
      }
    }
  }

  active_comm_pattern_keys = active_comm_pattern_keys_mirror;
  active_comp_pattern_keys = active_comp_pattern_keys_mirror;
  active_patterns = active_patterns_mirror;
}


void path::propagate_patterns(blocking& tracker, int rank){
  // Use info_receiver[last].second when deciding who to issue 3 broadcasts from
  // First need to broadcast the size of each of the 3 broadcasts so that the receiving buffers can prepare the size of their receiving buffers
  // Only the active kernels need propagating. Steady-state are treated differently depending on the communicator.

  bool true_eager_p2p = ((eager_p2p == 1) && (tracker.tag!=13) && (tracker.tag!=14));
  int size_array[3] = {0,0,0};
  if (rank == info_receiver[num_critical_path_measures-1].second){
    size_array[0] = active_comm_pattern_keys.size();
    size_array[1] = active_comp_pattern_keys.size();
    size_array[2] = active_patterns.size();
  }
  if (tracker.partner1 == -1){
    PMPI_Bcast(&size_array[0],3,MPI_INT,info_receiver[num_critical_path_measures-1].second,tracker.comm);
    if (rank != info_receiver[num_critical_path_measures-1].second){
        active_comm_pattern_keys.resize(size_array[0]);
        active_comp_pattern_keys.resize(size_array[1]);
        active_patterns.resize(size_array[2]);
    }
    PMPI_Bcast(&active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,info_receiver[num_critical_path_measures-1].second,tracker.comm);
    PMPI_Bcast(&active_comp_pattern_keys[0],size_array[1],comp_pattern_key_type,info_receiver[num_critical_path_measures-1].second,tracker.comm);
    PMPI_Bcast(&active_patterns[0],size_array[2],pattern_type,info_receiver[num_critical_path_measures-1].second,tracker.comm);
    // receivers update their maps and data structures. Yes, this is costly, but I guess this is what has to happen.
    // basically everything is set via the broadcasts themselves, except for the two maps.
  }
  else{
    // We leverage the fact that we know the path-defining process rank
    if (true_eager_p2p){
      if (tracker.is_sender){
        PMPI_Bsend(&size_array[0],3,MPI_INT,tracker.partner1,internal_tag2,tracker.comm);
        PMPI_Bsend(&active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,tracker.partner1,internal_tag2,tracker.comm);
        PMPI_Bsend(&active_comp_pattern_keys[0],size_array[1],comp_pattern_key_type,tracker.partner1,internal_tag2,tracker.comm);
        PMPI_Bsend(&active_patterns[0],size_array[2],pattern_type,tracker.partner1,internal_tag2,tracker.comm);
      } else{
        PMPI_Recv(&size_array[0],3,MPI_INT,tracker.partner1,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
        active_comm_pattern_keys.resize(size_array[0]);
        active_comp_pattern_keys.resize(size_array[1]);
        active_patterns.resize(size_array[2]);
        PMPI_Recv(&active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,tracker.partner1,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
        PMPI_Recv(&active_comp_pattern_keys[0],size_array[1],comp_pattern_key_type,tracker.partner1,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
        PMPI_Recv(&active_patterns[0],size_array[2],pattern_type,tracker.partner1,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      }
    }
    else{
      if (rank == info_receiver[num_critical_path_measures-1].second){
        PMPI_Send(&size_array[0],3,MPI_INT,tracker.partner1,internal_tag2,tracker.comm);
        PMPI_Send(&active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,tracker.partner1,internal_tag2,tracker.comm);
        PMPI_Send(&active_comp_pattern_keys[0],size_array[1],comp_pattern_key_type,tracker.partner1,internal_tag2,tracker.comm);
        PMPI_Send(&active_patterns[0],size_array[2],pattern_type,tracker.partner1,internal_tag2,tracker.comm);
      } else{
        PMPI_Recv(&size_array[0],3,MPI_INT,tracker.partner1,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
        active_comm_pattern_keys.resize(size_array[0]);
        active_comp_pattern_keys.resize(size_array[1]);
        active_patterns.resize(size_array[2]);
        PMPI_Recv(&active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,tracker.partner1,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
        PMPI_Recv(&active_comp_pattern_keys[0],size_array[1],comp_pattern_key_type,tracker.partner1,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
        PMPI_Recv(&active_patterns[0],size_array[2],pattern_type,tracker.partner1,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      }
      //TODO: Note that I may be screwing up the case in which tracker.partner1 != tracker.partner2
    }
  }

  // Lets have all processes update, even the root, so that they leave this routine (and subsequently leave the interception) at approximately the same time.
  for (auto i=0; i<active_comm_pattern_keys.size(); i++){
    comm_pattern_key id(active_comm_pattern_keys[i].pattern_index,active_comm_pattern_keys[i].tag,active_comm_pattern_keys[i].comm_size,
                              active_comm_pattern_keys[i].comm_color,active_comm_pattern_keys[i].msg_size,active_comm_pattern_keys[i].partner_offset); 
    if (comm_pattern_map.find(id) == comm_pattern_map.end()){
      comm_pattern_map[id] = pattern_key_id(true, i, active_comm_pattern_keys[i].pattern_index, true);
    } else{
      comm_pattern_map[id].is_active = true;	// always assumed true
      comm_pattern_map[id].key_index = i;
      comm_pattern_map[id].val_index = active_comm_pattern_keys[i].pattern_index;
      comm_pattern_map[id].is_updated = true;
    }
  }

  // Delete those keys that no longer lie along the critical path
  // TODO: I may be able to just skip this loop if I am not deleting anything.
  //   NO!!!!! I cannot just delete it. The problem is that these indices are no longer valid, so when I go to update, I will be updating the wrong kernel.
  for (auto it = comm_pattern_map.begin(); it != comm_pattern_map.end();){
    if (!it->second.is_active){ it++; continue;  }
    if (!it->second.is_updated){
      it = comm_pattern_map.erase(it);
    }
    else{
      it->second.is_updated=false;	// to prepare for next propagation
      it++;
    }
  }
  for (auto i=0; i<active_comp_pattern_keys.size(); i++){
    comp_pattern_key id(active_comp_pattern_keys[i].pattern_index,active_comp_pattern_keys[i].tag,active_comp_pattern_keys[i].flops,
                               active_comp_pattern_keys[i].param1,active_comp_pattern_keys[i].param2,active_comp_pattern_keys[i].param3,
                               active_comp_pattern_keys[i].param4,active_comp_pattern_keys[i].param5); 
    if (comp_pattern_map.find(id) == comp_pattern_map.end()){
      comp_pattern_map[id] = pattern_key_id(true, i, active_comp_pattern_keys[i].pattern_index,true);
    } else{
      comp_pattern_map[id].is_active = true;	// always assumed true
      comp_pattern_map[id].key_index = i;
      comp_pattern_map[id].val_index = active_comp_pattern_keys[i].pattern_index;
      comp_pattern_map[id].is_updated = true;
    }
  }
  // Delete those keys that no longer lie along the critical path
  for (auto it = comp_pattern_map.begin(); it != comp_pattern_map.end();){
    if (!it->second.is_active){ it++; continue;  }
    if (!it->second.is_updated){
      it = comp_pattern_map.erase(it);
    }
    else{
      it->second.is_updated=false;	// to prepare for next propagation
      it++;
    }
  }
}

void path::propagate(blocking& tracker){
  assert(tracker.comm != 0);
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  if ((rank == tracker.partner1) && (rank == tracker.partner2)) { return; } 
  bool true_eager_p2p = ((eager_p2p == 1) && (tracker.tag!=13) && (tracker.tag!=14));
  if (autotuning_mode>2){// Autotuning using critical path analysis requires knowledge of which rank determined the cp
    //TODO: Idea for 2-stage reduction: move this out of the mode>=2 if statement, and then after this, scan the critical_path_costs and zero out what is not defining a critical path and then post a MPI_Allreduce (via multi-root hack)
    for (int i=0; i<num_critical_path_measures; i++){
      info_sender[i].first = critical_path_costs[i];
      info_sender[i].second = rank;
    }
    if (tracker.partner1 == -1){
      PMPI_Allreduce(&info_sender[0].first, &info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, MPI_MAXLOC, tracker.comm);
      for (int i=0; i<num_critical_path_measures; i++){
        critical_path_costs[i] = info_receiver[i].first;
      }
    }
    else{
      if (!true_eager_p2p){
        PMPI_Sendrecv(&info_sender[0].first, num_critical_path_measures, MPI_DOUBLE_INT, tracker.partner1, internal_tag,
                      &info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, tracker.partner2, internal_tag, tracker.comm, MPI_STATUS_IGNORE);
        for (int i=0; i<num_critical_path_measures; i++){
          if (info_sender[i].first>info_receiver[i].first){info_receiver[i].second = rank;}
          else if (info_sender[i].first==info_receiver[i].first){ info_receiver[i].second = std::min(rank,tracker.partner1); }
          info_receiver[i].first = std::max(info_sender[i].first, info_receiver[i].first);
          critical_path_costs[i] = info_receiver[i].first;
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
            critical_path_costs[i] = info_receiver[i].first;
          }
        }
      }
      else{
        if (tracker.is_sender){
          PMPI_Bsend(&info_sender[0].first, num_critical_path_measures, MPI_DOUBLE_INT, tracker.partner1, internal_tag, tracker.comm);
        } else{
          PMPI_Recv(&info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, tracker.partner1, internal_tag, tracker.comm, MPI_STATUS_IGNORE);
          //TODO: Is there a bug here? Why is nothing being updated?
        }
      }
    }
  }
  else{
 
   // Exchange the tracked routine critical path data
    if (tracker.partner1 == -1){
      PMPI_Allreduce(MPI_IN_PLACE, &critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, MPI_MAX, tracker.comm);
    }
    else{
      // Note that a blocking sendrecv allows exchanges even when the other party issued a request via nonblocking communication, as the process with the nonblocking request posts both sends and receives.
      if (true_eager_p2p){
        if (tracker.is_sender){
          PMPI_Bsend(&critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, tracker.partner1, internal_tag2, tracker.comm);
        } else{
          PMPI_Recv(&new_cs[0], critical_path_costs.size(), MPI_DOUBLE, tracker.partner1, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
          update_critical_path(&new_cs[0],&critical_path_costs[0],critical_path_costs_size);
        }
      }
      else{
        PMPI_Sendrecv(&critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, tracker.partner1, internal_tag2, &new_cs[0], critical_path_costs.size(),
                      MPI_DOUBLE, tracker.partner2, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
        update_critical_path(&new_cs[0],&critical_path_costs[0],critical_path_costs_size);
        if (tracker.partner2 != tracker.partner1){
          // This if-statement will never be breached if 'true_eager_p2p'=true anyways.
          PMPI_Sendrecv(&critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, tracker.partner2, internal_tag2, &new_cs[0], critical_path_costs.size(), MPI_DOUBLE, tracker.partner1, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
          update_critical_path(&new_cs[0],&critical_path_costs[0],critical_path_costs_size);
        }
      }
    }
  }
}

void path::propagate(nonblocking& tracker){}

}
}
}
