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
  assert(new_comm_size>0);
  //int old_comm_rank; MPI_Comm_rank(oldcomm,&old_comm_rank);
  int world_comm_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_comm_rank);
  std::vector<int> gathered_ranks(new_comm_size,0);
  PMPI_Allgather(&world_comm_rank,1,MPI_INT,&gathered_ranks[0],1,MPI_INT,newcomm);
  // Now we detect the "color" (really the stride) via iteration
  // Step 1: subtract out the offset from 0 : assuming that the key arg to comm_split didn't re-shuffle
  //         I can try to use std::min_element vs. writing my own manual loop
  int offset = gathered_ranks[0];
  for (auto i=1; i<gathered_ranks.size(); i++) { offset = std::min(offset,gathered_ranks[i]); }
  for (auto i=0; i<gathered_ranks.size(); i++) { gathered_ranks[i] -= offset; }
  std::sort(gathered_ranks.begin(),gathered_ranks.end());
  comm_channel_node* channel = new comm_channel_node();
  channel->offset = offset;
  if (new_comm_size<=1){
    channel->id.push_back(std::make_pair(new_comm_size,world_comm_rank));
  }
  else{
    int stride = gathered_ranks[1]-gathered_ranks[0];
    int count = 0;
    int jump=1;
    int extra=0;
    int i=0;
    while (i < gathered_ranks.size()-1){
      if ((gathered_ranks[i+jump]-gathered_ranks[i]) != stride){
        channel->id.push_back(std::make_pair(count+extra+1,stride));
        stride = gathered_ranks[i+1]-gathered_ranks[0];
        i += jump;
        if (channel->id.size()==1){
          jump=count+extra+1;
        }
        else{
          jump = (count+extra+1)*channel->id[channel->id.size()-2].first;
        }
        extra=1;
        count = 0;
      } else{
        count++;
        i += jump;
      }
    }
    if (count != 0){
      channel->id.push_back(std::make_pair(count+extra+1,stride));
    }
  }
  spf.insert_node(channel);// This call will just fill in SPT via channel's parent/children members, and the members of related channels
  comm_channel_map[newcomm] = channel;

  int color = new_comm_size>1 ? gathered_ranks[1]-gathered_ranks[0] : 0;
  communicator_map[newcomm] = std::make_pair(new_comm_size,color);
  computation_timer = MPI_Wtime();
}

bool path::initiate_comp(size_t id, volatile double curtime, double flop_count, int param1, int param2, int param3, int param4, int param5){
  // accumulate computation time
  double save_comp_time = curtime - computation_timer;
  critical_path_costs[num_critical_path_measures-1] += save_comp_time;	// update critical path execution time
  critical_path_costs[num_critical_path_measures-2] += save_comp_time;	// update critical path computation time
  volume_costs[num_volume_measures-1]        += save_comp_time;		// update local execution time
  volume_costs[num_volume_measures-2]        += save_comp_time;		// update local computation time

  // Special exit if no kernels are to be scheduled -- the whole point is to get a reading on the total time it takes, which is
  //   to be attained with timers outside of critter..
  if (schedule_kernels==0){ return false; }
  volatile double overhead_start_time = MPI_Wtime();

  bool schedule_decision = true;
  comp_pattern_key key(-1,id,flop_count,param1,param2,param3,param4,param5);// '-1' argument is arbitrary, does not influence overloaded operators
  if (!(comp_pattern_map.find(key) == comp_pattern_map.end())){
    schedule_decision = should_schedule(comp_pattern_map[key])==1;
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

  comp_pattern_key key(active_patterns.size(),id,flop_count,param1,param2,param3,param4,param5);// 'active_patterns.size()' argument is arbitrary, does not influence overloaded operators
  if (comp_pattern_map.find(key) == comp_pattern_map.end()){
    active_comp_pattern_keys.emplace_back(key);
    active_patterns.emplace_back(pattern());
    comp_pattern_map[key] = pattern_key_id(true,active_comp_pattern_keys.size()-1,active_patterns.size()-1,false);
  }
  update_kernel_stats(comp_pattern_map[key],comp_analysis_param,comp_time,flop_count);
  // Note: 'get_estimate' must be called before setting the updated kernel state.
  if (should_schedule(comp_pattern_map[key]) == 0){
    comp_time = get_estimate(comp_pattern_map[key],comp_analysis_param,flop_count);
  }
  // Both non-optimized and optimized variants can update the local kernel state, but not the global kernel state
  // This gives unoptimized variant more license that that for a communication pattern (which cannot even update
  //   local steady state within a phase, due to potential of deadlock if done so).
  bool _is_steady = steady_test(key,comp_pattern_map[key],comp_analysis_param);
  set_kernel_state(comp_pattern_map[key],!_is_steady);

  critical_path_costs[num_critical_path_measures-1] += comp_time;
  critical_path_costs[num_critical_path_measures-2] += comp_time;
  critical_path_costs[num_critical_path_measures-3] += comp_time;
  volume_costs[num_volume_measures-1] += comp_time;
  volume_costs[num_volume_measures-2] += comp_time;
  volume_costs[num_volume_measures-3] += comp_time;

  comp_intercept_overhead += MPI_Wtime() - overhead_start_time;
  computation_timer = MPI_Wtime();
}

bool path::initiate_comm(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                         bool is_sender, int partner1, int partner2){
  // Save and accumulate the computation time between last communication routine as both execution-time and computation time
  //   into both the execution-time critical path data structures and the per-process data structures.
  tracker.comp_time = curtime - computation_timer;
  critical_path_costs[num_critical_path_measures-1] += tracker.comp_time;	// update critical path execution time
  critical_path_costs[num_critical_path_measures-2] += tracker.comp_time;	// update critical path computation time
  volume_costs[num_volume_measures-1]        += tracker.comp_time;		// update local runtime
  volume_costs[num_volume_measures-2]        += tracker.comp_time;		// update local runtime
  // Special exit if no kernels are to be scheduled -- the whole point is to get a reading on the total time it takes, which is
  //   to be attained with timers outside of critter..
  if (schedule_kernels==0){ return false; }

  // At this point, 'critical_path_costs' tells me the process's time up until now. A barrier won't suffice.
  int rank; MPI_Comm_rank(comm, &rank);
  volatile double overhead_start_time = MPI_Wtime();

  assert(comm_channel_map.find(comm) != comm_channel_map.end());// Any sub-communicator must have been registered via comm_split
  comm_channel_node* tree_node;
  if (tracker.partner1 != -1){// p2p
    auto world_rank = spf.translate_rank(comm,tracker.partner1);
    if (p2p_channel_map.find(world_rank) == p2p_channel_map.end()){
      comm_channel_node* node = new comm_channel_node();
      node->offset = world_rank;
      node->id.push_back(std::make_pair(1,0));
      spf.insert_node(node);
      p2p_channel_map[world_rank] = node;
    }
    tree_node = p2p_channel_map[world_rank];
  } else{
    tree_node = comm_channel_map[comm];
  }

  // We consider usage of Sendrecv variants to forfeit usage of eager internal communication.
  // Note that the reason we can't force user Bsends to be 'true_eager_p2p' is because the corresponding Receives would be expecting internal communications
  bool true_eager_p2p = ((eager_p2p == 1) && (tracker.tag!=13) && (tracker.tag!=14));
  if (true_eager_p2p){ MPI_Buffer_attach(&eager_pad[0],eager_pad.size()); }

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
  tracker.barrier_time = 0.;// might get updated below

  // Non-optimized variant will always post barriers, although of course, just as with the optimized variant, the barriers only remove idle time
  //   from corrupting communication time measurements. The process that enters barrier last is not necessarily the critical path root. The
  //     critical path root is decided based on a reduction using 'critical_path_costs'.
  bool post_barrier = true;
  bool schedule_decision = true;
  double reduced_info[5] = {critical_path_costs[num_critical_path_measures-4],critical_path_costs[num_critical_path_measures-3],
                            critical_path_costs[num_critical_path_measures-2],critical_path_costs[num_critical_path_measures-1],0};
  double reduced_info_foreign[5] = {0,0,0,0,0};
  comm_pattern_key key(rank,-1,tracker.tag,communicator_map[tracker.comm].first,communicator_map[tracker.comm].second,tracker.nbytes,tracker.partner1);
  if (!(comm_pattern_map.find(key) == comm_pattern_map.end())){
    if (is_optimized){ post_barrier = should_schedule_global(comm_pattern_map[key])==1; }
    schedule_decision = should_schedule(comm_pattern_map[key])==1;
    reduced_info[4] = (double)should_schedule(comm_pattern_map[key]);
  }

  // Both unoptimized and optimized variants will reduce an execution time and a schedule_decision.
  // Note that formally, the unoptimized variant does not need to reduce a schedule_decision, but I allow it becase the
  //   procedure is synchronization bound anyway.
  // The unoptimized variant will always perform this reduction of execution time, regardless of its schedule_decision.
  // The optimized variant will only perform this reduction until it detects steady state.
  // For the optimized variant, if kernel is not in global steady state, in effort to accelerate convergence
  //   to both steady state and global steady state (i.e. early detection of steady state convergence at expense of extra communication)
  //   This early detection optimization is not foolproof, hence why the non-optimized variant does not consider it.
  if (post_barrier==true){
    assert(partner1 != MPI_ANY_SOURCE); if ((tracker.tag == 13) || (tracker.tag == 14)){ assert(partner2 != MPI_ANY_SOURCE); }
    if (partner1 == -1){
      PMPI_Allreduce(MPI_IN_PLACE, &reduced_info[0], 5, MPI_DOUBLE, MPI_MAX, tracker.comm);
      schedule_decision = reduced_info[4] > 0;
    }
    else{
      if ((true_eager_p2p) && (rank != tracker.partner1)){
        assert(0);// Not implemented yet. Does eager protocol even make sense, given these methods?
/*
        if (tracker.is_sender){
          PMPI_Bsend(&schedule_decision_int, 2, MPI_DOUBLE, tracker.partner1, internal_tag2, tracker.comm);
        }
        else{
          PMPI_Recv(&schedule_decision_foreign_int, 1, MPI_DOUBLE, tracker.partner1, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
          schedule_decision = (schedule_decision_foreign_int > 0 ? true : false);
        }
*/
      }
      else if (!true_eager_p2p){
        PMPI_Sendrecv(&reduced_info[0], 5, MPI_DOUBLE, tracker.partner1, internal_tag2, &reduced_info_foreign[0], 5,
                      MPI_DOUBLE, tracker.partner2, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
        for (auto i=0; i<4; i++){ reduced_info[i] = std::max(reduced_info[i],reduced_info_foreign[i]); }
        schedule_decision = (reduced_info[4] > 0) || (reduced_info_foreign[4]>0);
        if (tracker.partner2 != tracker.partner1){
          // This if-statement will never be breached if 'true_eager_p2p'=true anyways.
          PMPI_Sendrecv(&reduced_info[0], 5, MPI_DOUBLE, tracker.partner2, internal_tag2, &reduced_info_foreign[0], 5,
                        MPI_DOUBLE, tracker.partner1, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
          for (auto i=0; i<4; i++){ reduced_info[i] = std::max(reduced_info[i],reduced_info_foreign[i]); }
          schedule_decision = schedule_decision || (reduced_info[4] > 0) || (reduced_info_foreign[4]>0);
        }
      }
    }
    tracker.barrier_time = reduced_info[3] - critical_path_costs[num_critical_path_measures-1];

    // If local steady_state has been reached, and we find out the other processes have reached the same, then we can set global_steady_state=1
    //   and never have to use another internal collective to check again.
    if (is_optimized){
      if (!(comm_pattern_map.find(key) == comm_pattern_map.end())){
        assert(should_schedule_global(comm_pattern_map[key])==1);
        set_kernel_state(comm_pattern_map[key],schedule_decision);
        if (analysis_mode == 3){
          set_kernel_state_global(comm_pattern_map[key],true);// prevents critical path analysis from changing global steady_state, even in optimized variant
        } else{
          set_kernel_state_global(comm_pattern_map[key],schedule_decision);
          if (!schedule_decision) { flush_pattern(key); }
        }
      }
    }
    // If kernel is about to be scheduled, post one more barrier for safety if collective,
    //   because the AllReduce posted above may allow ranks to leave early, thus corrupting the sample measurement.
    if (schedule_decision && tracker.partner1 == -1){ PMPI_Barrier(tracker.comm); }
  }

  comm_intercept_overhead_stage1 += MPI_Wtime() - overhead_start_time;
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

  bool should_propagate = true;
  bool should_update_idle = true;
  comm_pattern_key key(rank,active_patterns.size(),tracker.tag,communicator_map[tracker.comm].first,communicator_map[tracker.comm].second,tracker.nbytes,tracker.partner1);
  assert(communicator_map.find(tracker.comm) != communicator_map.end());
  if (comm_pattern_map.find(key) == comm_pattern_map.end()){
    active_comm_pattern_keys.emplace_back(key);
    active_patterns.emplace_back(pattern());
    comm_pattern_map[key] = pattern_key_id(true,active_comm_pattern_keys.size()-1,active_patterns.size()-1,false);
  }
  if (comm_pattern_pair_map.find(std::make_pair(previous_comm_key,key)) == comm_pattern_pair_map.end()){
    comm_pattern_pair_map[std::make_pair(previous_comm_key,key)] = idle_pattern();
  }
  if (should_schedule(comm_pattern_map[key])==0){
    comm_time = get_estimate(comm_pattern_map[key],comm_analysis_param,tracker.nbytes);
    if (is_optimized==1){
/*
      if (tracker.barrier_time == -1) { should_update_idle = false;
                                        tracker.barrier_time = comm_pattern_pair_map[std::make_pair(previous_comm_key,key)].M1;
                                      }
*/
      should_propagate = false;
    }
  }
  update_kernel_stats(comm_pattern_map[key],comm_analysis_param,comm_time,tracker.nbytes);
  //update_kernel_stats(comm_pattern_pair_map[std::make_pair(previous_comm_key,key)],!should_update_idle,tracker.barrier_time);
  if (is_optimized==1){
    bool is_steady = steady_test(key,comm_pattern_map[key],comm_analysis_param);
    set_kernel_state(comm_pattern_map[key],!is_steady);
  }

  critical_path_costs[num_critical_path_measures-1] += (comm_time + tracker.barrier_time);
  critical_path_costs[num_critical_path_measures-4] += (comm_time + tracker.barrier_time);
  volume_costs[num_volume_measures-1] += (comm_time + tracker.barrier_time);
  volume_costs[num_volume_measures-4] += (comm_time + tracker.barrier_time);

  // Propogate critical paths for all processes in communicator based on what each process has seen up until now (not including this communication)
  if (should_propagate){
    bool is_world_communication = (tracker.comm == MPI_COMM_WORLD) && (tracker.partner1 == -1);
    if ((rank == tracker.partner1) && (rank == tracker.partner2)) { ; }
    else{
      propagate(tracker);
      if (analysis_mode==1){
        // Note: in optimized version, exchange_patterns and flush_patterns would never need to be called
        if (is_optimized==0){
          if (is_world_communication && !is_key_skipable(key)){//Note: for practical reasons, we force more constraints on when profile exchanges take place
            exchange_patterns_per_process(tracker);
            flush_patterns();
          }
        }
      }
      else if (analysis_mode==2){
        // Note: in both unoptimized and optimized versions, exchange_patterns and flush_patterns are useful for merging many samples than possib;e
        //       with per-process analysis. The goal here, unlike with per-process, is not simply to detect global steady state.
        if (is_world_communication && !is_key_skipable(key)){
          exchange_patterns_volumetric(tracker);
          //flush_patterns();
        }
      }
      else if (analysis_mode==3){
        propagate_patterns(tracker,key,rank);
        // check for world communication, in which case we can flush the steady-state kernels out of the active buffers for more efficient propagation
        if (is_world_communication && !is_key_skipable(key)){
          flush_patterns();
        }
      }
    }
  }
  if (true_eager_p2p){
    void* temp_buf; int temp_size;
    // Forces buffered messages to send. Ideally we should wait till the next invocation of 'path::initiate(blocking&,...)' to call this,
    //   but to be safe and avoid stalls caused by MPI implementation not sending until this routine is called, we call it here.
    MPI_Buffer_detach(&temp_buf,&temp_size);
  }

  comm_intercept_overhead_stage2 += MPI_Wtime() - overhead_start_time;
  // Prepare to leave interception and re-enter user code by restarting computation timers.
  previous_comm_key = key;
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

void path::flush_pattern(comm_pattern_key key){
  return;
  // 'flush_pattern' is used for the sole purpose of keeping the invariant that an inactive kernel's statistical data
  //    belongs in the steady_state buffers. This procedure need only occur in 'flush_patterns' for critical path analysis,
  //      but as in per-process/volumetric analysis the global state can change in 'initiate_comm', this routine serves as a way to
  //        achieve that.
  assert(comm_pattern_map[key].is_active);
  assert(analysis_mode>=1 && analysis_mode<=2);
  assert(is_optimized);
  assert(active_patterns[comm_pattern_map[key].val_index].steady_state==1);
  assert(active_patterns[comm_pattern_map[key].val_index].global_steady_state==1);
  if (active_comp_pattern_keys.size()>0){
    //TODO: note: this assert is not valid, because after a few flushes, the last keys might point to intermediate patterns.
    assert((active_comm_pattern_keys[active_comm_pattern_keys.size()-1].pattern_index == active_patterns.size()-1) ||
           (active_comp_pattern_keys[active_comp_pattern_keys.size()-1].pattern_index == active_patterns.size()-1));
  }

  // Flush to steady_state_patterns
  steady_state_patterns.push_back(active_patterns[key.pattern_index]);
  steady_state_comm_pattern_keys.push_back(key);

  // First check for special corner case that the key we want to remove is the last key in the active key array. The reason this
  //   is a case to check because it prevents the need to swap with a different active comm key. If this case does not pass,
  //   then we know we will have to swap with a different active comm key (the last one in the active comm key array)
  if (comm_pattern_map[key].key_index == (active_comm_pattern_keys.size()-1)){
    active_comm_pattern_keys.pop_back();
    if (key.pattern_index == (active_patterns.size()-1)){ active_patterns.pop_back(); }
    else{
      // We must swap the key's pattern with the last active pattern, which must be a comp, because the comm key we are operating on is the last
      // We can leverage the fact that in this case, there must be at least one active comp pattern, because it must be the one stored in the last entry
      //   of active_patterns. We should include some asserts to check for that
      assert(active_comp_pattern_keys.size()>0);
      assert(active_comp_pattern_keys[active_comp_pattern_keys.size()-1].pattern_index == active_patterns.size()-1);
      // Update active_patterns and active_comm_pattern_keys to remove that "hole" in each array
      auto lra_active_pattern = active_patterns[active_patterns.size()-1];
      auto lra_comp_key = active_comp_pattern_keys[active_comp_pattern_keys.size()-1];
      active_patterns[key.pattern_index] = lra_active_pattern;
      active_patterns.pop_back();
      // last comp pattern gets swapped into place of newly steady comm pattern (above), and pattern_index must be updated
      active_comp_pattern_keys[active_comp_pattern_keys.size()-1].pattern_index = key.pattern_index;
      comp_pattern_map[lra_comp_key].val_index = key.pattern_index;
    }
  }
  else{
    // Find the key of the last entry in active_patterns (assumed to be the last entry in either
    //   active_comm_pattern_keys or active_comp_pattern_keys
    auto lra_active_pattern = active_patterns[active_patterns.size()-1];
    auto lra_comm_key = active_comm_pattern_keys[active_comm_pattern_keys.size()-1];
    auto lra_comp_key = active_comp_pattern_keys[active_comp_pattern_keys.size()-1];
    int dir;
    if (active_comp_pattern_keys.size()==0){ dir = 1; }
    else if (active_comm_pattern_keys[active_comm_pattern_keys.size()-1].pattern_index <
             active_comp_pattern_keys[active_comp_pattern_keys.size()-1].pattern_index){ dir = 0; }
    else{ dir = 1; }
    // Update active_patterns and active_comm_pattern_keys to remove that "hole" in each array
    active_patterns[key.pattern_index] = lra_active_pattern;
    active_patterns.pop_back();
    if (dir==0){
      // last comp pattern gets swapped into place of newly steady comm pattern (above), and pattern_index must be updated
      active_comp_pattern_keys[active_comp_pattern_keys.size()-1].pattern_index = key.pattern_index;
      comp_pattern_map[lra_comp_key].val_index = key.pattern_index;
      // last comm key gets swapped into place of newly steady comm pattern's key, and all relevant indices must be updated
      active_comm_pattern_keys[comm_pattern_map[key].key_index] = lra_comm_key;
      active_comm_pattern_keys.pop_back();
      comm_pattern_map[lra_comm_key].key_index = comm_pattern_map[key].key_index;
    } else{
      // last comm pattern gets swapped into place of newly steady comm pattern (above), and all relevant indices must be updated
      active_comm_pattern_keys[comm_pattern_map[key].key_index] = lra_comm_key;
      active_comm_pattern_keys.pop_back();
      active_comm_pattern_keys[comm_pattern_map[key].key_index].pattern_index = key.pattern_index;
      comm_pattern_map[lra_comm_key].val_index = comm_pattern_map[key].val_index;
      comm_pattern_map[lra_comm_key].key_index = comm_pattern_map[key].key_index;
    }
  }

  // Update the val/key indices to reflect the change in buffers
  comm_pattern_map[key].val_index = steady_state_patterns.size()-1;
  comm_pattern_map[key].key_index = steady_state_comm_pattern_keys.size()-1;
  comm_pattern_map[key].is_active = false;	// final update to make sure its known that this pattern resides in the steady-state buffers
  //steady_state_patterns[steady_state_patterns.size()-1].global_steady_state=1;// prevents any more schedules in subsequent phases
  steady_state_comm_pattern_keys[steady_state_comm_pattern_keys.size()-1].pattern_index = steady_state_patterns.size()-1;
  // Note: As this routine is to be used for optimized per-process/volumetric statistical analysis, no special care is needed for p2p.
}

void path::flush_patterns(){

  // Iterate over all computation and communication kernel patterns and
  //   flush steady-state patterns currently residing in active buffers into steady-state buffers to avoid propagation cost,
  //     in subsequent exchanges or propagations, as these patterns are no longer being scheduled,
  //     and thus their arithmetic mean is fixed for the rest of the program.
  // As I iterate over the entries, I will mark "is_updated=true" for those that were flushed.
  // A 3rd and final loop over active_patterns will collapse the array and update the key and value indices in the corresponding map values
  std::vector<comm_pattern_key> active_comm_pattern_keys_mirror;
  std::vector<comp_pattern_key> active_comp_pattern_keys_mirror;
  std::vector<pattern> active_patterns_mirror;
  std::map<comm_pattern_key,pattern_key_id> p2p_map;
  for (auto it : comm_pattern_map){
    // We only care about those patterns that are stil active. If its not, the profile data will belong to the steady state buffers
    if (it.second.is_active){
      // Separate out p2p from collectives. p2p coordination is more difficult to manage
      if (it.first.tag == 16){
        auto key_copy = it.first; key_copy.tag = 17;
        if (comm_envelope_param == 0) { key_copy.partner_offset *= (-1); }
        else if (comm_envelope_param == 1) { key_copy.partner_offset = abs(key_copy.partner_offset); }
        else { key_copy.partner_offset=-1; }
        if (comm_pattern_map.find(key_copy) == comm_pattern_map.end()){
          active_patterns[it.second.val_index].steady_state = 0;
          active_patterns[it.second.val_index].global_steady_state = 0;
        }
        else{// debug check
          //assert(it.second.is_active == comm_pattern_map[key_copy].is_active);// debug
        }
        p2p_map[it.first] = it.second;
      }
      else if (it.first.tag == 17){
        bool is_steady = true;
        auto key_copy = it.first; key_copy.tag = 16;
        if (comm_envelope_param == 0) { key_copy.partner_offset *= (-1); }
        else if (comm_envelope_param == 1) { key_copy.partner_offset = abs(key_copy.partner_offset); }
        else { key_copy.partner_offset=-1; }
        if (p2p_map.find(key_copy) == p2p_map.end()){
          is_steady = false;
        }
        else{
          //assert(it.second.is_active == comm_pattern_map[key_copy].is_active);// debug
          is_steady = (active_patterns[it.second.val_index].steady_state == 1) && (active_patterns[p2p_map[key_copy].val_index].steady_state == 1);
        }
        // If the corresponding kernel reached steady state in the last phase, flush it.
        if (is_steady){
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
          // corresponding send will get flushed below
        } else{
          // Update active_patterns and active_comm_pattern_keys to remove that "hole" in each array
          active_comm_pattern_keys_mirror.push_back(active_comm_pattern_keys[it.second.key_index]);
          active_patterns_mirror.push_back(active_patterns[it.second.val_index]);
          comm_pattern_map[it.first].val_index = active_patterns_mirror.size()-1;
          comm_pattern_map[it.first].key_index = active_comm_pattern_keys_mirror.size()-1;
          active_comm_pattern_keys_mirror[active_comm_pattern_keys_mirror.size()-1].pattern_index = active_patterns_mirror.size()-1;
          // must mark corresponding send, if available, to not be in steady_state
          if (p2p_map.find(key_copy) != p2p_map.end()){
            active_patterns[p2p_map[key_copy].val_index].steady_state = 0;
            active_patterns[p2p_map[key_copy].val_index].global_steady_state = 0;
          }
        }
      }
      else{
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
  }
  for (auto it : p2p_map){
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


void path::propagate_patterns(blocking& tracker, comm_pattern_key comm_key, int rank){
  // Use info_receiver[last].second when deciding who to issue 3 broadcasts from
  // First need to broadcast the size of each of the 3 broadcasts so that the receiving buffers can prepare the size of their receiving buffers
  // Only the active kernels need propagating. Steady-state are treated differently depending on the communicator.

  auto local_pattern = active_patterns[comm_pattern_map[comm_key].val_index];

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
        //TODO: I want to keep both variants. The first is easy to rewrite, just take from above.
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
  bool key_match = ( (comm_key.tag <= 12) || (comm_key.tag >= 19) ? true : false);// Only start as false if dealing with p2p
  for (auto i=0; i<active_comm_pattern_keys.size(); i++){
    comm_pattern_key id(active_comm_pattern_keys[i].pattern_index,active_comm_pattern_keys[i].tag,active_comm_pattern_keys[i].comm_size,
                        active_comm_pattern_keys[i].comm_color,active_comm_pattern_keys[i].msg_size,active_comm_pattern_keys[i].partner_offset); 
    if (comm_pattern_map.find(id) == comm_pattern_map.end()){
      comm_pattern_map[id] = pattern_key_id(true, i, active_comm_pattern_keys[i].pattern_index, false);
    } else{
      comm_pattern_map[id].is_active = true;	// always assumed true
      comm_pattern_map[id].key_index = i;
      comm_pattern_map[id].val_index = active_comm_pattern_keys[i].pattern_index;
      comm_pattern_map[id].is_updated = true;
    }
    if (rank != info_receiver[num_critical_path_measures-1].second){
      if ((id.tag > 12) && (id.tag < 19)){
        if (id == comm_key){
          key_match = true;
          active_patterns[comm_pattern_map[id].val_index] = local_pattern;// update with local data that must have been added once before.
        }
      }
    }
  }
  // If not key match, then we need to add the key ourselves
  if (!key_match){
    // It is assumed that 'comm_pattern_map' stores our local 'comm_key' as a key, as we don't delete keys until loop below
    active_patterns.push_back(local_pattern);
    comm_key.pattern_index = active_patterns.size()-1;
    active_comm_pattern_keys.push_back(comm_key);
    comm_pattern_map[comm_key].is_active = true;	// always assumed true
    comm_pattern_map[comm_key].key_index = active_comm_pattern_keys.size()-1;
    comm_pattern_map[comm_key].val_index = comm_key.pattern_index;
    comm_pattern_map[comm_key].is_updated = true;
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
      comp_pattern_map[id] = pattern_key_id(true, i, active_comp_pattern_keys[i].pattern_index,false);
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
  if (analysis_mode==3){// Autotuning using critical path analysis requires knowledge of which rank determined the cp
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
}

void path::propagate(nonblocking& tracker){}

// No overload on decltype(tracker) because this exchange occurs only with global communication
// Note: this exchange does not need to exchange computation kernels in the active pathset, because
//       it does not prevent feasibility of the idea. Deadlock is not possible.
// Note: entire data structures, rather than just the is_active member, are required to be exchanged,
//       as otherwise, it wouldn't be feasible: receiver would have no way of verifying whether the sent
//       order matches its own order.
// Note: per-process statistical analysis of kernels occurs only within a phase.
//       at phase-end, we exchange to pdate kernel state in a global way to prevent deadlock in scheduling comm patterns.
void path::exchange_patterns_per_process(blocking& tracker){
  assert(analysis_mode==1);// Only per-process statistical kernel data allowed.
  assert(tracker.comm == MPI_COMM_WORLD);
  int size; MPI_Comm_size(tracker.comm,&size);
  int rank; MPI_Comm_rank(tracker.comm,&rank);

  // Update the state of each communication kernel. The computation kernels were able to update their state eagerly.
  for (auto it : comm_pattern_map){
    bool is_steady = steady_test(it.first,it.second,comm_analysis_param);
    set_kernel_state(it.second,!is_steady);// Only set state, not global state!
  }
  std::vector<comm_pattern_key> foreign_active_comm_pattern_keys = active_comm_pattern_keys;
  std::vector<pattern> foreign_active_patterns = active_patterns;

  size_t active_size = size;
  size_t active_rank = rank;
  size_t active_mult = 1;
  while (active_size>1){
    if (active_rank % 2 == 1){
      int partner = (active_rank-1)*active_mult;
      int size_array[2] = {foreign_active_comm_pattern_keys.size(),foreign_active_patterns.size()};
      // Send sizes before true message so that receiver can be aware of the array sizes for subsequent communication
      PMPI_Send(&size_array[0],2,MPI_INT,partner,internal_tag,tracker.comm);
      // Send active patterns with keys
      PMPI_Send(&foreign_active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,partner,internal_tag2,tracker.comm);
      PMPI_Send(&foreign_active_patterns[0],size_array[1],pattern_type,partner,internal_tag2,tracker.comm);
      break;// Incredibely important. Senders must not update {active_size,active_rank,active_mult}
    }
    else if ((active_rank % 2 == 0) && (active_rank < (active_size-1))){
      int partner = (active_rank+1)*active_mult;
      int size_array[2] = {0,0};
      // Recv sizes of arrays to create buffers for subsequent communication
      PMPI_Recv(&size_array[0],2,MPI_INT,partner,internal_tag,tracker.comm,MPI_STATUS_IGNORE);
      // Recv partner's active patterns with keys
      foreign_active_comm_pattern_keys.resize(size_array[0]);
      foreign_active_patterns.resize(size_array[1]);
      PMPI_Recv(&foreign_active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      PMPI_Recv(&foreign_active_patterns[0],size_array[1],pattern_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      /* Iterate over all active patterns and simply perform an AND operation on whether a pattern is in steady state.
         If just one is active across the world communicator, the kernel must remain active.
         If kernel does not exist among the sent patterns, it does not count as active. The logical operation is a trivial (AND 1) */
      for (auto i=0; i<foreign_active_comm_pattern_keys.size(); i++){
        comm_pattern_key id(foreign_active_comm_pattern_keys[i].pattern_index,foreign_active_comm_pattern_keys[i].tag,foreign_active_comm_pattern_keys[i].comm_size,
                            foreign_active_comm_pattern_keys[i].comm_color,foreign_active_comm_pattern_keys[i].msg_size,foreign_active_comm_pattern_keys[i].partner_offset); 
        if (comm_pattern_map.find(id) != comm_pattern_map.end()){
          // Note: I'd like an assert here that both my local kernel is active and that matches the received kernel of the same signature, but I have no access to that foreign data member.
          // Must check whether both the foreign process's kernel data and the recv process's kernel data are sufficiently predictable
          //   (but still set to active, because active -> inactive (global steady state) occurs officially in flush_patterns)
          // I think I can still get away with simply checking each pattern's 'steady_state' member
          /* Iterate over all foreign active patterns and simply update the local active pattern steady_state members for each
            If just one is active across the world communicator, the kernel must remain active.
            If kernel does not exist among the sent patterns, it does not count as active. The logical operation is a trivial (AND 1) */
          foreign_active_patterns[foreign_active_comm_pattern_keys[i].pattern_index].steady_state &= active_patterns[comm_pattern_map[id].val_index].steady_state;
          comm_pattern_map[id].is_updated=true;
        }
      }
      // Add any local patterns that were not in foreign_active_patterns (use set notation)
      for (auto& it : comm_pattern_map){
        if (!it.second.is_updated){
          // Add the unknown kernel into our data structures to keep the loop invariant, as receiver may become sender in a later iteration
          foreign_active_patterns.push_back(active_patterns[it.second.val_index]);
          foreign_active_comm_pattern_keys.push_back(it.first);
          foreign_active_comm_pattern_keys[foreign_active_comm_pattern_keys.size()-1].pattern_index = foreign_active_patterns.size()-1;
        } else{
          it.second.is_updated=false;// reset
        }
      }
    }
    active_size = active_size/2 + active_size%2;
    active_rank /= 2;
    active_mult *= 2;
  }
  // Broadcast final exchanged kernel statistics
  int size_array[2] = {rank==0 ? foreign_active_comm_pattern_keys.size() : 0,
                       rank==0 ? foreign_active_patterns.size() : 0};
  PMPI_Bcast(&size_array[0],2,MPI_INT,0,tracker.comm);
  if (rank != 0){
    foreign_active_comm_pattern_keys.resize(size_array[0]);
    foreign_active_patterns.resize(size_array[1]);
  }
  PMPI_Bcast(&foreign_active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,0,tracker.comm);
  PMPI_Bcast(&foreign_active_patterns[0],size_array[1],pattern_type,0,tracker.comm);
  for (auto i=0; i<foreign_active_comm_pattern_keys.size(); i++){
    comm_pattern_key id(foreign_active_comm_pattern_keys[i].pattern_index,foreign_active_comm_pattern_keys[i].tag,foreign_active_comm_pattern_keys[i].comm_size,
                        foreign_active_comm_pattern_keys[i].comm_color,foreign_active_comm_pattern_keys[i].msg_size,foreign_active_comm_pattern_keys[i].partner_offset); 
    // If this process doesn't have a key matching that of the foreigner's, then we simply add in the key to our local data structures.
    // The logic here is that not doing so would allow potential deadlock in the next phase, if say there was a new communication routine with a matching signature to that
    //   already specified as a key. Some processes might recognize the key while others might not, causing deadlock in some cases.
    // Further, there is no harm in doing so, as per-process statistical analysis is still observed. We just get some other process's recorded kernel data.
    //   In a subsequent phase, a kernel's data might contain multiple process's data, but not really concerned with that. Its rare that the deadlock case would occur anyway.
    if (comm_pattern_map.find(id) == comm_pattern_map.end()){
      active_patterns.push_back(foreign_active_patterns[foreign_active_comm_pattern_keys[i].pattern_index]);
      foreign_active_comm_pattern_keys[i].pattern_index = active_patterns.size()-1;
      active_comm_pattern_keys.push_back(foreign_active_comm_pattern_keys[i]);
      comm_pattern_map[id] = pattern_key_id(true, active_comm_pattern_keys.size()-1, active_patterns.size()-1, false);
    } else{
      // Simply update the steady_state member (global_steady_state member and is_active cannot be updated until flush_patterns for non-optimized per-process statistical analysis
      active_patterns[comm_pattern_map[id].val_index].steady_state = foreign_active_patterns[foreign_active_comm_pattern_keys[i].pattern_index].steady_state;
    }
  }
}

// No overload on decltype(tracker) because this exchange occurs only with global communication
// Note: this exchange does need to exchange computation kernels in the active pathset, because
//       beyond detecting and changing kernel state, it combines statistical data (volumetric statistical analysis) to accelerate convergence to steady state
// Note: entire data structures, rather than just the is_active member, are required to be exchanged,
//       as otherwise, it wouldn't be feasible: receiver would have no way of verifying whether the sent
//       order matches its own order.
// Note: volumetric statistical analysis of kernels occurs only within a phase.
//       at phase-end, we exchange to pdate kernel state in a global way to prevent deadlock in scheduling comm patterns.
void path::exchange_patterns_volumetric(blocking& tracker){
  assert(analysis_mode==2);// Only volumetric statistical kernel data allowed.
  assert(tracker.comm == MPI_COMM_WORLD);

  int size; MPI_Comm_size(tracker.comm,&size);
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  std::vector<comm_pattern_key> foreign_active_comm_pattern_keys = active_comm_pattern_keys;
  std::vector<comp_pattern_key> foreign_active_comp_pattern_keys = active_comp_pattern_keys;
  std::vector<pattern> foreign_active_patterns = active_patterns;

  size_t active_size = size;
  size_t active_rank = rank;
  size_t active_mult = 1;
  while (active_size>1){
    if (active_rank % 2 == 1){
      int partner = (active_rank-1)*active_mult;
      int size_array[3] = {foreign_active_comm_pattern_keys.size(),foreign_active_comp_pattern_keys.size(),foreign_active_patterns.size()};
      // Send sizes before true message so that receiver can be aware of the array sizes for subsequent communication
      PMPI_Send(&size_array[0],3,MPI_INT,partner,internal_tag,tracker.comm);
      // Send active patterns with keys
      PMPI_Send(&foreign_active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,partner,internal_tag2,tracker.comm);
      PMPI_Send(&foreign_active_comp_pattern_keys[0],size_array[1],comp_pattern_key_type,partner,internal_tag2,tracker.comm);
      PMPI_Send(&foreign_active_patterns[0],size_array[2],pattern_type,partner,internal_tag2,tracker.comm);
      break;// Incredibely important. Senders must not update {active_size,active_rank,active_mult}
    }
    else if ((active_rank % 2 == 0) && (active_rank < (active_size-1))){
      int partner = (active_rank+1)*active_mult;
      int size_array[3] = {0,0,0};
      // Recv sizes of arrays to create buffers for subsequent communication
      PMPI_Recv(&size_array[0],3,MPI_INT,partner,internal_tag,tracker.comm,MPI_STATUS_IGNORE);
      // Recv partner's active patterns with keys
      foreign_active_comm_pattern_keys.resize(size_array[0]);
      foreign_active_comp_pattern_keys.resize(size_array[1]);
      foreign_active_patterns.resize(size_array[2]);
      PMPI_Recv(&foreign_active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      PMPI_Recv(&foreign_active_comp_pattern_keys[0],size_array[1],comp_pattern_key_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      PMPI_Recv(&foreign_active_patterns[0],size_array[2],pattern_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      /* Iterate over all active patterns and simply perform an AND operation on whether a pattern is in steady state.
         If just one is active across the world communicator, the kernel must remain active.
         If kernel does not exist among the sent patterns, it does not count as active. The logical operation is a trivial (AND 1) */
      for (auto i=0; i<foreign_active_comm_pattern_keys.size(); i++){
        comm_pattern_key id(foreign_active_comm_pattern_keys[i].pattern_index,foreign_active_comm_pattern_keys[i].tag,foreign_active_comm_pattern_keys[i].comm_size,
                            foreign_active_comm_pattern_keys[i].comm_color,foreign_active_comm_pattern_keys[i].msg_size,foreign_active_comm_pattern_keys[i].partner_offset); 
        if (comm_pattern_map.find(id) != comm_pattern_map.end()){
          update_kernel_stats(foreign_active_patterns[foreign_active_comm_pattern_keys[i].pattern_index],active_patterns[comm_pattern_map[id].val_index],comm_analysis_param);
          bool is_steady = steady_test(id,foreign_active_patterns[foreign_active_comm_pattern_keys[i].pattern_index],comm_analysis_param);
          set_kernel_state(foreign_active_patterns[foreign_active_comm_pattern_keys[i].pattern_index],!is_steady);
          comm_pattern_map[id].is_updated=true;
        }
      }
      for (auto i=0; i<foreign_active_comp_pattern_keys.size(); i++){
        comp_pattern_key id(foreign_active_comp_pattern_keys[i].pattern_index,foreign_active_comp_pattern_keys[i].tag,foreign_active_comp_pattern_keys[i].flops,
                            foreign_active_comp_pattern_keys[i].param1,foreign_active_comp_pattern_keys[i].param2,foreign_active_comp_pattern_keys[i].param3,
                            foreign_active_comp_pattern_keys[i].param4,foreign_active_comp_pattern_keys[i].param5); 
        if (comp_pattern_map.find(id) != comp_pattern_map.end()){
          update_kernel_stats(foreign_active_patterns[foreign_active_comp_pattern_keys[i].pattern_index],active_patterns[comp_pattern_map[id].val_index],comp_analysis_param);
          bool is_steady = steady_test(id,foreign_active_patterns[foreign_active_comp_pattern_keys[i].pattern_index],comp_analysis_param);
          set_kernel_state(foreign_active_patterns[foreign_active_comp_pattern_keys[i].pattern_index],!is_steady);
          comp_pattern_map[id].is_updated=true;
        }
      }
      // Add any local patterns that were not in foreign_active_patterns (use set notation)
      for (auto& it : comm_pattern_map){
        if (!it.second.is_updated){
          // Add the unknown kernel into our data structures to keep the loop invariant, as receiver may become sender in a later iteration
          foreign_active_patterns.push_back(active_patterns[it.second.val_index]);
          foreign_active_comm_pattern_keys.push_back(it.first);
          foreign_active_comm_pattern_keys[foreign_active_comm_pattern_keys.size()-1].pattern_index = foreign_active_patterns.size()-1;
        } else{
          it.second.is_updated=false;// reset
        }
      }
      for (auto& it : comp_pattern_map){
        if (!it.second.is_updated){
          // Add the unknown kernel into our data structures to keep the loop invariant, as receiver may become sender in a later iteration
          foreign_active_patterns.push_back(active_patterns[it.second.val_index]);
          foreign_active_comp_pattern_keys.push_back(it.first);
          foreign_active_comp_pattern_keys[foreign_active_comp_pattern_keys.size()-1].pattern_index = foreign_active_patterns.size()-1;
        } else{
          it.second.is_updated=false;// reset
        }
      }
    }
    active_size = active_size/2 + active_size%2;
    active_rank /= 2;
    active_mult *= 2;
  }
  // Broadcast final exchanged kernel statistics
  int size_array[3] = {rank==0 ? foreign_active_comm_pattern_keys.size() : 0,
                       rank==0 ? foreign_active_comp_pattern_keys.size() : 0,
                       rank==0 ? foreign_active_patterns.size() : 0};
  PMPI_Bcast(&size_array[0],3,MPI_INT,0,tracker.comm);
  if (rank != 0){
    foreign_active_comm_pattern_keys.resize(size_array[0]);
    foreign_active_comp_pattern_keys.resize(size_array[1]);
    foreign_active_patterns.resize(size_array[2]);
  }
  PMPI_Bcast(&foreign_active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,0,tracker.comm);
  PMPI_Bcast(&foreign_active_comp_pattern_keys[0],size_array[1],comp_pattern_key_type,0,tracker.comm);
  PMPI_Bcast(&foreign_active_patterns[0],size_array[2],pattern_type,0,tracker.comm);
  for (auto i=0; i<foreign_active_comm_pattern_keys.size(); i++){
    comm_pattern_key id(foreign_active_comm_pattern_keys[i].pattern_index,foreign_active_comm_pattern_keys[i].tag,foreign_active_comm_pattern_keys[i].comm_size,
                        foreign_active_comm_pattern_keys[i].comm_color,foreign_active_comm_pattern_keys[i].msg_size,foreign_active_comm_pattern_keys[i].partner_offset); 
    // If this process doesn't have a key matching that of the foreigner's, then we simply add in the key to our local data structures.
    // The logic here is that not doing so would allow potential deadlock in the next phase, if say there was a new communication routine with a matching signature to that
    //   already specified as a key. Some processes might recognize the key while others might not, causing deadlock in some cases.
    // Further, there is no harm in doing so, as volumetric statistical analysis is still observed.
    if (comm_pattern_map.find(id) == comm_pattern_map.end()){
      active_patterns.push_back(foreign_active_patterns[foreign_active_comm_pattern_keys[i].pattern_index]);
      foreign_active_comm_pattern_keys[i].pattern_index = active_patterns.size()-1;
      active_comm_pattern_keys.push_back(foreign_active_comm_pattern_keys[i]);
      comm_pattern_map[id] = pattern_key_id(true, active_comm_pattern_keys.size()-1, active_patterns.size()-1, false);
    } else{
      // Update existing entry.
      active_patterns[comm_pattern_map[id].val_index] = foreign_active_patterns[foreign_active_comm_pattern_keys[i].pattern_index];
    }
  }
  for (auto i=0; i<foreign_active_comp_pattern_keys.size(); i++){
    comp_pattern_key id(foreign_active_comp_pattern_keys[i].pattern_index,foreign_active_comp_pattern_keys[i].tag,foreign_active_comp_pattern_keys[i].flops,
                        foreign_active_comp_pattern_keys[i].param1,foreign_active_comp_pattern_keys[i].param2,foreign_active_comp_pattern_keys[i].param3,
                        foreign_active_comp_pattern_keys[i].param4,foreign_active_comp_pattern_keys[i].param5); 
    // If this process doesn't have a key matching that of the foreigner's, then we simply add in the key to our local data structures.
    // The logic here is that not doing so would allow potential deadlock in the next phase, if say there was a new communication routine with a matching signature to that
    //   already specified as a key. Some processes might recognize the key while others might not, causing deadlock in some cases.
    // Further, there is no harm in doing so, as volumetric statistical analysis is still observed.
    if (comp_pattern_map.find(id) == comp_pattern_map.end()){
      active_patterns.push_back(foreign_active_patterns[foreign_active_comp_pattern_keys[i].pattern_index]);
      foreign_active_comp_pattern_keys[i].pattern_index = active_patterns.size()-1;
      active_comp_pattern_keys.push_back(foreign_active_comp_pattern_keys[i]);
      comp_pattern_map[id] = pattern_key_id(true, active_comp_pattern_keys.size()-1, active_patterns.size()-1, false);
    } else{
      // Update existing entry.
      active_patterns[comp_pattern_map[id].val_index] = foreign_active_patterns[foreign_active_comp_pattern_keys[i].pattern_index];
    }
  }
}

}
}
}
