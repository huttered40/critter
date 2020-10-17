#include "path.h"
#include "../container/symbol_tracker.h"
#include "../util/util.h"
#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

void path::exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm){
  // Save and accumulate the computation time between last communication routine as both execution-time and computation time
  //   into both the execution-time critical path data structures and the per-process data structures.
  double save_comp_time = MPI_Wtime() - computation_timer;
  critical_path_costs[num_critical_path_measures-1] += save_comp_time;	// update critical path execution time
  critical_path_costs[num_critical_path_measures-3] += save_comp_time;	// update critical path computation time
  volume_costs[num_volume_measures-1]        += save_comp_time;		// update local execution time
  volume_costs[num_volume_measures-3]        += save_comp_time;		// update local computation time

  generate_aggregate_channels(oldcomm,newcomm);
  PMPI_Barrier(oldcomm);
  computation_timer = MPI_Wtime();
}

bool path::initiate_comp(std::vector<intptr_t>& user_array, size_t id, volatile double curtime, double flop_count, int param1, int param2, int param3, int param4, int param5){
  // Save and accumulate the computation time between last communication routine as both execution-time and computation time
  //   into both the execution-time critical path data structures and the per-process data structures.
  double save_comp_time = curtime - computation_timer;
  critical_path_costs[num_critical_path_measures-1] += save_comp_time;	// update critical path execution time
  critical_path_costs[num_critical_path_measures-3] += save_comp_time;	// update critical path computation time
  volume_costs[num_volume_measures-1]        += save_comp_time;		// update local execution time
  volume_costs[num_volume_measures-3]        += save_comp_time;		// update local computation time
  // Special exit if no kernels are to be scheduled -- the goal is to track the total overhead time (no comp/comm kernels), which should
  //   be attained with timers outside of critter.
  if (schedule_kernels==0){ return false; }
  volatile double overhead_start_time = MPI_Wtime();

  bool schedule_decision = true;
  if (!(tuning_delta==0 || tuning_delta==2)){// these tuning_deltas (0, 2) signify unconditional scheduling of computation kernels
    comp_pattern_key key(-1,id,flop_count,param1,param2,param3,param4,param5);// '-1' argument is arbitrary, does not influence overloaded operators
    // Below, the idea is that key doesn't exist in comp_pattern_map iff the key hasn't been seen before. If the key has been seen, we automatically
    //   create an entry in comp_pattern_key, although it will be empty.
    if (comp_pattern_map.find(key) != comp_pattern_map.end()){
      schedule_decision = should_schedule(comp_pattern_map[key])==1;
    }
    if (schedule_decision && comp_kernel_buffer_id){
      for (auto it : user_array){
        if (comp_pattern_map.find(key) != comp_pattern_map.end()){// Note I could also put an assert here
          if (skip_ptr_set.find(it) != skip_ptr_set.end()){
            schedule_decision = false;
            set_kernel_state(comp_pattern_map[key],false);
            if (comp_state_aggregation_mode==0) set_kernel_state_global(comp_pattern_map[key],false);
            skip_ptr_set.erase(it);
          }
        }
      }
    }
  }
  intercept_overhead[0] += MPI_Wtime() - overhead_start_time;
  // start compunication timer for compunication routine
  comp_start_time = MPI_Wtime();
  return schedule_decision;
}

void path::complete_comp(std::vector<intptr_t>& user_array, size_t id, double flop_count, int param1, int param2, int param3, int param4, int param5){
  volatile double comp_time = MPI_Wtime( ) - comp_start_time;	// complete computation time
  // Special exit if no kernels are to be scheduled -- the goal is to track the total overhead time (no comp/comm kernels), which should
  //   be attained with timers outside of critter.
  if (schedule_kernels==0){ return; }
  volatile double overhead_start_time = MPI_Wtime();

  if (!(tuning_delta==0 || tuning_delta==2)){// these deltas (0, 2) signify unconditional scheduling of computation kernels
    comp_pattern_key key(active_patterns.size(),id,flop_count,param1,param2,param3,param4,param5);// 'active_patterns.size()' argument is arbitrary, does not influence overloaded operators
    // Below, the idea is that key doesn't exist in comp_pattern_map iff the key hasn't been seen before. If the key has been seen, we automatically
    //   create an entry in comp_pattern_key, although it will be empty.
    if (comp_pattern_map.find(key) == comp_pattern_map.end()){
      active_comp_pattern_keys.push_back(key);
      active_patterns.emplace_back();
      comp_pattern_map[key] = pattern_key_id(true,active_comp_pattern_keys.size()-1,active_patterns.size()-1,false);
    }
    if ((comp_sample_aggregation_mode==0) || (should_schedule(comp_pattern_map[key]) == 0)){
      // Note that we can enter here even if comp_state_aggregation_mode>0
      // Because of the second case in the branch, this will be called even if in sample_aggregation_mode>=1 as long as the kernel is no longer being scheduled (i.e. is in steady state)
      update_kernel_stats(comp_pattern_map[key],comp_analysis_param,comp_time,flop_count);
    }
    else if (comp_sample_aggregation_mode == 1){
      // Its assumed that this branch is entered iff the kernel was scheduled and comm_time is reliable
      if (comp_batch_map.find(key) == comp_batch_map.end()){
        comp_batch_map[key].emplace_back();// initialize an empty pattern_batch
      }

      bool found_batch=false;
      auto& batch_list = comp_batch_map[key];
      int index=0;
      while (index < batch_list.size()){
        // Look for the first inactive batch and aggregate new sample with it
        if (batch_list[index].hash_id==0){// inactive state
          found_batch=true;
          break;
        } else{
          index++;
        }
      }
      if (!found_batch){
        batch_list.emplace_back();
        index = batch_list.size()-1;
      }
      update_kernel_stats(batch_list[index],comp_analysis_param,comp_time,flop_count);
    }
    // Note: 'get_estimate' must be called before setting the updated kernel state. If kernel was not scheduled, comp_time set below overwrites 'comp_time'
    if (should_schedule(comp_pattern_map[key]) == 0){
      if (comp_sample_aggregation_mode == 1) assert(comp_batch_map[key].size() == 0);
      comp_time = get_estimate(comp_pattern_map[key],comp_analysis_param,flop_count);
    } else{
      bool is_steady = steady_test(key,comp_pattern_map[key],comp_analysis_param);
      set_kernel_state(comp_pattern_map[key],!is_steady);
      if (comp_state_aggregation_mode==0) set_kernel_state_global(comp_pattern_map[key],!is_steady);
    }
  }

  critical_path_costs[num_critical_path_measures-1] += comp_time;	// execution time
  critical_path_costs[num_critical_path_measures-2] += comp_time;	// computation kernel time
  critical_path_costs[num_critical_path_measures-3] += comp_time;	// computational time
  volume_costs[num_volume_measures-1] += comp_time;			// execution time
  volume_costs[num_volume_measures-2] += comp_time;			// computation kernel time
  volume_costs[num_volume_measures-3] += comp_time;			// computation time

  intercept_overhead[0] += MPI_Wtime() - overhead_start_time;
  computation_timer = MPI_Wtime();
}

bool path::initiate_comm(std::vector<intptr_t>& user_msg, blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                         bool is_sender, int partner1, int partner2){
  // Save and accumulate the computation time between last communication routine as both execution-time and computation time
  //   into both the execution-time critical path data structures and the per-process data structures.
  tracker.comp_time = curtime - computation_timer;
  critical_path_costs[num_critical_path_measures-1] += tracker.comp_time;	// update critical path execution time
  critical_path_costs[num_critical_path_measures-3] += tracker.comp_time;	// update critical path computation time
  volume_costs[num_volume_measures-1]        += tracker.comp_time;		// update local execution time
  volume_costs[num_volume_measures-3]        += tracker.comp_time;		// update local computation time
  // Special exit if no kernels are to be scheduled -- the goal is to track the total overhead time (no comp/comm kernels), which should
  //   be attained with timers outside of critter.
  if (schedule_kernels==0){ return false; }

  // At this point, 'critical_path_costs' tells me the process's time up until now. A barrier won't suffice when kernels are conditionally scheduled.
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  int rank; MPI_Comm_rank(comm, &rank);
  volatile double overhead_start_time = MPI_Wtime();

  // We consider usage of Sendrecv variants to forfeit usage of eager internal communication.
  // Note that the reason we can't force user Bsends to be 'true_eager_p2p' is because the corresponding Recvs would be expecting internal communications
  bool true_eager_p2p = ((eager_p2p == 1) && (tracker.tag!=13) && (tracker.tag!=14));
  if (true_eager_p2p){ MPI_Buffer_attach(&eager_pad[0],eager_pad.size()); }

  // Save caller communication attributes into reference object for use in corresponding static method 'complete_comm'
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
  tracker.aggregate_comp_patterns=false;
  tracker.aggregate_comm_patterns=false;

  // Non-optimized variant will always post barriers, although of course, just as with the optimized variant, the barriers only remove idle time
  //   from corrupting communication time measurements. The process that enters barrier last is not necessarily the critical path root. The
  //     critical path root is decided based on a reduction using 'critical_path_costs'. Therefore, no explicit barriers are invoked, instead relying on Allreduce
  bool schedule_decision = true;
  double reduced_info[17] = {critical_path_costs[num_critical_path_measures-4],critical_path_costs[num_critical_path_measures-3],
                            critical_path_costs[num_critical_path_measures-2],critical_path_costs[num_critical_path_measures-1],0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,bsp_counter};
  double reduced_info_foreign[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  // Assume that the communicator of either collective/p2p is registered via comm_split, and that its described using a max of 3 dimension tuples.
  assert(comm_channel_map.find(tracker.comm) != comm_channel_map.end());
  int comm_sizes[2]={0,0}; int comm_strides[2]={0,0};
  for (auto i=0; i<comm_channel_map[tracker.comm]->id.size(); i++){
    comm_sizes[i]=comm_channel_map[tracker.comm]->id[i].first;
    comm_strides[i]=comm_channel_map[tracker.comm]->id[i].second;
  }
  // Below, the idea is that key doesn't exist in comm_pattern_map iff the key hasn't been seen before. If the key has been seen, we automatically
  //   create an entry in comm_pattern_key, although it will be empty.
  comm_pattern_key key(rank,-1,tracker.tag,comm_sizes,comm_strides,tracker.nbytes,tracker.partner1);
  if (tuning_delta > 1){// tuning_delta of 1 signifies that communication kernel scheduling is to be unconditional
    if (comm_pattern_map.find(key) != comm_pattern_map.end()){
      schedule_decision = should_schedule(comm_pattern_map[key])==1;
      reduced_info[4] = (!schedule_decision ? 1 : 0);// set to 1 if kernel is globally steady.
      if (!schedule_decision){
        // If this particular kernel is globally steady, meaning it has exhausted its state aggregation channels,
        //   then we can overwrite the '-1' with the sample mean of the globally-steady kernel
        reduced_info[7] = active_patterns[comm_pattern_map[key].val_index].hash_id;
        reduced_info[8] = active_patterns[comm_pattern_map[key].val_index].num_schedules;
        reduced_info[9] = active_patterns[comm_pattern_map[key].val_index].num_local_schedules;
        reduced_info[10] = active_patterns[comm_pattern_map[key].val_index].num_scheduled_units;
        reduced_info[11] = active_patterns[comm_pattern_map[key].val_index].num_local_scheduled_units;
        reduced_info[12] = active_patterns[comm_pattern_map[key].val_index].M1;
        reduced_info[13] = active_patterns[comm_pattern_map[key].val_index].M2;
        reduced_info[14] = active_patterns[comm_pattern_map[key].val_index].total_exec_time;
        reduced_info[15] = active_patterns[comm_pattern_map[key].val_index].total_local_exec_time;
      }
    }
  }
  // Register the p2p channel
  if (tracker.partner1 != -1){// p2p
    int my_world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&my_world_rank);
    auto world_partner_rank = channel::translate_rank(tracker.comm,tracker.partner1);
    if (p2p_channel_map.find(world_partner_rank) == p2p_channel_map.end()){
      solo_channel* node = new solo_channel();
      node->tag = world_partner_rank;//key.partner_offset;
      node->offset = std::min(my_world_rank,world_partner_rank);
      node->id.push_back(std::make_pair(2,abs(my_world_rank-world_partner_rank)));
      std::string local_channel_hash_str = ".." + std::to_string(node->id[0].first) + "." + std::to_string(node->id[0].second);
      std::string global_channel_hash_str = ".." + std::to_string(node->id[0].first) + "." + std::to_string(node->id[0].second);
      node->local_hash_tag = std::hash<std::string>()(local_channel_hash_str);
      node->global_hash_tag = std::hash<std::string>()(global_channel_hash_str);
      //spf.insert_node(node);
      p2p_channel_map[world_partner_rank] = node;
      // Build up aggregates: DO NOT ADD THE P2P CHANNEL ALONE
      // This channel is basically going to be defunct: lots of logic is pulled out because, for now at least, an aggregate with a p2p channel embedded
      //   cannot add another p2p channel. This holds for single-stage aggregation, but not sure if it does for multi-stage aggregation
      // Check if 'node' is a sibling of all existing aggregates already formed. Note that we do not include p2p aggregates, nor p2p+comm aggregates.
      // Note this loop assumes that the local_hash_tags of each aggregate across new_comm are in the same sorted order (hence the assert below)
      for (auto& it : aggregate_channel_map){
        if (it.second->offset == -1) continue; // For now at least, p2p channels cannot form aggregates with aggregates already built with a p2p channel on top.
        // 0. Check that each process in newcomm is processing the same aggregate -------> SKIP
        // 1. Check if 'node' is a child of 'aggregate'
        bool is_child_1 = channel::verify_ancestor_relation(it.second,node);
        // 2. Check if 'aggregate' is a child of 'node'
        bool is_child_2 = channel::verify_ancestor_relation(node,it.second);
        // 3. Check if 'node'+'aggregate' form a sibling
        bool is_sibling = channel::verify_sibling_relation(it.second,node);
        if (is_sibling && !is_child_1 && !is_child_2){
          // DO NOT RESET If current aggregate forms a larger one with 'node', as we did when registering communicators
          // it.second->is_final = false;
          int new_local_hash_tag = it.second->local_hash_tag ^ node->local_hash_tag;
          int new_global_hash_tag = it.second->global_hash_tag ^ node->global_hash_tag;
          auto new_aggregate_channel = new aggregate_channel(it.second->id,new_local_hash_tag,new_global_hash_tag,0,it.second->num_channels+1);// '0' gets updated below
          // Set the hashes of each communicator.
          new_aggregate_channel->channels.insert(node->local_hash_tag);
          for (auto it_2 : it.second->channels){
            new_aggregate_channel->channels.insert(it_2);
          }
          // Do not communicate to attain the minimum offset of all process in newcomm's aggregate channel -> simply set offset to -1 to recognize an aggregate that is not to be built on top of.
          new_aggregate_channel->offset = -1;
        }
      }
    }
  }

  // Use pathsets, not batches, to check if kernel can leverage an aggregation. Such a kernel must be locally steady (i.e.
  //   from its own schedules, its steady), and must be able to aggregate across the channel associated with 'tracker.comm'
  if ((tracker.partner1 == -1) && (comp_state_aggregation_mode>0)){
    for (auto& it : comp_pattern_map){
      if (!((active_patterns[it.second.val_index].steady_state==1) && (should_schedule(it.second)==1))) continue;
      // Any global communicator can fast-track a communication kernel to being in global steady state. No need to match up hashes (that would only be necessary for sample aggregation)
      if (aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->is_final){
        tracker.save_comp_key.push_back(it.first);
        reduced_info[5]++;
      }
      else{
        if (active_patterns[it.second.val_index].registered_channels.find(comm_channel_map[tracker.comm]) != active_patterns[it.second.val_index].registered_channels.end()) continue;
        // TODO: Not exactly sure whether to use global_hash_id below or local_hash_id
        if (aggregate_channel_map.find(active_patterns[it.second.val_index].hash_id ^ aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag) == aggregate_channel_map.end()) continue;
        tracker.save_comp_key.push_back(it.first);
        reduced_info[5]++;
      }
    }
  }
  if ((tracker.partner1 == -1) && (comm_state_aggregation_mode>0)){
    for (auto& it : comm_pattern_map){
      if (!((active_patterns[it.second.val_index].steady_state==1) && (should_schedule(it.second)==1))) continue;
      if (aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->is_final){
        tracker.save_comm_key.push_back(it.first);
        reduced_info[6]++;
      }
      else{
        if (active_patterns[it.second.val_index].registered_channels.find(comm_channel_map[tracker.comm]) != active_patterns[it.second.val_index].registered_channels.end()) continue;
        // TODO: Not exactly sure whether to use global_hash_id below or local_hash_id
        if (aggregate_channel_map.find(active_patterns[it.second.val_index].hash_id ^ aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag) == aggregate_channel_map.end()) continue;
        tracker.save_comm_key.push_back(it.first);
        reduced_info[6]++;
      }
    }
  }

  tracker.should_propagate = false;
  assert(partner1 != MPI_ANY_SOURCE); if ((tracker.tag == 13) || (tracker.tag == 14)){ assert(partner2 != MPI_ANY_SOURCE); }
  if (partner1 == -1){
    PMPI_Allreduce(MPI_IN_PLACE, &reduced_info[0], 17, MPI_DOUBLE, MPI_MAX, tracker.comm);
    tracker.barrier_time = reduced_info[3] - critical_path_costs[num_critical_path_measures-1];
    for (auto i=0; i<num_critical_path_measures; i++){ critical_path_costs[i] = reduced_info[i]; }
    schedule_decision = (reduced_info[4] == 0 ? true : false);
    tracker.aggregate_comp_patterns = reduced_info[5]>0;
    tracker.aggregate_comm_patterns = reduced_info[6]>0;
    tracker.should_propagate = reduced_info[5]>0 || reduced_info[6]>0;
    bsp_counter = reduced_info[16];
    if (comm_pattern_map.find(key) != comm_pattern_map.end()){
      if (!schedule_decision){
        if (should_schedule(comm_pattern_map[key])){
          if (comm_kernel_transfer_id == 1){
            // Enter here if a process's local comm kernel is not globally steady, yet its found that at least one of the processors in its communicator is.
            // Completely swap out its pattern statistics for the elements reduced. Note that I am avoiding an extra broadcast, thus the reduced members might
            //   each be from different processors. Likely though, just one is globally steady.
            active_patterns[comm_pattern_map[key].val_index].hash_id = reduced_info[7];
            active_patterns[comm_pattern_map[key].val_index].num_schedules = reduced_info[8];
            active_patterns[comm_pattern_map[key].val_index].num_local_schedules = reduced_info[9];
            active_patterns[comm_pattern_map[key].val_index].num_scheduled_units = reduced_info[10];
            active_patterns[comm_pattern_map[key].val_index].num_local_scheduled_units = reduced_info[11];
            active_patterns[comm_pattern_map[key].val_index].M1 = reduced_info[12];
            active_patterns[comm_pattern_map[key].val_index].M2 = reduced_info[13];
            active_patterns[comm_pattern_map[key].val_index].total_exec_time = reduced_info[14];
            active_patterns[comm_pattern_map[key].val_index].total_local_exec_time  = reduced_info[15];
          }
        }
        set_kernel_state(comm_pattern_map[key],false);
        set_kernel_state_global(comm_pattern_map[key],false);
      }
    }
  }
  else{
    if ((true_eager_p2p) && (rank != tracker.partner1)){
      assert(0);// Not implemented yet. Does eager protocol even make sense, given these methods?
    }
    else if (!true_eager_p2p){
      PMPI_Sendrecv(&reduced_info[0], 17, MPI_DOUBLE, tracker.partner1, internal_tag2, &reduced_info_foreign[0], 17,
                    MPI_DOUBLE, tracker.partner2, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
      tracker.barrier_time = std::max(reduced_info[3],reduced_info_foreign[3]) - critical_path_costs[num_critical_path_measures-1];
      for (auto i=0; i<num_critical_path_measures; i++){ critical_path_costs[i] = std::max(reduced_info[i],reduced_info_foreign[i]); }
      schedule_decision = (reduced_info[4] == 0) && (reduced_info_foreign[4] == 0);// If either are 1, then we don't schedule
      tracker.should_propagate = false;// We currently don't use p2p channels to aggregate
      bsp_counter = std::max(reduced_info[16],reduced_info_foreign[16]);
      if (comm_pattern_map.find(key) != comm_pattern_map.end()){
        if (!schedule_decision){
          if (should_schedule(comm_pattern_map[key])){
            // Enter here if a process's local comm kernel is not globally steady, yet its found that at least one of the processors in its communicator is.
            // Completely swap out its pattern statistics for the elements reduced. Note that I am avoiding an extra broadcast, thus the reduced members might
            //   each be from different processors. Likely though, just one is globally steady.
            active_patterns[comm_pattern_map[key].val_index].hash_id = reduced_info_foreign[7];
            active_patterns[comm_pattern_map[key].val_index].num_schedules = reduced_info_foreign[8];
            active_patterns[comm_pattern_map[key].val_index].num_local_schedules = reduced_info_foreign[9];
            active_patterns[comm_pattern_map[key].val_index].num_scheduled_units = reduced_info_foreign[10];
            active_patterns[comm_pattern_map[key].val_index].num_local_scheduled_units = reduced_info_foreign[11];
            active_patterns[comm_pattern_map[key].val_index].M1 = reduced_info_foreign[12];
            active_patterns[comm_pattern_map[key].val_index].M2 = reduced_info_foreign[13];
            active_patterns[comm_pattern_map[key].val_index].total_exec_time = reduced_info_foreign[14];
            active_patterns[comm_pattern_map[key].val_index].total_local_exec_time  = reduced_info_foreign[15];
          }
          set_kernel_state(comm_pattern_map[key],false);
          set_kernel_state_global(comm_pattern_map[key],false);
        }
      }
      if (tracker.partner2 != tracker.partner1){
        assert(0);// TODO: Fix later
        // This if-statement will never be breached if 'true_eager_p2p'=true anyways.
        for (auto i=0; i<num_critical_path_measures; i++){ reduced_info[i] = critical_path_costs[i]; }
        PMPI_Sendrecv(&reduced_info[0], 17, MPI_DOUBLE, tracker.partner2, internal_tag2, &reduced_info_foreign[0], 17,
                      MPI_DOUBLE, tracker.partner1, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
        tracker.barrier_time = std::max(reduced_info[3],reduced_info_foreign[3]) - critical_path_costs[num_critical_path_measures-1];
        for (auto i=0; i<num_critical_path_measures; i++){ critical_path_costs[i] = std::max(reduced_info[i],reduced_info_foreign[i]); }
        schedule_decision = schedule_decision && (reduced_info[4] == 0) && (reduced_info_foreign[4] == 0);
        tracker.should_propagate = tracker.should_propagate || reduced_info[5]>0 || reduced_info_foreign[5]>0;
      }
    }
  }


  // If kernel is about to be scheduled, post one more barrier for safety if collective,
  //   because the AllReduce posted above may allow ranks to leave early, thus corrupting the sample measurement.
  if (schedule_decision && tracker.partner1 == -1) PMPI_Barrier(tracker.comm);
  else if (schedule_decision){
    char ch1='H',ch2='H';
    PMPI_Sendrecv(&ch1, 1, MPI_CHAR, tracker.partner1, internal_tag2, &ch2, 1,
                  MPI_CHAR, tracker.partner2, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
  }

  intercept_overhead[1] += MPI_Wtime() - overhead_start_time;
  if (!schedule_decision){
    for (auto it : user_msg){
      skip_ptr_set.insert(it);
    }
  }
  // start communication timer for communication routine
  tracker.start_time = MPI_Wtime();
  return schedule_decision;
}

// Used only for p2p communication. All blocking collectives use sychronous protocol
void path::complete_comm(std::vector<intptr_t>& user_msg, blocking& tracker, int recv_source){
  volatile double comm_time = MPI_Wtime() - tracker.start_time;	// complete communication time
  // Special exit if no kernels are to be scheduled -- the goal is to track the total overhead time (no comp/comm kernels), which should
  //   be attained with timers outside of critter.
  if (schedule_kernels==0){ return; }
  volatile double overhead_start_time = MPI_Wtime();

  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  bool true_eager_p2p = ((eager_p2p == 1) && (tracker.tag!=13) && (tracker.tag!=14));
  // We handle wildcard sources (for MPI_Recv variants) only after the user communication.
  if (recv_source != -1){
    if ((tracker.tag == 13) || (tracker.tag == 14)){ tracker.partner2=recv_source; }
    else{ assert(tracker.tag==17); tracker.partner1=recv_source; }
  }

  if (tuning_delta > 1){// tuning_delta of 1 signifies that communication kernel scheduling is to be unconditional
    int comm_sizes[2]={0,0}; int comm_strides[2]={0,0};
    for (auto i=0; i<comm_channel_map[tracker.comm]->id.size(); i++){
      comm_sizes[i]=comm_channel_map[tracker.comm]->id[i].first;
      comm_strides[i]=comm_channel_map[tracker.comm]->id[i].second;
    }
    // Below, the idea is that key doesn't exist in comm_pattern_map iff the key hasn't been seen before. If the key has been seen, we automatically
    //   create an entry in comm_pattern_key, although it will be empty.
    comm_pattern_key key(rank,active_patterns.size(),tracker.tag,comm_sizes,comm_strides,tracker.nbytes,tracker.partner1);
    //if (world_rank == 8) std::cout << "******" << key.tag << " " << key.dim_sizes[0] << " " << key.dim_strides[0] << " " << key.msg_size << " " << key.partner_offset << ")\n";
    if (comm_pattern_map.find(key) == comm_pattern_map.end()){
      active_comm_pattern_keys.push_back(key);
      if (tracker.partner1 != -1){
        auto world_partner_rank = channel::translate_rank(tracker.comm,tracker.partner1);
        assert(p2p_channel_map.find(world_partner_rank) != p2p_channel_map.end());
        active_patterns.emplace_back(p2p_channel_map[world_partner_rank]);
      } else{
        assert(comm_channel_map.find(tracker.comm) != comm_channel_map.end());
        active_patterns.emplace_back(comm_channel_map[tracker.comm]);
      }
      comm_pattern_map[key] = pattern_key_id(true,active_comm_pattern_keys.size()-1,active_patterns.size()-1,false);
    }
    int comm_hash_tag;
    if ((comm_sample_aggregation_mode==0) || (should_schedule(comm_pattern_map[key]) == 0)){
      update_kernel_stats(comm_pattern_map[key],comm_analysis_param,comm_time,tracker.nbytes);
    }
    else if (comm_sample_aggregation_mode >= 1){
      assert(should_schedule(comm_pattern_map[key]) == 1);
      // Note: its common for a key's entry into the batch_map to be deleted many times (especially aggregation strategy #1)
      if (comm_batch_map.find(key) == comm_batch_map.end()){
        // Register the channel's batches
        if (tracker.partner1 != -1){
          auto world_partner_rank = channel::translate_rank(tracker.comm,tracker.partner1);
          assert(p2p_channel_map.find(world_partner_rank) != p2p_channel_map.end());
          comm_batch_map[key].emplace_back(p2p_channel_map[world_partner_rank]);
        } else{
          assert(comm_channel_map.find(tracker.comm) != comm_channel_map.end());
          comm_batch_map[key].emplace_back(comm_channel_map[tracker.comm]);
        }
      }
      bool found_batch=false;
      auto& batch_list = comm_batch_map[key];
      int index=0;
      while (index < batch_list.size()){
        // Look for the first inactive batch and aggregate new sample with it
        if (batch_list[index].channel_count == 0){// inactive state
          found_batch=true;
          break;
        } else{
          index++;
        }
      }
      if (!found_batch){
        if (tracker.partner1 != -1){
          auto world_partner_rank = channel::translate_rank(tracker.comm,tracker.partner1);
          batch_list.push_back(pattern_batch(p2p_channel_map[world_partner_rank]));
        } else{
          batch_list.push_back(pattern_batch(comm_channel_map[tracker.comm]));
        }
        index = batch_list.size()-1;
      }
      update_kernel_stats(batch_list[index],comm_analysis_param,comm_time,tracker.nbytes);
    }
    if (should_schedule(comm_pattern_map[key])==0){
      // Note if this is true, the corresponding entry in the batch map must be cleared. However, I think I delete the entire map in aggregation mode 1, so asserting
      //   on this is difficult.
      // Note: I think branching on aggregation mode is not needed. The pathset should contain all batch samples and the batch should be cleared.
      comm_time = get_estimate(comm_pattern_map[key],comm_analysis_param,tracker.nbytes);
      if (comm_sample_aggregation_mode == 1) assert(comm_batch_map[key].size() == 0);
    } else{
      bool is_steady = steady_test(key,comm_pattern_map[key],comm_analysis_param);
      set_kernel_state(comm_pattern_map[key],!is_steady);
      if (sample_constraint_mode == -1) set_kernel_state_global(comm_pattern_map[key],!is_steady);// Force global state to steady.
      // The line below will force p2ps to change from unsteady to globally steady immediately, thus preventing need to aggregate.
      // We can allow this temporarily, but cannot allow this for collectives.
      if ((tracker.partner1 != -1) || (comm_state_aggregation_mode == 0)) { set_kernel_state_global(comm_pattern_map[key],!is_steady); }
    }
  }

  critical_path_costs[num_critical_path_measures-1] += comm_time;
  critical_path_costs[num_critical_path_measures-4] += comm_time;
  volume_costs[num_volume_measures-1] += (comm_time + tracker.barrier_time);
  volume_costs[num_volume_measures-4] += comm_time;

  //volume_costs[num_volume_measures-4] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
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
  if (tracker.should_propagate && tracker.partner1 == -1){
    bool is_world_communication = (tracker.comm == MPI_COMM_WORLD) && (tracker.partner1 == -1);
    if ((rank == tracker.partner1) && (rank == tracker.partner2)) { ; }
    else{
      if (comm_sample_aggregation_mode==1){
        if (comm_state_aggregation_mode==1){
          single_stage_sample_aggregation(tracker);
        } else if (comm_state_aggregation_mode==2){
          multi_stage_sample_aggregation(tracker);
        }
      }
      else{
        if (tracker.aggregate_comm_patterns) comm_state_aggregation(tracker);
        if (tracker.aggregate_comp_patterns) comp_state_aggregation(tracker);
      }
    }
  }
  if (true_eager_p2p){
    void* temp_buf; int temp_size;
    // Forces buffered messages to send. Ideally we should wait till the next invocation of 'path::initiate(blocking&,...)' to call this,
    //   but to be safe and avoid stalls caused by MPI implementation not sending until this routine is called, we call it here.
    MPI_Buffer_detach(&temp_buf,&temp_size);
  }

  tracker.should_propagate = false;
  tracker.aggregate_comp_patterns = false;
  tracker.aggregate_comm_patterns = false;
  bsp_counter++;
  intercept_overhead[2] += MPI_Wtime() - overhead_start_time;
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

void path::comp_state_aggregation(blocking& tracker){
  int size; MPI_Comm_size(tracker.comm,&size);
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  std::vector<pattern_propagate> foreign_active_patterns;
  std::map<comp_pattern_key,pattern_propagate> save_comp_kernels;

  // First save the kernels we want to contribute to the aggregation (because they are steady)
  for (auto& it : tracker.save_comp_key){
    save_comp_kernels[it] = active_patterns[comp_pattern_map[it].val_index];
  }

  size_t active_size = size;
  size_t active_rank = rank;
  size_t active_mult = 1;
  while (active_size>1){
    if (active_rank % 2 == 1){
      // Fill-in the associated comp pattern
      tracker.save_comp_key.clear();
      foreign_active_patterns.clear();
      for (auto& it : save_comp_kernels){
        tracker.save_comp_key.push_back(it.first);
        foreign_active_patterns.push_back(it.second);
      }

      int partner = (active_rank-1)*active_mult;
      int size_array[2] = {tracker.save_comp_key.size(),tracker.save_comp_key.size()};
      // Send sizes before true message so that receiver can be aware of the array sizes for subsequent communication
      PMPI_Send(&size_array[0],2,MPI_INT,partner,internal_tag,tracker.comm);
      // Send active patterns with keys
      PMPI_Send(&tracker.save_comp_key[0],size_array[0],comp_pattern_key_type,partner,internal_tag2,tracker.comm);
      PMPI_Send(&foreign_active_patterns[0],size_array[1],pattern_type,partner,internal_tag2,tracker.comm);
      break;// Incredibely important. Senders must not update {active_size,active_rank,active_mult}
    }
    else if ((active_rank % 2 == 0) && (active_rank < (active_size-1))){
      int partner = (active_rank+1)*active_mult;
      int size_array[2] = {0,0};
      // Recv sizes of arrays to create buffers for subsequent communication
      PMPI_Recv(&size_array[0],2,MPI_INT,partner,internal_tag,tracker.comm,MPI_STATUS_IGNORE);
      // Recv partner's active patterns with keys
      tracker.save_comp_key.resize(size_array[0]);
      foreign_active_patterns.resize(size_array[1]);
      PMPI_Recv(&tracker.save_comp_key[0],size_array[0],comp_pattern_key_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      PMPI_Recv(&foreign_active_patterns[0],size_array[1],pattern_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      // Iterate over all active patterns and simply perform an AND operation on whether a pattern is in steady state.
      //   If just one is active across the world communicator, the kernel must remain active.
      //   If kernel does not exist among the sent patterns, it does not count as active. The logical operation is a trivial (AND 1)
      for (auto i=0; i<tracker.save_comp_key.size(); i++){
        auto& key = tracker.save_comp_key[i];
        if (save_comp_kernels.find(key) != save_comp_kernels.end()){
          double ci_local = get_error_estimate(save_comp_kernels[key],comp_analysis_param);
          double ci_foreign = get_error_estimate(foreign_active_patterns[i],comp_analysis_param);
          if (ci_foreign < ci_local){
            save_comp_kernels[key] = foreign_active_patterns[i];
          }
        } else{
          save_comp_kernels[key] = foreign_active_patterns[i];
        }
      }
    }
    active_size = active_size/2 + active_size%2;
    active_rank /= 2;
    active_mult *= 2;
  }
  // Broadcast final exchanged kernel statistics
  if (rank==0){
    tracker.save_comp_key.clear();
    foreign_active_patterns.clear();
    for (auto& it : save_comp_kernels){
      tracker.save_comp_key.push_back(it.first);
      foreign_active_patterns.push_back(it.second);
    }
  }
  int size_array[2] = {rank==0 ? tracker.save_comp_key.size() : 0,
                       rank==0 ? tracker.save_comp_key.size() : 0};
  PMPI_Bcast(&size_array[0],2,MPI_INT,0,tracker.comm);
  if (rank != 0){
    tracker.save_comp_key.resize(size_array[0]);
    foreign_active_patterns.resize(size_array[1]);
  }
  PMPI_Bcast(&tracker.save_comp_key[0],size_array[0],comp_pattern_key_type,0,tracker.comm);
  PMPI_Bcast(&foreign_active_patterns[0],size_array[1],pattern_type,0,tracker.comm);
  for (auto i=0; i<tracker.save_comp_key.size(); i++){
    auto& key = tracker.save_comp_key[i];
    if (comp_pattern_map.find(key) != comp_pattern_map.end()){
      if (comp_kernel_transfer_id==0){
        active_patterns[comp_pattern_map[key].val_index].hash_id = foreign_active_patterns[i].hash_id;
      }
      else if (comp_kernel_transfer_id==1){
        active_patterns[comp_pattern_map[key].val_index] = foreign_active_patterns[i];
      }
    } else{
      // Add new entry.
      active_comp_pattern_keys.push_back(key);
      active_patterns.emplace_back(foreign_active_patterns[i]);
      comp_pattern_map[key] = pattern_key_id(true,active_comp_pattern_keys.size()-1,active_patterns.size()-1,false);
    }
    set_kernel_state(comp_pattern_map[key],false);
    if (comp_state_aggregation_mode==1){
      set_kernel_state_global(comp_pattern_map[key],false);
    } else if (comp_state_aggregation_mode==2){
      if (aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->is_final){
        active_patterns[comp_pattern_map[key].val_index].hash_id = aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag;
        active_patterns[comp_pattern_map[key].val_index].registered_channels.clear();
        active_patterns[comp_pattern_map[key].val_index].registered_channels.insert(comm_channel_map[tracker.comm]);// Add the solo channel, not the aggregate
        assert(aggregate_channel_map.find(active_patterns[comp_pattern_map[key].val_index].hash_id) != aggregate_channel_map.end());
        set_kernel_state_global(comp_pattern_map[key],false);
      }
      else{
        active_patterns[comp_pattern_map[key].val_index].hash_id ^= aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag;
        active_patterns[comp_pattern_map[key].val_index].registered_channels.insert(comm_channel_map[tracker.comm]);// Add the solo channel, not the aggregate
        assert(aggregate_channel_map.find(active_patterns[comp_pattern_map[key].val_index].hash_id) != aggregate_channel_map.end());
        if (aggregate_channel_map[active_patterns[comp_pattern_map[key].val_index].hash_id]->is_final){
          set_kernel_state_global(comp_pattern_map[key],false);
        }
      }
    }
  }
  tracker.save_comp_key.clear();
}

void path::comm_state_aggregation(blocking& tracker){
  int size; MPI_Comm_size(tracker.comm,&size);
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  std::vector<pattern_propagate> foreign_active_patterns;
  std::map<comm_pattern_key,pattern_propagate> save_comm_kernels;

  // First save the kernels we want to contribute to the aggregation (because they are steady)
  for (auto& it : tracker.save_comm_key){
    save_comm_kernels[it] = active_patterns[comm_pattern_map[it].val_index];
  }

  size_t active_size = size;
  size_t active_rank = rank;
  size_t active_mult = 1;
  while (active_size>1){
    if (active_rank % 2 == 1){
      // Fill-in the associated comm pattern
      tracker.save_comm_key.clear();
      foreign_active_patterns.clear();
      for (auto& it : save_comm_kernels){
        tracker.save_comm_key.push_back(it.first);
        foreign_active_patterns.push_back(it.second);
      }

      int partner = (active_rank-1)*active_mult;
      int size_array[2] = {tracker.save_comm_key.size(),tracker.save_comm_key.size()};
      // Send sizes before true message so that receiver can be aware of the array sizes for subsequent communication
      PMPI_Send(&size_array[0],2,MPI_INT,partner,internal_tag,tracker.comm);
      // Send active patterns with keys
      PMPI_Send(&tracker.save_comm_key[0],size_array[0],comm_pattern_key_type,partner,internal_tag2,tracker.comm);
      PMPI_Send(&foreign_active_patterns[0],size_array[1],pattern_type,partner,internal_tag2,tracker.comm);
      break;// Incredibely important. Senders must not update {active_size,active_rank,active_mult}
    }
    else if ((active_rank % 2 == 0) && (active_rank < (active_size-1))){
      int partner = (active_rank+1)*active_mult;
      int size_array[2] = {0,0};
      // Recv sizes of arrays to create buffers for subsequent communication
      PMPI_Recv(&size_array[0],2,MPI_INT,partner,internal_tag,tracker.comm,MPI_STATUS_IGNORE);
      // Recv partner's active patterns with keys
      tracker.save_comm_key.resize(size_array[0]);
      foreign_active_patterns.resize(size_array[1]);
      PMPI_Recv(&tracker.save_comm_key[0],size_array[0],comm_pattern_key_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      PMPI_Recv(&foreign_active_patterns[0],size_array[1],pattern_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      // Iterate over all active patterns and simply perform an AND operation on whether a pattern is in steady state.
      //   If just one is active across the world communicator, the kernel must remain active.
      //   If kernel does not exist among the sent patterns, it does not count as active. The logical operation is a trivial (AND 1)
      for (auto i=0; i<tracker.save_comm_key.size(); i++){
        auto& key = tracker.save_comm_key[i];
        if (save_comm_kernels.find(key) != save_comm_kernels.end()){
          double ci_local = get_error_estimate(save_comm_kernels[key],comm_analysis_param);
          double ci_foreign = get_error_estimate(foreign_active_patterns[i],comm_analysis_param);
          if (ci_foreign < ci_local){
            save_comm_kernels[key] = foreign_active_patterns[i];
          }
        } else{
          save_comm_kernels[key] = foreign_active_patterns[i];
        }
      }
    }
    active_size = active_size/2 + active_size%2;
    active_rank /= 2;
    active_mult *= 2;
  }
  // Broadcast final exchanged kernel statistics
  if (rank==0){
    tracker.save_comm_key.clear();
    foreign_active_patterns.clear();
    for (auto& it : save_comm_kernels){
      tracker.save_comm_key.push_back(it.first);
      foreign_active_patterns.push_back(it.second);
    }
  }
  int size_array[2] = {rank==0 ? tracker.save_comm_key.size() : 0,
                       rank==0 ? tracker.save_comm_key.size() : 0};
  PMPI_Bcast(&size_array[0],2,MPI_INT,0,tracker.comm);
  if (rank != 0){
    tracker.save_comm_key.resize(size_array[0]);
    foreign_active_patterns.resize(size_array[1]);
  }
  PMPI_Bcast(&tracker.save_comm_key[0],size_array[0],comm_pattern_key_type,0,tracker.comm);
  PMPI_Bcast(&foreign_active_patterns[0],size_array[1],pattern_type,0,tracker.comm);
  for (auto i=0; i<tracker.save_comm_key.size(); i++){
    auto& key = tracker.save_comm_key[i];
    if (comm_pattern_map.find(key) != comm_pattern_map.end()){
      if (comm_kernel_transfer_id==0){
        active_patterns[comm_pattern_map[key].val_index].hash_id = foreign_active_patterns[i].hash_id;
      }
      else if (comm_kernel_transfer_id==1){
        active_patterns[comm_pattern_map[key].val_index] = foreign_active_patterns[i];
      }
    } else{
      // Add new entry.
      active_comm_pattern_keys.push_back(key);
      active_patterns.emplace_back(foreign_active_patterns[i]);
      comm_pattern_map[key] = pattern_key_id(true,active_comm_pattern_keys.size()-1,active_patterns.size()-1,false);
    }
    set_kernel_state(comm_pattern_map[key],false);
    if (comm_state_aggregation_mode==1){
      set_kernel_state_global(comm_pattern_map[key],false);
    } else if (comm_state_aggregation_mode==2){
      if (aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->is_final){
        active_patterns[comm_pattern_map[key].val_index].hash_id = aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag;
        active_patterns[comm_pattern_map[key].val_index].registered_channels.clear();
        active_patterns[comm_pattern_map[key].val_index].registered_channels.insert(comm_channel_map[tracker.comm]);// Add the solo channel, not the aggregate
        assert(aggregate_channel_map.find(active_patterns[comm_pattern_map[key].val_index].hash_id) != aggregate_channel_map.end());
        set_kernel_state_global(comm_pattern_map[key],false);
      }
      else{
        active_patterns[comm_pattern_map[key].val_index].hash_id ^= aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag;
        active_patterns[comm_pattern_map[key].val_index].registered_channels.insert(comm_channel_map[tracker.comm]);// Add the solo channel, not the aggregate
        assert(aggregate_channel_map.find(active_patterns[comm_pattern_map[key].val_index].hash_id) != aggregate_channel_map.end());
        if (aggregate_channel_map[active_patterns[comm_pattern_map[key].val_index].hash_id]->is_final){
          set_kernel_state_global(comm_pattern_map[key],false);
        }
      }
    }
  }
  tracker.save_comm_key.clear();
}

void path::single_stage_sample_aggregation(blocking& tracker){
  assert(0);// Shut off for now
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  int comm_rank; MPI_Comm_rank(tracker.comm,&comm_rank);
  int comm_size; MPI_Comm_size(tracker.comm,&comm_size);
  std::vector<comm_pattern_key> foreign_active_comm_pattern_keys;
  std::vector<comp_pattern_key> foreign_active_comp_pattern_keys;
  std::vector<pattern_batch_propagate> foreign_active_batches;

  size_t active_size = comm_size;
  size_t active_rank = comm_rank;
  size_t active_mult = 1;
  while (active_size>1){
    if (active_rank % 2 == 1){
      foreign_active_comm_pattern_keys.clear();
      foreign_active_comp_pattern_keys.clear();
      foreign_active_batches.clear();
      // Local batches can only be added into batch list once. As a process never sends twice, below is a reasonable place.
      for (auto& it : comm_batch_map){
        assert(it.second.size() < 2);// Number of states/batches of a particular key must be < 2. This is almost trivially satisfied.
        if (it.second.size() == 1){// Some comm keys will have no active batches/states, so this check is mandatory
          // Three checks below:
          //   1. Does this propagation channel match the channel of the particular batch?
          //   2. Does this propagation channel form an aggregate with the channel of the particular batch?
          //   3. Does this batch have a stride of 0? This indicates a trivial p2p that does not need sending
          if (it.second[0].registered_channels.find(comm_channel_map[tracker.comm]) != it.second[0].registered_channels.end()) continue;
          // TODO: Not exactly sure whether to use global_hash_id below or local_hash_id
          if (aggregate_channel_map.find(it.second[0].hash_id ^ aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag) == aggregate_channel_map.end()) continue;
          // if ((it.second[0].hash_id == 1) && (it.second[0].id[0].second == 0)) continue;
          foreign_active_comm_pattern_keys.push_back(it.first);
          foreign_active_batches.emplace_back(it.second[0]);
        }
      }
      for (auto& it : comp_batch_map){
        assert(it.second.size() < 2);
        if (it.second.size() == 1){
          foreign_active_comp_pattern_keys.push_back(it.first);
          foreign_active_batches.emplace_back(it.second[0]);
        }
      }
      int partner = (active_rank-1)*active_mult;
      int size_array[3] = {foreign_active_comm_pattern_keys.size(),foreign_active_comp_pattern_keys.size(),foreign_active_batches.size()};
      // Send sizes before true message so that receiver can be aware of the array sizes for subsequent communication
      PMPI_Send(&size_array[0],3,MPI_INT,partner,internal_tag,tracker.comm);
      // Send active patterns with keys
      PMPI_Send(&foreign_active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,partner,internal_tag2,tracker.comm);
      PMPI_Send(&foreign_active_comp_pattern_keys[0],size_array[1],comp_pattern_key_type,partner,internal_tag2,tracker.comm);
      PMPI_Send(&foreign_active_batches[0],size_array[2],batch_type,partner,internal_tag2,tracker.comm);
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
      foreign_active_batches.resize(size_array[2]);
      PMPI_Recv(&foreign_active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      PMPI_Recv(&foreign_active_comp_pattern_keys[0],size_array[1],comp_pattern_key_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      PMPI_Recv(&foreign_active_batches[0],size_array[2],batch_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      // Iterate over each received batch and incproporate it into its running list of batches in the other arrays
      for (auto i=0; i<foreign_active_comm_pattern_keys.size(); i++){
        comm_pattern_key id(foreign_active_comm_pattern_keys[i].pattern_index,foreign_active_comm_pattern_keys[i].tag,foreign_active_comm_pattern_keys[i].dim_sizes,
                            foreign_active_comm_pattern_keys[i].dim_strides,foreign_active_comm_pattern_keys[i].msg_size,foreign_active_comm_pattern_keys[i].partner_offset); 
        pattern_batch temp; temp.M1 = foreign_active_batches[i].M1; temp.M2 = foreign_active_batches[i].M2; temp.hash_id = foreign_active_batches[i].hash_id;
                            temp.channel_count = foreign_active_batches[i].channel_count; temp.num_schedules = foreign_active_batches[i].num_schedules;
                            temp.num_scheduled_units = foreign_active_batches[i].num_scheduled_units; temp.total_exec_time = foreign_active_batches[i].total_exec_time;
        if (comm_batch_map.find(id) != comm_batch_map.end()){
          assert(comm_batch_map[id].size() < 2);
          if (comm_batch_map[id].size() == 0){
            // No need updating the registered_channels member because this aggregation strategy liquidates batches immediately.
            comm_batch_map[id].push_back(temp);
          } else{
            update_kernel_stats(comm_batch_map[id][0],temp,comm_analysis_param);
          }
        } else{
          // No need updating the registered_channels member because this aggregation strategy liquidates batches immediately.
          comm_batch_map[id].push_back(temp);
        }
      }
      for (auto i=0; i<foreign_active_comp_pattern_keys.size(); i++){
        comp_pattern_key id(foreign_active_comp_pattern_keys[i].pattern_index,foreign_active_comp_pattern_keys[i].tag,foreign_active_comp_pattern_keys[i].flops,
                            foreign_active_comp_pattern_keys[i].param1,foreign_active_comp_pattern_keys[i].param2,foreign_active_comp_pattern_keys[i].param3,
                            foreign_active_comp_pattern_keys[i].param4,foreign_active_comp_pattern_keys[i].param5); 
        int j = size_array[0]+i;
        pattern_batch temp; temp.M1 = foreign_active_batches[j].M1; temp.M2 = foreign_active_batches[j].M2; temp.hash_id = foreign_active_batches[j].hash_id;
                            temp.channel_count = foreign_active_batches[j].channel_count; temp.num_schedules = foreign_active_batches[j].num_schedules;
                            temp.num_scheduled_units = foreign_active_batches[j].num_scheduled_units; temp.total_exec_time = foreign_active_batches[j].total_exec_time;
        if (comp_batch_map.find(id) != comp_batch_map.end()){
          assert(comp_batch_map[id].size() < 2);
          if (comp_batch_map[id].size() == 0){
            // No need updating the registered_channels member because this aggregation strategy liquidates batches immediately.
            comp_batch_map[id].push_back(temp);
          } else{
            update_kernel_stats(comp_batch_map[id][0],temp,comp_analysis_param);
          }
        } else{
          // No need updating the registered_channels member because this aggregation strategy liquidates batches immediately.
          comp_batch_map[id].push_back(temp);
        }
      }
    }
    active_size = active_size/2 + active_size%2;
    active_rank /= 2;
    active_mult *= 2;
  }
  // Broadcast final exchanged kernel statistics
  if (comm_rank==0){
    // Clear these three arrays before pushing back with new entries
    foreign_active_comm_pattern_keys.clear();
    foreign_active_comp_pattern_keys.clear();
    foreign_active_batches.clear(); 
    for (auto& it : comm_batch_map){
      assert(it.second.size() < 2);
      if (it.second.size() == 1){
        // We still need special checks here, because this root will likely not have sent to anyone.
        // Three checks below:
        //   1. Does this propagation channel match the channel of the particular batch?
        //   2. Does this propagation channel form an aggregate with the channel of the particular batch?
        //   3. Does this batch have a stride of 0? This indicates a trivial p2p that does not need sending
        if (it.second[0].registered_channels.find(comm_channel_map[tracker.comm]) != it.second[0].registered_channels.end()) continue;
        // TODO: Not exactly sure whether to use global_hash_id below or local_hash_id
        if (aggregate_channel_map.find(it.second[0].hash_id ^ aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag) == aggregate_channel_map.end()) continue;
        // if ((it.second[0].hash_id == 1) && (it.second[0].id[0].second == 0)) continue;
        foreign_active_comm_pattern_keys.push_back(it.first);
        foreign_active_batches.push_back(it.second[0]);
      }
    }
    for (auto& it : comp_batch_map){
      assert(it.second.size() < 2);
      if (it.second.size() == 1){
        foreign_active_comp_pattern_keys.push_back(it.first);
        foreign_active_batches.push_back(it.second[0]);
      }
    }
  }
  int size_array[3] = {comm_rank==0 ? foreign_active_comm_pattern_keys.size() : 0,
                       comm_rank==0 ? foreign_active_comp_pattern_keys.size() : 0,
                       comm_rank==0 ? foreign_active_batches.size() : 0};
  PMPI_Bcast(&size_array[0],3,MPI_INT,0,tracker.comm);
  if (comm_rank != 0){
    foreign_active_comm_pattern_keys.resize(size_array[0]);
    foreign_active_comp_pattern_keys.resize(size_array[1]);
    foreign_active_batches.resize(size_array[2]);
  }
  PMPI_Bcast(&foreign_active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,0,tracker.comm);
  PMPI_Bcast(&foreign_active_comp_pattern_keys[0],size_array[1],comp_pattern_key_type,0,tracker.comm);
  PMPI_Bcast(&foreign_active_batches[0],size_array[2],batch_type,0,tracker.comm);
  // Now incorporate this not into comm_batch_map, but into comm_pattern_map. This is the reason why we allow root rank 0 to enter these loops
  for (auto i=0; i<foreign_active_comm_pattern_keys.size(); i++){
    comm_pattern_key id(foreign_active_comm_pattern_keys[i].pattern_index,foreign_active_comm_pattern_keys[i].tag,foreign_active_comm_pattern_keys[i].dim_sizes,
                        foreign_active_comm_pattern_keys[i].dim_strides,foreign_active_comm_pattern_keys[i].msg_size,foreign_active_comm_pattern_keys[i].partner_offset); 
    pattern_batch temp; temp.M1 = foreign_active_batches[i].M1; temp.M2 = foreign_active_batches[i].M2; temp.hash_id = foreign_active_batches[i].hash_id;
                        temp.channel_count = foreign_active_batches[i].channel_count; temp.num_schedules = foreign_active_batches[i].num_schedules;
                        temp.num_scheduled_units = foreign_active_batches[i].num_scheduled_units; temp.total_exec_time = foreign_active_batches[i].total_exec_time;
    if (comm_pattern_map.find(id) == comm_pattern_map.end()){
      active_patterns.emplace_back(temp);
      foreign_active_comm_pattern_keys[i].pattern_index = active_patterns.size()-1;
      active_comm_pattern_keys.push_back(foreign_active_comm_pattern_keys[i]);
      comm_pattern_map[id] = pattern_key_id(true, active_comm_pattern_keys.size()-1, active_patterns.size()-1, false);
    }
    if (comm_batch_map.find(id) != comm_batch_map.end()){
      if (comm_batch_map[id].size() > 0){
        temp.num_local_schedules = comm_batch_map[id][0].num_local_schedules;
        temp.num_local_scheduled_units = comm_batch_map[id][0].num_local_scheduled_units;
        temp.total_local_exec_time = comm_batch_map[id][0].total_local_exec_time;
      }
      comm_batch_map[id].clear();// the samples residing here have been propagated and the batch's journey is complete.
    }
    // Update existing entry.
    bool is_steady = steady_test(id,comm_pattern_map[id],comm_analysis_param);
    bool is_global_steady = should_schedule(comm_pattern_map[id]);
    if (!is_steady && is_global_steady){
      update_kernel_stats(active_patterns[comm_pattern_map[id].val_index], temp, comm_analysis_param);
      bool is_steady = steady_test(id,comm_pattern_map[id],comm_analysis_param);
      set_kernel_state(comm_pattern_map[id],!is_steady);
    }
  }
  for (auto i=0; i<foreign_active_comp_pattern_keys.size(); i++){
    comp_pattern_key id(foreign_active_comp_pattern_keys[i].pattern_index,foreign_active_comp_pattern_keys[i].tag,foreign_active_comp_pattern_keys[i].flops,
                        foreign_active_comp_pattern_keys[i].param1,foreign_active_comp_pattern_keys[i].param2,foreign_active_comp_pattern_keys[i].param3,
                        foreign_active_comp_pattern_keys[i].param4,foreign_active_comp_pattern_keys[i].param5); 
    int j = size_array[0]+i;
    pattern_batch temp; temp.M1 = foreign_active_batches[j].M1; temp.M2 = foreign_active_batches[j].M2; temp.hash_id = foreign_active_batches[j].hash_id;
                        temp.channel_count = foreign_active_batches[j].channel_count; temp.num_schedules = foreign_active_batches[j].num_schedules;
                        temp.num_scheduled_units = foreign_active_batches[j].num_scheduled_units; temp.total_exec_time = foreign_active_batches[j].total_exec_time;
    if (comp_pattern_map.find(id) == comp_pattern_map.end()){
      active_patterns.emplace_back(temp);
      foreign_active_comp_pattern_keys[i].pattern_index = active_patterns.size()-1;
      active_comp_pattern_keys.push_back(foreign_active_comp_pattern_keys[i]);
      comp_pattern_map[id] = pattern_key_id(true, active_comp_pattern_keys.size()-1, active_patterns.size()-1, false);
    }
    if (comp_batch_map.find(id) != comp_batch_map.end()){
      if (comp_batch_map[id].size() > 0){
        temp.num_local_schedules = comp_batch_map[id][0].num_local_schedules;
        temp.num_local_scheduled_units = comp_batch_map[id][0].num_local_scheduled_units;
        temp.total_local_exec_time = comp_batch_map[id][0].total_local_exec_time;
      }
      comp_batch_map[id].clear();// the samples residing here have been propagated and the batch's journey is complete.
    }
    // Update existing entry.
    bool is_steady = steady_test(id,comp_pattern_map[id],comp_analysis_param);
    bool is_global_steady = should_schedule(comp_pattern_map[id]);
    if (!is_steady && is_global_steady){
      update_kernel_stats(active_patterns[comp_pattern_map[id].val_index], temp, comp_analysis_param);
      bool is_steady = steady_test(id,comp_pattern_map[id],comp_analysis_param);
      set_kernel_state(comp_pattern_map[id],!is_steady);
    }
  }
}
void path::multi_stage_sample_aggregation(blocking& tracker){
  assert(0);// Shut off for now
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  int comm_rank; MPI_Comm_rank(tracker.comm,&comm_rank);
  int comm_size; MPI_Comm_size(tracker.comm,&comm_size);
  std::vector<comm_pattern_key> foreign_active_comm_pattern_keys;
  std::vector<comp_pattern_key> foreign_active_comp_pattern_keys;
  std::vector<pattern_batch_propagate> foreign_active_batches;

  std::map<comm_pattern_key,std::vector<pattern_batch>> temp_comm_batch_map;
  std::map<comp_pattern_key,std::vector<pattern_batch>> temp_comp_batch_map;
  // Fill up temporary maps with those local batches fit to aggregate across this channel 'tracker.comm'
  for (auto& it : comm_batch_map){
    for (auto& batch_state : it.second){
      // Three checks below:
      //   1. Does this propagation channel match the channel of the particular batch, or has it been used before?
      //   2. Does this propagation channel form an aggregate with the channel of the particular batch?
      //   3. Does this batch have a stride of 0? This indicates a trivial p2p that does not need sending
      if (batch_state.registered_channels.find(comm_channel_map[tracker.comm]) != batch_state.registered_channels.end()) continue;
      // TODO: Not exactly sure whether to use global_hash_id below or local_hash_id
      if (aggregate_channel_map.find(batch_state.hash_id ^ aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag) == aggregate_channel_map.end()) continue;
//      if ((it.second[0].hash_id == 1) && (it.second[0].id[0].second == 0)) continue;
      temp_comm_batch_map[it.first].push_back(batch_state);
    }
  }
  for (auto& it : comp_batch_map){
    for (auto& batch_state : it.second){
      // Three checks below:
      //   1. Has this propagation channel been used before?
      //   2. Does this propagation channel form an aggregate with the channel of the particular batch?
      if (batch_state.registered_channels.find(comm_channel_map[tracker.comm]) != batch_state.registered_channels.end()) continue;
      // TODO: Not exactly sure whether to use global_hash_id below or local_hash_id
      if (aggregate_channel_map.find(batch_state.hash_id ^ aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag) == aggregate_channel_map.end()) continue;
      temp_comp_batch_map[it.first].push_back(batch_state);
    }
  }

  size_t active_size = comm_size;
  size_t active_rank = comm_rank;
  size_t active_mult = 1;
  while (active_size>1){
    if (active_rank % 2 == 1){
      foreign_active_comm_pattern_keys.clear();
      foreign_active_comp_pattern_keys.clear();
      foreign_active_batches.clear();
      // Local batches can only be added into batch list once. As a process never sends twice, below is a reasonable place.
      for (auto& it : temp_comm_batch_map){
        for (auto& batch_state : it.second){
          foreign_active_comm_pattern_keys.push_back(it.first);
          foreign_active_batches.emplace_back(batch_state);
        }
      }
      for (auto& it : temp_comp_batch_map){
        for (auto& batch_state : it.second){
          foreign_active_comp_pattern_keys.push_back(it.first);
          foreign_active_batches.emplace_back(batch_state);
        }
      }
      int partner = (active_rank-1)*active_mult;
      int size_array[3] = {foreign_active_comm_pattern_keys.size(),foreign_active_comp_pattern_keys.size(),foreign_active_batches.size()};
      // Send sizes before true message so that receiver can be aware of the array sizes for subsequent communication
      PMPI_Send(&size_array[0],3,MPI_INT,partner,internal_tag,tracker.comm);
      // Send active patterns with keys
      PMPI_Send(&foreign_active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,partner,internal_tag2,tracker.comm);
      PMPI_Send(&foreign_active_comp_pattern_keys[0],size_array[1],comp_pattern_key_type,partner,internal_tag2,tracker.comm);
      PMPI_Send(&foreign_active_batches[0],size_array[2],batch_type,partner,internal_tag2,tracker.comm);
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
      foreign_active_batches.resize(size_array[2]);
      PMPI_Recv(&foreign_active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      PMPI_Recv(&foreign_active_comp_pattern_keys[0],size_array[1],comp_pattern_key_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      PMPI_Recv(&foreign_active_batches[0],size_array[2],batch_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      // Iterate over each received batch and incproporate it into its running list of batches in the other arrays
      for (auto i=0; i<foreign_active_comm_pattern_keys.size(); i++){
        comm_pattern_key id(foreign_active_comm_pattern_keys[i].pattern_index,foreign_active_comm_pattern_keys[i].tag,foreign_active_comm_pattern_keys[i].dim_sizes,
                            foreign_active_comm_pattern_keys[i].dim_strides,foreign_active_comm_pattern_keys[i].msg_size,foreign_active_comm_pattern_keys[i].partner_offset); 
        pattern_batch temp; temp.M1 = foreign_active_batches[i].M1; temp.M2 = foreign_active_batches[i].M2; temp.hash_id = foreign_active_batches[i].hash_id;
                            temp.channel_count = foreign_active_batches[i].channel_count; temp.num_schedules = foreign_active_batches[i].num_schedules;
                            temp.num_scheduled_units = foreign_active_batches[i].num_scheduled_units; temp.total_exec_time = foreign_active_batches[i].total_exec_time;
        bool found_match = false;
        if (temp_comm_batch_map.find(id) != temp_comm_batch_map.end()){
          // Search for an entry in temp_comm_batch_map with the same hash_id.
          for (auto& it : temp_comm_batch_map[id]){
            if (it.hash_id == temp.hash_id){
              found_match = true;
              update_kernel_stats(it,temp,comm_analysis_param);
              break;
            }
          }
        }
        if (!found_match){
          // No need updating the registered_channels member because this aggregation strategy liquidates batches immediately.
          temp_comm_batch_map[id].push_back(temp);
        }
      }
      for (auto i=0; i<foreign_active_comp_pattern_keys.size(); i++){
        comp_pattern_key id(foreign_active_comp_pattern_keys[i].pattern_index,foreign_active_comp_pattern_keys[i].tag,foreign_active_comp_pattern_keys[i].flops,
                            foreign_active_comp_pattern_keys[i].param1,foreign_active_comp_pattern_keys[i].param2,foreign_active_comp_pattern_keys[i].param3,
                            foreign_active_comp_pattern_keys[i].param4,foreign_active_comp_pattern_keys[i].param5); 
        int j = size_array[0]+i;
        pattern_batch temp; temp.M1 = foreign_active_batches[j].M1; temp.M2 = foreign_active_batches[j].M2; temp.hash_id = foreign_active_batches[j].hash_id;
                            temp.channel_count = foreign_active_batches[j].channel_count; temp.num_schedules = foreign_active_batches[j].num_schedules;
                            temp.num_scheduled_units = foreign_active_batches[j].num_scheduled_units; temp.total_exec_time = foreign_active_batches[j].total_exec_time;
        bool found_match = false;
        if (temp_comp_batch_map.find(id) != temp_comp_batch_map.end()){
          // Search for an entry in temp_comp_batch_map with the same hash_id.
          for (auto& it : temp_comp_batch_map[id]){
            if (it.hash_id == temp.hash_id){
              found_match = true;
              update_kernel_stats(it,temp,comp_analysis_param);
              break;
            }
          }
        }
        if (!found_match){
          // No need updating the registered_channels member because this aggregation strategy liquidates batches immediately.
          temp_comp_batch_map[id].push_back(temp);
        }
      }
    }
    active_size = active_size/2 + active_size%2;
    active_rank /= 2;
    active_mult *= 2;
  }
  // Broadcast final exchanged kernel statistics
  if (comm_rank==0){
    // Clear these three arrays before pushing back with new entries
    foreign_active_comm_pattern_keys.clear();
    foreign_active_comp_pattern_keys.clear();
    foreign_active_batches.clear(); 
    for (auto& it : temp_comm_batch_map){
      for (auto& batch_state : it.second){
        // We still need special checks here, because this root will likely not have sent to anyone.
        // Three checks below:
        //   1. Does this propagation channel match the channel of the particular batch?
        //   2. Does this propagation channel form an aggregate with the channel of the particular batch?
        //   3. Does this batch have a stride of 0? This indicates a trivial p2p that does not need sending
        if (batch_state.registered_channels.find(comm_channel_map[tracker.comm]) != batch_state.registered_channels.end()) continue;
        // TODO: Not exactly sure whether to use global_hash_id below or local_hash_id
        if (aggregate_channel_map.find(batch_state.hash_id ^ aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag) == aggregate_channel_map.end()) continue;
        // if ((it.second[0].hash_id == 1) && (it.second[0].id[0].second == 0)) continue;
        foreign_active_comm_pattern_keys.push_back(it.first);
        foreign_active_batches.emplace_back(batch_state);
      }
    }
    for (auto& it : temp_comp_batch_map){
      for (auto& batch_state : it.second){
        // Three checks below:
        //   1. Has this propagation channel been used before?
        //   2. Does this propagation channel form an aggregate with the channel of the particular batch?
        if (batch_state.registered_channels.find(comm_channel_map[tracker.comm]) != batch_state.registered_channels.end()) continue;
        // TODO: Not exactly sure whether to use global_hash_id below or local_hash_id
        if (aggregate_channel_map.find(batch_state.hash_id ^ aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag) == aggregate_channel_map.end()) continue;
        foreign_active_comp_pattern_keys.push_back(it.first);
        foreign_active_batches.emplace_back(batch_state);
      }
    }
  }
  // Clear the temporary batch maps so that they can be reset solely with batches of matching states
  temp_comm_batch_map.clear();
  temp_comp_batch_map.clear();

  int size_array[3] = {comm_rank==0 ? foreign_active_comm_pattern_keys.size() : 0,
                       comm_rank==0 ? foreign_active_comp_pattern_keys.size() : 0,
                       comm_rank==0 ? foreign_active_batches.size() : 0};
  PMPI_Bcast(&size_array[0],3,MPI_INT,0,tracker.comm);
  if (comm_rank != 0){
    foreign_active_comm_pattern_keys.resize(size_array[0]);
    foreign_active_comp_pattern_keys.resize(size_array[1]);
    foreign_active_batches.resize(size_array[2]);
  }
  PMPI_Bcast(&foreign_active_comm_pattern_keys[0],size_array[0],comm_pattern_key_type,0,tracker.comm);
  PMPI_Bcast(&foreign_active_comp_pattern_keys[0],size_array[1],comp_pattern_key_type,0,tracker.comm);
  PMPI_Bcast(&foreign_active_batches[0],size_array[2],batch_type,0,tracker.comm);
  // Now incorporate this not into comm_batch_map, but into comm_pattern_map. This is the reason why we allow root rank 0 to enter these loops
  for (auto i=0; i<foreign_active_comm_pattern_keys.size(); i++){
    comm_pattern_key id(foreign_active_comm_pattern_keys[i].pattern_index,foreign_active_comm_pattern_keys[i].tag,foreign_active_comm_pattern_keys[i].dim_sizes,
                        foreign_active_comm_pattern_keys[i].dim_strides,foreign_active_comm_pattern_keys[i].msg_size,foreign_active_comm_pattern_keys[i].partner_offset); 
    // Unlike with single-stage, we are only interested in those batches that match with ours.
    // Therefore, it makes sense to check 'comm_batch_map' rather than 'temp_comm_batch_map'
    if (comm_batch_map.find(id) != comm_batch_map.end()){
      // Here, I want to search and replace the existing batch entry with matching state
      for (auto& it : comm_batch_map[id]){
        if (it.hash_id == foreign_active_batches[i].hash_id){
          it.num_schedules = foreign_active_batches[i].num_schedules;
          it.num_scheduled_units = foreign_active_batches[i].num_scheduled_units;
          it.total_exec_time = foreign_active_batches[i].total_exec_time;
          it.M1 = foreign_active_batches[i].M1;
          it.M2 = foreign_active_batches[i].M2;
          // The three local members need not be updated. They wouldn't have changed.
          break;
        }
      }
    }
  }
  for (auto i=0; i<foreign_active_comp_pattern_keys.size(); i++){
    comp_pattern_key id(foreign_active_comp_pattern_keys[i].pattern_index,foreign_active_comp_pattern_keys[i].tag,foreign_active_comp_pattern_keys[i].flops,
                        foreign_active_comp_pattern_keys[i].param1,foreign_active_comp_pattern_keys[i].param2,foreign_active_comp_pattern_keys[i].param3,
                        foreign_active_comp_pattern_keys[i].param4,foreign_active_comp_pattern_keys[i].param5); 
    int j = size_array[0]+i;
    if (comp_batch_map.find(id) != comp_batch_map.end()){
      // Here, I want to search and replace the existing batch entry with matching state
      for (auto& it : comp_batch_map[id]){
        if (it.hash_id == foreign_active_batches[j].hash_id){
          it.num_schedules = foreign_active_batches[j].num_schedules;
          it.num_scheduled_units = foreign_active_batches[j].num_scheduled_units;
          it.total_exec_time = foreign_active_batches[j].total_exec_time;
          it.M1 = foreign_active_batches[j].M1;
          it.M2 = foreign_active_batches[j].M2;
          // The three local members need not be updated. They wouldn't have changed.
          break;
        }
      }
    }
  }

  // Final update loop: iterate over the batch maps and update the batch maps.
  for (auto& it : comm_batch_map){
    for (auto& batch_state : it.second){
      // Apply the same checks on each batch's state to verify whether or not to update it's envelope
      if (batch_state.registered_channels.find(comm_channel_map[tracker.comm]) != batch_state.registered_channels.end()) continue;
      if (aggregate_channel_map.find(batch_state.hash_id ^ aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag) == aggregate_channel_map.end()) continue;

      batch_state.hash_id ^= aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag;
      batch_state.channel_count++;
      batch_state.registered_channels.insert(comm_channel_map[tracker.comm]);// Add the solo channel, not the aggregate
    }
    //assert(comm_pattern_map.find(it.first) != comm_pattern_map.end());
    if ((it.second.size() > 0) && (it.first.tag >= 13) && (it.first.tag <= 19)){
      assert(it.second.size() == 1);
      // TODO: liquidate p2p batch into its pathset
      bool is_steady = steady_test(it.first,comm_pattern_map[it.first],comm_analysis_param);
      bool is_global_steady = should_schedule(comm_pattern_map[it.first]);
      if (!is_steady && is_global_steady){
        // No need to update kernel statistics, because its already been done in the loops above
        set_kernel_state(comm_pattern_map[it.first],!is_steady);
      }
      it.second.clear();
      continue;
    }
    merge_batches(it.second,comm_analysis_param);
    bool is_steady = steady_test(it.first,comm_pattern_map[it.first],comm_analysis_param);
    bool is_global_steady = should_schedule(comm_pattern_map[it.first]);
    if (!is_steady && is_global_steady){
      // No need to update kernel statistics, because its already been done in the loops above
      set_kernel_state(comm_pattern_map[it.first],!is_steady);
/*
      .. i also want to update global kernel state here somewhere
      if ((!schedule_decision) && (comm_sample_aggregation_mode == 1)){
        //assert(comm_batch_map.find(key) != comm_batch_map.end());
        if (comm_batch_map[key].size() > 0){
          auto stats = intermediate_stats(comm_pattern_map[key],comm_batch_map[key]);
          update_kernel_stats(comm_pattern_map[key],stats);
          comm_batch_map[key].clear();
        }
      }
*/
    }
  }
  for (auto& it : comp_batch_map){
    for (auto& batch_state : it.second){
      // Apply the same checks on each batch's state to verify whether or not to update it's envelope
      if (batch_state.registered_channels.find(comm_channel_map[tracker.comm]) != batch_state.registered_channels.end()) continue;
      if (aggregate_channel_map.find(batch_state.hash_id ^ aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag) == aggregate_channel_map.end()) continue;

      batch_state.hash_id ^= aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag;
      batch_state.channel_count++;
      batch_state.registered_channels.insert(comm_channel_map[tracker.comm]);// Add the solo channel, not the aggregate
    }
    //assert(comp_pattern_map.find(it.first) != comp_pattern_map.end());
    merge_batches(it.second,comp_analysis_param);
    bool is_steady = steady_test(it.first,comp_pattern_map[it.first],comp_analysis_param);
    bool is_global_steady = should_schedule(comp_pattern_map[it.first]);
    if (!is_steady && is_global_steady){
      // No need to update kernel statistics, because its already been done in the loops above
      set_kernel_state(comp_pattern_map[it.first],!is_steady);
/*
      .. i also want to update global kernel state here somewhere
      .. make sure if 'comp_state_aggregation_mode==1', that we update the global state somewhere...
      // If steady, and sample_aggregation_mode==1, then liquidate all batch states into the pathset and clear the corresponding batch entry in map.
      //   Only need to perform this once.
      if (is_steady && (comp_sample_aggregation_mode == 1)){
        auto stats = intermediate_stats(comp_pattern_map[key],comp_batch_map[key]);
        update_kernel_stats(comp_pattern_map[key],stats);
        comp_batch_map[key].clear();
        //TODO: Why don't we update global kernel state here?
      }
*/
    }
  }
//  if (world_rank == 8) std::cout << "LEAVE\n";
}

}
}
}
