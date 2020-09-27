#include "path.h"
#include "../container/symbol_tracker.h"
#include "../../replay/path/path.h"
#include "../../discretization/util/util.h"
#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

void path::exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm){
  int world_comm_size; MPI_Comm_size(MPI_COMM_WORLD,&world_comm_size);
  int old_comm_size; MPI_Comm_size(oldcomm,&old_comm_size);
  int new_comm_size; MPI_Comm_size(newcomm,&new_comm_size);
  int world_comm_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_comm_rank);
  int new_comm_rank; MPI_Comm_rank(newcomm,&new_comm_rank);

  // There is no way to know whether each generated newcomm is of equal size without global knowledge of all colors specified (i.e. an Allgather)
  //   Therefore, I will just set the below assert inside the branch, which is motivated by the fact that I am not yet sure how to handle stride-specification for channels with a single process.
  if (world_comm_size<=1) return;
  //int old_comm_rank; MPI_Comm_rank(oldcomm,&old_comm_rank);
  solo_channel* node = new solo_channel();
  std::vector<int> gathered_info(new_comm_size,0);
  PMPI_Allgather(&world_comm_rank,1,MPI_INT,&gathered_info[0],1,MPI_INT,newcomm);
  // Now we detect the "color" (really the stride) via iteration
  // Step 1: subtract out the offset from 0 : assuming that the key arg to comm_split didn't re-shuffle
  //         I can try to use std::min_element vs. writing my own manual loop
  std::sort(gathered_info.begin(),gathered_info.end());
  node->offset = gathered_info[0];
  for (auto i=0; i<gathered_info.size(); i++) { gathered_info[i] -= node->offset; }
  node->id = channel::generate_tuple(gathered_info,new_comm_size);
  node->tag = communicator_count++;
  // Local hash_str will include offset and (size,stride) of each tuple
  // Global hash_str will include (size,stride) of each tuple -> not this may be changed (break example would be diagonal+row)
  // Note: I want to incorporate 'node->offset' into local_channel_hash_str', yet the issue is that when I go to iterate over aggregates,
  //   there is no guarantee that the global hash tags will be in same sorted order as local_hash_tag.
  std::string local_channel_hash_str = "";//std::to_string(node->offset);
  std::string global_channel_hash_str = "";
  for (auto i=0; i<node->id.size(); i++){
    local_channel_hash_str += ".." + std::to_string(node->id[i].first) + "." + std::to_string(node->id[i].second);
    global_channel_hash_str += ".." + std::to_string(node->id[i].first) + "." + std::to_string(node->id[i].second);
  }
  node->local_hash_tag = std::hash<std::string>()(local_channel_hash_str);// will avoid any local overlap.
  node->global_hash_tag = std::hash<std::string>()(global_channel_hash_str);// will avoid any global overlap.
  spf.insert_node(node);// This call will just fill in SPT via node's parent/children members, and the members of related channels
  comm_channel_map[newcomm] = node;

  // Recursively build up other legal aggregate channels that include 'node'
  std::vector<int> local_hash_array;
  std::vector<aggregate_channel*> new_aggregate_channels;
  int max_sibling_node_size=0;
  std::vector<int> save_max_indices;
  // Check if 'node' is a sibling of all existing aggregates already formed. Note that we do not include p2p aggregates, nor p2p+comm aggregates.
  // Note this loop assumes that the local_hash_tags of each aggregate across new_comm are in the same sorted order (hence the assert below)
  for (auto it : aggregate_channel_map){
    // 0. Check that each process in newcomm is processing the same aggregate.
    int verify_global_agg_hash;
    PMPI_Allreduce(&it.second->global_hash_tag,&verify_global_agg_hash,1,MPI_INT,MPI_MIN,newcomm);
    //if (verify_global_agg_hash != it.second->global_hash_tag) std::cout << "Verify - " << verify_global_agg_hash << " " << it.second->global_hash_tag << std::endl;
    assert(verify_global_agg_hash == it.second->global_hash_tag);
    // 1. Check if 'node' is a child of 'aggregate'
    bool is_child_1 = channel::verify_ancestor_relation(it.second,node);
    // 2. Check if 'aggregate' is a child of 'node'
    bool is_child_2 = channel::verify_ancestor_relation(node,it.second);
    // 3. Check if 'node'+'aggregate' form a sibling
    bool is_sibling = channel::verify_sibling_relation(it.second,node);
    if (is_sibling && !is_child_1 && !is_child_2){
      // If current aggregate forms a larger one with 'node', reset its 'is_final' member to be false, and always set a new aggregate's 'is_final' member to true
      it.second->is_final = false;
      int new_local_hash_tag = it.second->local_hash_tag ^ node->local_hash_tag;
      int new_global_hash_tag = it.second->global_hash_tag ^ node->global_hash_tag;
      auto new_aggregate_channel = new aggregate_channel(it.second->id,new_local_hash_tag,new_global_hash_tag,0,it.second->num_channels+1);// '0' gets updated below
      // Set the hashes of each communicator.
      new_aggregate_channel->channels.insert(node->local_hash_tag);
      for (auto it_2 : it.second->channels){
        new_aggregate_channel->channels.insert(it_2);
      }
      // Communicate to attain the minimum offset of all process in newcomm's aggregate channel.
      PMPI_Allgather(&it.second->offset,1,MPI_INT,&gathered_info[0],1,MPI_INT,newcomm);
      std::sort(gathered_info.begin(),gathered_info.end());
      new_aggregate_channel->offset = gathered_info[0];
      assert(new_aggregate_channel->offset <= it.second->offset);
      for (auto i=0; i<gathered_info.size(); i++) { gathered_info[i] -= new_aggregate_channel->offset; }
      auto tuple_list = channel::generate_tuple(gathered_info,new_comm_size);
      // Generate IR for new aggregate by replacing newcomm's tuple with that of the offsets of its distinct aggregates.
      for (auto it_2 : tuple_list){
        new_aggregate_channel->id.push_back(it_2);
      }
      std::sort(new_aggregate_channel->id.begin(),new_aggregate_channel->id.end(),[](const std::pair<int,int>& p1, const std::pair<int,int>& p2){return p1.second < p2.second;});
      channel::contract_tuple(new_aggregate_channel->id);
      new_aggregate_channels.push_back(new_aggregate_channel);
      local_hash_array.push_back(new_local_hash_tag);
      if (new_aggregate_channels[new_aggregate_channels.size()-1]->num_channels > max_sibling_node_size){
        max_sibling_node_size = new_aggregate_channels[new_aggregate_channels.size()-1]->num_channels;
        save_max_indices.clear();
        save_max_indices.push_back(new_aggregate_channels.size()-1);
      }
      else if (new_aggregate_channels[new_aggregate_channels.size()-1]->num_channels == max_sibling_node_size){
        save_max_indices.push_back(new_aggregate_channels.size()-1);
      }
    }
  }
  // Populate the aggregate_channel_map with the saved pointers that were created in the loop above.
  int index_window=0;
  for (auto i=0; i<new_aggregate_channels.size(); i++){
    // Update is_final to true iff its the largest subset size that includes 'node' (or if there are multiple)
    if ((index_window < save_max_indices.size()) && (save_max_indices[index_window]==i)){ 
      new_aggregate_channels[i]->is_final=true;
      index_window++;
    }
    assert(aggregate_channel_map.find(new_aggregate_channels[i]->local_hash_tag) == aggregate_channel_map.end());
    aggregate_channel_map[new_aggregate_channels[i]->local_hash_tag] = new_aggregate_channels[i];
    if (world_comm_rank == 8){
      auto str1 = channel::generate_tuple_string(new_aggregate_channels[i]);
      auto str2 = aggregate_channel::generate_hash_history(new_aggregate_channels[i]);
      std::cout << "Process " << world_comm_rank << " has aggregate " << str1 << " " << str2 << " with hashes (" << new_aggregate_channels[i]->local_hash_tag << " " << new_aggregate_channels[i]->global_hash_tag << "), num_channels - " << new_aggregate_channels[i]->num_channels << std::endl;
    }
  }

  // Verify that the aggregates are build with the same hashes
  int local_sibling_size = local_hash_array.size();
  // Always treat 1-communicator channels as trivial aggregate channels.
  aggregate_channel* agg_node = new aggregate_channel(node->id,node->local_hash_tag,node->global_hash_tag,node->offset,1);
  agg_node->channels.insert(node->local_hash_tag);
  assert(aggregate_channel_map.find(node->local_hash_tag) == aggregate_channel_map.end());
  aggregate_channel_map[node->local_hash_tag] = agg_node;
  if (world_comm_rank == 8){
    auto str1 = channel::generate_tuple_string(agg_node);
    auto str2 = aggregate_channel::generate_hash_history(agg_node);
    std::cout << "Process " << world_comm_rank << " has aggregate " << str1 << " " << str2 << " with hashes (" << agg_node->local_hash_tag << " " << agg_node->global_hash_tag << "), num_channels - " << agg_node->num_channels << std::endl;
  }
  
  if (local_sibling_size==0){// Only if 'node' exists as the smallest trivial aggregate should it be considered final. Think of 'node==world' of the very first registered channel
    aggregate_channel_map[node->local_hash_tag]->is_final=true;
  }
  PMPI_Barrier(oldcomm);
  computation_timer = MPI_Wtime();
}

bool path::initiate_comp(size_t id, volatile double curtime, double flop_count, int param1, int param2, int param3, int param4, int param5){
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
  comp_pattern_key key(-1,id,flop_count,param1,param2,param3,param4,param5);// '-1' argument is arbitrary, does not influence overloaded operators
  // Below, the idea is that key doesn't exist in comp_pattern_map iff the key hasn't been seen before. If the key has been seen, we automatically
  //   create an entry in comp_pattern_key, although it will be empty.
  if (comp_pattern_map.find(key) != comp_pattern_map.end()){
    schedule_decision = should_schedule(comp_pattern_map[key])==1;
  }

  comp_intercept_overhead += MPI_Wtime() - overhead_start_time;
  // start compunication timer for compunication routine
  comp_start_time = MPI_Wtime();
  return schedule_decision;
}

void path::complete_comp(size_t id, double flop_count, int param1, int param2, int param3, int param4, int param5){
  volatile double comp_time = MPI_Wtime( ) - comp_start_time;	// complete computation time
  // Special exit if no kernels are to be scheduled -- the goal is to track the total overhead time (no comp/comm kernels), which should
  //   be attained with timers outside of critter.
  if (schedule_kernels==0){ return; }
  volatile double overhead_start_time = MPI_Wtime();

  comp_pattern_key key(active_patterns.size(),id,flop_count,param1,param2,param3,param4,param5);// 'active_patterns.size()' argument is arbitrary, does not influence overloaded operators
  // Below, the idea is that key doesn't exist in comp_pattern_map iff the key hasn't been seen before. If the key has been seen, we automatically
  //   create an entry in comp_pattern_key, although it will be empty.
  if (comp_pattern_map.find(key) == comp_pattern_map.end()){
    active_comp_pattern_keys.push_back(key);
    active_patterns.emplace_back();
    comp_pattern_map[key] = pattern_key_id(true,active_comp_pattern_keys.size()-1,active_patterns.size()-1,false);
  }
  if ((aggregation_mode==0) || (should_schedule(comp_pattern_map[key]) == 0)){
    // Because of the second case in the branch, this will be called even if in aggregation_mode>=1 as long as the kernel is no longer being scheduled (i.e. is in steady state)
    update_kernel_stats(comp_pattern_map[key],comp_analysis_param,comp_time,flop_count);
  }
  else if (aggregation_mode >= 1){
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
    assert(comp_batch_map[key].size() == 0);
    comp_time = get_estimate(comp_pattern_map[key],comp_analysis_param,flop_count);
  } else{
    // Both non-optimized and optimized variants can update the local kernel state and the global kernel state (no chance of deadlock)
    bool is_steady = steady_test(key,comp_pattern_map[key],comp_analysis_param);
    set_kernel_state(comp_pattern_map[key],!is_steady);
    set_kernel_state_global(comp_pattern_map[key],!is_steady);
    // If steady, and aggregation_mode>=1, then liquidate all batch states into the pathset and clear the corresponding batch entry in map.
    //   Only need to perform this once. Note that this is allowed, whereas it is not allowed in 'complete_comm', because
    //   a processor's kernel may be steady, yet its not globally steady.
    if (is_steady && (aggregation_mode >= 1)){
      auto stats = intermediate_stats(comp_pattern_map[key],comp_batch_map[key]);
      update_kernel_stats(comp_pattern_map[key],stats);
      comp_batch_map[key].clear();
    }
    if (is_steady) { flush_pattern(key); }
  }

  critical_path_costs[num_critical_path_measures-1] += comp_time;	// execution time
  critical_path_costs[num_critical_path_measures-2] += comp_time;	// computation kernel time
  critical_path_costs[num_critical_path_measures-3] += comp_time;	// computational time
  volume_costs[num_volume_measures-1] += comp_time;			// execution time
  volume_costs[num_volume_measures-2] += comp_time;			// computation kernel time
  volume_costs[num_volume_measures-3] += comp_time;			// computation time

  comp_intercept_overhead += MPI_Wtime() - overhead_start_time;
  computation_timer = MPI_Wtime();
}

bool path::initiate_comm(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
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

  // Non-optimized variant will always post barriers, although of course, just as with the optimized variant, the barriers only remove idle time
  //   from corrupting communication time measurements. The process that enters barrier last is not necessarily the critical path root. The
  //     critical path root is decided based on a reduction using 'critical_path_costs'. Therefore, no explicit barriers are invoked, instead relying on Allreduce
  bool post_barrier = true; bool schedule_decision = true;
  double reduced_info[5] = {critical_path_costs[num_critical_path_measures-4],critical_path_costs[num_critical_path_measures-3],
                            critical_path_costs[num_critical_path_measures-2],critical_path_costs[num_critical_path_measures-1],0};
  double reduced_info_foreign[5] = {0,0,0,0,0};

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
  auto world_key = key; int world_rank;
  if (tracker.partner1 != -1){
    world_rank = channel::translate_rank(tracker.comm,tracker.partner1); 
    world_key.partner_offset = tracker.partner1;//world_rank;
/*
    if (p2p_global_state_override.find(world_key) == p2p_global_state_override.end()){
      p2p_global_state_override[world_key] = true;// not in global steady state
    }
*/
  }
  if (comm_pattern_map.find(key) != comm_pattern_map.end()){
    if (is_optimized){
      post_barrier = should_schedule_global(comm_pattern_map[key])==1;
      if (tracker.partner1 != -1){
        //post_barrier = post_barrier || p2p_global_state_override[world_key];
      }
    }
    schedule_decision = should_schedule(comm_pattern_map[key])==1;
    reduced_info[4] = (double)schedule_decision;
  }
  // Register the p2p channel
  if (aggregation_mode >= 1){
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
        spf.insert_node(node);
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
  }
      
  // Both unoptimized and optimized variants will reduce an execution time and a schedule_decision if post_barrier==true.
  // Note that formally, the unoptimized variant does not need to reduce a schedule_decision, but I allow it becase the
  //   procedure is synchronization bound anyway and adding an extra double won't matter.
  // The unoptimized variant will always perform this reduction of execution time, regardless of its schedule_decision.
  // The optimized variant will only perform this reduction until it detects global steady state.

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

    // Below, the idea is that key doesn't exist in comm_pattern_map iff the key hasn't been seen before. If the key has been seen, we automatically
    //   create an entry in comm_pattern_key, although it will be empty.
    if (comm_pattern_map.find(key) != comm_pattern_map.end()){
      // If local steady_state has been reached, and we find out the other processes have reached the same, then we can set global_steady_state=1
      //   and never have to use another internal collective to check again.
      if (is_optimized){
        // Do not update kernel state after its been set to steady. This situation can occur with use of 'p2p_global_state_override' to prevent deadlock
        // That also means we won't update our statistics with the scheduled comm_time, for better or for worse.
        //if (should_schedule_global(comm_pattern_map[key])==1){
          set_kernel_state(comm_pattern_map[key],schedule_decision);
          set_kernel_state_global(comm_pattern_map[key],schedule_decision);
        //}
        if (tracker.partner1 != -1){
          //TODO: Why is this not working???
          //p2p_global_state_override[world_key] = schedule_decision;
        }
        // TODO: The statement below worries me in the presence of 'p2p_global_state_override'
        if (!schedule_decision) { flush_pattern(key); }
      } else{
        // Not sure if this is needed. It will keep setting to steady if all processors in tracker.comm have a steady kernel, yet never
        //   in global steady state (i.e. we will always need a reduction to make sure).
        set_kernel_state(comm_pattern_map[key],schedule_decision);
      }
      // Note: once a kernel is steady, we want to liquidate batch information into the static pathset. Only need to do this once
      if ((!schedule_decision) && (aggregation_mode >= 1)){
        assert(comm_batch_map.find(key) != comm_batch_map.end());
        if (comm_batch_map[key].size() > 0){
          auto stats = intermediate_stats(comm_pattern_map[key],comm_batch_map[key]);
          update_kernel_stats(comm_pattern_map[key],stats);
          comm_batch_map[key].clear();
        }
      }
    }
    // If kernel is about to be scheduled, post one more barrier for safety if collective,
    //   because the AllReduce posted above may allow ranks to leave early, thus corrupting the sample measurement.
    if (schedule_decision && tracker.partner1 == -1){ PMPI_Barrier(tracker.comm); }
    //TODO: Might consider a sendrecv here if tracker.partner1 != -1
  }
  else{
    if (comm_pattern_map.find(key) != comm_pattern_map.end()){
      if (comm_batch_map.find(key) != comm_batch_map.end()){
        assert(comm_batch_map[key].size()==0);
      }
    }
  }

  comm_intercept_overhead_stage1 += MPI_Wtime() - overhead_start_time;
  // start communication timer for communication routine
  tracker.start_time = MPI_Wtime();
  return schedule_decision;
}

// Used only for p2p communication. All blocking collectives use sychronous protocol
void path::complete_comm(blocking& tracker, int recv_source){
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

  bool should_propagate = true;
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
    active_patterns.emplace_back();
    comm_pattern_map[key] = pattern_key_id(true,active_comm_pattern_keys.size()-1,active_patterns.size()-1,false);
  }
  int comm_hash_tag;
  if (should_schedule_global(comm_pattern_map[key]) == 0) assert(should_schedule(comm_pattern_map[key])==0);
  if ((aggregation_mode==0) || (should_schedule(comm_pattern_map[key]) == 0)){
    update_kernel_stats(comm_pattern_map[key],comm_analysis_param,comm_time,tracker.nbytes);
  }
  else if (aggregation_mode >= 1){
    assert(should_schedule_global(comm_pattern_map[key]) == 1);
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
    assert(comm_batch_map[key].size() == 0);
    if (is_optimized==1){ should_propagate = false; }
  } else{
    assert(should_schedule_global(comm_pattern_map[key]) == 1);
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
      if (aggregation_mode==0){
        // Note: in optimized version, exchange_patterns and flush_patterns would never need to be called
        if (is_optimized==0){
          if (is_world_communication && !is_key_skipable(key)){//Note: for practical reasons, we force more constraints on when profile exchanges take place
            exchange_patterns_per_process(tracker);
            flush_patterns();
          }
        }
      }
      else if (aggregation_mode >= 1){
        // TODO: Note: I expect to update each kernel's state (not global state though, and only if kernel's state is active, as specified in the complete_pathset)
        //         in comm_pattern_map and comp_pattern_map. Remember, the state specified in the complete_pathset does not entirely reflect the data set in the complete_pathset,
        //           as it takes into account data still in active batches.
        // TODO: Note: post-propagation, check if any batch has an open_count==0. If so, add it to the complete pathset.
        //         Then figure out whether to update its values to indicate no samples, or remove it altogether. After thinking about it,
        //           its probably safest to swap with the last pattern_batch, and then pop_back.
        if (tracker.partner1 == -1){
          if (aggregation_mode==1){
            single_stage_sample_aggregation(tracker);
          } else if (aggregation_mode==2){
            multi_stage_sample_aggregation(tracker);
          } else { assert(0); }
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

void path::flush_pattern(comp_pattern_key key){
  return;
}

void path::flush_pattern(comm_pattern_key key){
  return;
  // 'flush_pattern' is used for the sole purpose of keeping the invariant that an inactive kernel's statistical data
  //    belongs in the steady_state buffers. This procedure need only occur in 'flush_patterns' for critical path analysis,
  //      but as in per-process/volumetric analysis the global state can change in 'initiate_comm', this routine serves as a way to
  //        achieve that.
  assert(comm_pattern_map[key].is_active);
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
    comm_pattern_key id(active_comm_pattern_keys[i].pattern_index,active_comm_pattern_keys[i].tag,active_comm_pattern_keys[i].dim_sizes,
                        active_comm_pattern_keys[i].dim_strides,active_comm_pattern_keys[i].msg_size,active_comm_pattern_keys[i].partner_offset); 
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


// No overload on decltype(tracker) because this exchange occurs only with global communication
// Note: this exchange does not need to exchange computation kernels in the active pathset, because
//       it does not prevent feasibility of the idea. Deadlock is not possible.
// Note: entire data structures, rather than just the is_active member, are required to be exchanged,
//       as otherwise, it wouldn't be feasible: receiver would have no way of verifying whether the sent
//       order matches its own order.
// Note: per-process statistical analysis of kernels occurs only within a phase.
//       at phase-end, we exchange to pdate kernel state in a global way to prevent deadlock in scheduling comm patterns.
void path::exchange_patterns_per_process(blocking& tracker){
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
        comm_pattern_key id(foreign_active_comm_pattern_keys[i].pattern_index,foreign_active_comm_pattern_keys[i].tag,foreign_active_comm_pattern_keys[i].dim_sizes,
                            foreign_active_comm_pattern_keys[i].dim_strides,foreign_active_comm_pattern_keys[i].msg_size,foreign_active_comm_pattern_keys[i].partner_offset); 
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
    comm_pattern_key id(foreign_active_comm_pattern_keys[i].pattern_index,foreign_active_comm_pattern_keys[i].tag,foreign_active_comm_pattern_keys[i].dim_sizes,
                        foreign_active_comm_pattern_keys[i].dim_strides,foreign_active_comm_pattern_keys[i].msg_size,foreign_active_comm_pattern_keys[i].partner_offset); 
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

void path::single_stage_sample_aggregation(blocking& tracker){
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  int comm_rank; MPI_Comm_rank(tracker.comm,&comm_rank);
  int comm_size; MPI_Comm_size(tracker.comm,&comm_size);
  std::vector<comm_pattern_key> foreign_active_comm_pattern_keys;
  std::vector<comp_pattern_key> foreign_active_comp_pattern_keys;
  std::vector<pattern_batch_propagate> foreign_active_batches;
/*
  if (world_rank == 0) std::cout << "CHANNEL - " << comm_channel_map[tracker.comm]->global_hash_tag << std::endl;

  // Debug
      for (auto& it : comm_batch_map){
        assert(it.second.size() < 2);// Number of states/batches of a particular key must be < 2. This is almost trivially satisfied.
        if (it.second.size() == 1){// Some comm keys will have no active batches/states, so this check is mandatory
          if (world_rank==0){
             std::cout << "BEFORE -- What is this comm " << " (" << it.first.tag << "," << it.first.dim_sizes[0] << ","
                       << it.first.dim_sizes[1] << "," << it.first.dim_strides[0] << "," << it.first.dim_strides[1] << "," << it.first.msg_size << "," << it.first.partner_offset << ") - " << comm_batch_map[it.first][0].num_local_schedules << " " << comm_batch_map[it.first][0].num_schedules << std::endl;
          }
        }
      }
      for (auto& it : comp_batch_map){
        assert(it.second.size() < 2);
        if (it.second.size() == 1){
          if (world_rank==0){
            std::cout << "BEFORE -- What is this comp " << " (" << it.first.tag << "," << it.first.flops << ","
                      << it.first.param1 << "," << it.first.param2 << "," << it.first.param3 << ") - " << comp_batch_map[it.first][0].num_local_schedules << " " << comp_batch_map[it.first][0].num_schedules << std::endl;
          }
        }
      }
*/

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
/*
          if (world_rank==0){
             std::cout << "CONTRIBUTED -- What is this comm " << " (" << it.first.tag << "," << it.first.dim_sizes[0] << ","
                       << it.first.dim_sizes[1] << "," << it.first.dim_strides[0] << "," << it.first.dim_strides[1] << "," << it.first.msg_size << "," << it.first.partner_offset << ") - " << comm_batch_map[it.first][0].num_local_schedules << " " << comm_batch_map[it.first][0].num_schedules << std::endl;
          }
*/
          foreign_active_comm_pattern_keys.push_back(it.first);
          foreign_active_batches.emplace_back(it.second[0]);
        }
      }
      for (auto& it : comp_batch_map){
        assert(it.second.size() < 2);
        if (it.second.size() == 1){
/*
          if (world_rank==0){
            std::cout << "CONTRIBUTED -- What is this comp " << " (" << it.first.tag << "," << it.first.flops << ","
                      << it.first.param1 << "," << it.first.param2 << "," << it.first.param3 << ") - " << comp_batch_map[it.first][0].num_local_schedules << " " << comp_batch_map[it.first][0].num_schedules << std::endl;
          }
*/
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
/*
        if (world_rank==0){
           std::cout << "CONTRIBUTED -- What is this comm " << " (" << it.first.tag << "," << it.first.dim_sizes[0] << ","
                     << it.first.dim_sizes[1] << "," << it.first.dim_strides[0] << "," << it.first.dim_strides[1] << "," << it.first.msg_size << "," << it.first.partner_offset << ") - " << comm_batch_map[it.first][0].num_local_schedules << " " << comm_batch_map[it.first][0].num_schedules << std::endl;
        }
*/
        foreign_active_comm_pattern_keys.push_back(it.first);
        foreign_active_batches.push_back(it.second[0]);
      }
    }
    for (auto& it : comp_batch_map){
      assert(it.second.size() < 2);
      if (it.second.size() == 1){
/*
        if (world_rank==0){
          std::cout << "CONTRIBUTED -- What is this comp " << " (" << it.first.tag << "," << it.first.flops << ","
                    << it.first.param1 << "," << it.first.param2 << "," << it.first.param3 << ") - " << comp_batch_map[it.first][0].num_local_schedules << " " << comp_batch_map[it.first][0].num_schedules << std::endl;
        }
*/
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
/*
        if (world_rank==0){
          std::cout << "AFTER -- What is this comm - " << i << " " << comm_pattern_map[id].val_index << " (" << id.tag << "," << id.dim_sizes[0] << ","
                    << id.dim_sizes[1] << "," << id.dim_strides[0] << "," << id.dim_strides[1] << "," << id.msg_size << "," << id.partner_offset << ") - "
                    << comm_batch_map[id][0].num_local_schedules << " " << comm_batch_map[id][0].num_schedules << std::endl;
        }
*/
        temp.num_local_schedules = comm_batch_map[id][0].num_local_schedules;
        temp.num_local_scheduled_units = comm_batch_map[id][0].num_local_scheduled_units;
        temp.total_local_exec_time = comm_batch_map[id][0].total_local_exec_time;
      }
      comm_batch_map[id].clear();// the samples residing here have been propagated and the batch's journey is complete.
    }
    // Update existing entry.
    bool is_steady = steady_test(id,comm_pattern_map[id],comm_analysis_param);
    bool is_global_steady = should_schedule_global(comm_pattern_map[id]);
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
/*
        if (world_rank==0){
          std::cout << "AFTER -- What is this comp - " << i << " " << comp_pattern_map[id].val_index << " (" << id.tag << "," << id.flops << ","
                    << id.param1 << "," << id.param2 << "," << id.param3 << ") - "
                    << comp_batch_map[id][0].num_local_schedules << " " << comp_batch_map[id][0].num_schedules << std::endl;
        }
*/
        temp.num_local_schedules = comp_batch_map[id][0].num_local_schedules;
        temp.num_local_scheduled_units = comp_batch_map[id][0].num_local_scheduled_units;
        temp.total_local_exec_time = comp_batch_map[id][0].total_local_exec_time;
      }
      comp_batch_map[id].clear();// the samples residing here have been propagated and the batch's journey is complete.
    }
    // Update existing entry.
    bool is_steady = steady_test(id,comp_pattern_map[id],comp_analysis_param);
    bool is_global_steady = should_schedule_global(comp_pattern_map[id]);
    if (!is_steady && is_global_steady){
      update_kernel_stats(active_patterns[comp_pattern_map[id].val_index], temp, comp_analysis_param);
      bool is_steady = steady_test(id,comp_pattern_map[id],comp_analysis_param);
      set_kernel_state(comp_pattern_map[id],!is_steady);
    }
  }
//  if (world_rank == 0) std::cout << "LEAVE\n";
}
void path::multi_stage_sample_aggregation(blocking& tracker){
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  //if (world_rank == 8) std::cout << "\n\nChannel - " << comm_channel_map[tracker.comm]->id[0].first << " " << comm_channel_map[tracker.comm]->id[0].second << " " << comm_channel_map[tracker.comm]->global_hash_tag << " " << comm_channel_map[tracker.comm] << " with msg size " << tracker.nbytes << std::endl;
  int comm_rank; MPI_Comm_rank(tracker.comm,&comm_rank);
  int comm_size; MPI_Comm_size(tracker.comm,&comm_size);
  std::vector<comm_pattern_key> foreign_active_comm_pattern_keys;
  std::vector<comp_pattern_key> foreign_active_comp_pattern_keys;
  std::vector<pattern_batch_propagate> foreign_active_batches;

  std::map<comm_pattern_key,std::vector<pattern_batch>> temp_comm_batch_map;
  std::map<comp_pattern_key,std::vector<pattern_batch>> temp_comp_batch_map;
  // Fill up temporary maps with those local batches fit to aggregate across this channel 'tracker.comm'
  for (auto& it : comm_batch_map){
//    if (world_rank == 8) std::cout << "comm key " << it.first.tag << " " << it.first.dim_sizes[0] << " " << it.first.dim_strides[0] << " " << it.first.msg_size << " " << it.first.partner_offset << ") has " << it.second.size() << " active batches and " << active_patterns[comm_pattern_map[it.first].val_index].num_schedules << " total samples and " << active_patterns[comm_pattern_map[it.first].val_index].num_local_schedules << " local samples\n";
    for (auto& batch_state : it.second){
//      if (world_rank == 8) std::cout << "\tWith state " << batch_state.hash_id << ", num_schedules " << batch_state.num_schedules << ", num_local_schedules - " << batch_state.num_local_schedules << ", M1 - " << batch_state.M1 << ", M2 - " << batch_state.M2 << " " << get_error_estimate(it.first,comm_pattern_map[it.first],comm_analysis_param) << std::endl;
      // Three checks below:
      //   1. Does this propagation channel match the channel of the particular batch, or has it been used before?
      //   2. Does this propagation channel form an aggregate with the channel of the particular batch?
      //   3. Does this batch have a stride of 0? This indicates a trivial p2p that does not need sending
      if (batch_state.registered_channels.find(comm_channel_map[tracker.comm]) != batch_state.registered_channels.end()) continue;
      // TODO: Not exactly sure whether to use global_hash_id below or local_hash_id
      if (aggregate_channel_map.find(batch_state.hash_id ^ aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag) == aggregate_channel_map.end()) continue;
//      if ((it.second[0].hash_id == 1) && (it.second[0].id[0].second == 0)) continue;
      temp_comm_batch_map[it.first].push_back(batch_state);
/*
      if (world_rank == 8){
        std::cout << "\t\tContributed - " << batch_state.channel_count << " " << batch_state.registered_channels.size();
        for (auto& blah : batch_state.registered_channels){ std::cout << " (" << blah->global_hash_tag << "," << blah << ")"; }
        std::cout << "\n";
      }
*/
    }
  }
  for (auto& it : comp_batch_map){
//    if (world_rank == 8) std::cout << "comp key " << it.first.tag << " " << it.first.flops << " " << it.first.param1 << " " << it.first.param2 << " " << it.first.param3 << ") has " << it.second.size() << " active batches and " << active_patterns[comp_pattern_map[it.first].val_index].num_schedules << " total samples\n";
    for (auto& batch_state : it.second){
//      if (world_rank == 8) std::cout << "\tWith state " << batch_state.hash_id << ", num_schedules " << batch_state.num_schedules << ", num_local_schedules - " << batch_state.num_local_schedules << ", M1 - " << batch_state.M1 << ", M2 - " << batch_state.M2 << " " << get_error_estimate(it.first,comp_pattern_map[it.first],comp_analysis_param) << std::endl;
      // Three checks below:
      //   1. Has this propagation channel been used before?
      //   2. Does this propagation channel form an aggregate with the channel of the particular batch?
      if (batch_state.registered_channels.find(comm_channel_map[tracker.comm]) != batch_state.registered_channels.end()) continue;
      // TODO: Not exactly sure whether to use global_hash_id below or local_hash_id
      if (aggregate_channel_map.find(batch_state.hash_id ^ aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag) == aggregate_channel_map.end()) continue;
      temp_comp_batch_map[it.first].push_back(batch_state);
/*
      if (world_rank == 8){
        std::cout << "\t\tContributed - " << batch_state.channel_count << " " << batch_state.registered_channels.size();
        for (auto& blah : batch_state.registered_channels){ std::cout << " (" << blah->global_hash_tag << "," << blah << ")"; }
        std::cout << "\n";
      }
*/
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
      bool is_global_steady = should_schedule_global(comm_pattern_map[it.first]);
      if (!is_steady && is_global_steady){
        // No need to update kernel statistics, because its already been done in the loops above
        set_kernel_state(comm_pattern_map[it.first],!is_steady);
      }
      it.second.clear();
      continue;
    }
    merge_batches(it.second,comm_analysis_param);
    bool is_steady = steady_test(it.first,comm_pattern_map[it.first],comm_analysis_param);
    bool is_global_steady = should_schedule_global(comm_pattern_map[it.first]);
    if (!is_steady && is_global_steady){
      // No need to update kernel statistics, because its already been done in the loops above
      set_kernel_state(comm_pattern_map[it.first],!is_steady);
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
    bool is_global_steady = should_schedule_global(comp_pattern_map[it.first]);
    if (!is_steady && is_global_steady){
      // No need to update kernel statistics, because its already been done in the loops above
      set_kernel_state(comp_pattern_map[it.first],!is_steady);
    }
  }
//  if (world_rank == 8) std::cout << "LEAVE\n";
}

}
}
}
