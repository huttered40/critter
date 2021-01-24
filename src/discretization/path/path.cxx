#include "path.h"
#include "../container/symbol_tracker.h"
#include "../util/util.h"
#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

static void update_frequency(float* in, float* inout, size_t len){
  assert(len == cp_costs_size);	// this assert prevents user from obtaining wrong output if MPI implementation cuts up the message.
  if (in[num_cp_measures-1] > inout[num_cp_measures-1]){
    std::memcpy(inout+18,in+18,(len-18)*sizeof(float));
  }
  for (auto i=0; i<18; i++) inout[i] = std::max(in[i],inout[i]);
}
static void propagate_cp_op(float* in, float* inout, int* len, MPI_Datatype* dtype){
  update_frequency(in,inout,static_cast<size_t>(*len));
}

static void kernel_update(float* read_ptr){
  // Leave the un-updated kernels alone, just use the per-process count
  // Lets have all processes update, even the root, so that they leave this routine (and subsequently leave the interception) at approximately the same time.
  for (auto i=0; i<comp_kernel_select_count; i++){
    auto offset = i*9;
    if (read_ptr[offset] == -1) break;
    comp_kernel_key id(-1,(int)read_ptr[offset],read_ptr[offset+7],
                          (int)read_ptr[offset+1],(int)read_ptr[offset+2],(int)read_ptr[offset+3],
                          (int)read_ptr[offset+4],(int)read_ptr[offset+5]);
    // Don't bother adding new kernels unseen by the current processor.
    // Kernels un-updated will use per-process count as an approximation
    if (comp_kernel_map.find(id) != comp_kernel_map.end()){
      active_kernels[comp_kernel_map[id].val_index].num_local_schedules = read_ptr[offset+8];
    }
  }
  for (auto i=0; i<comm_kernel_select_count; i++){
    auto offset = comp_kernel_select_count*9+i*9;
    if (read_ptr[offset] == -1) break;
    comm_kernel_key id(-1,(int)read_ptr[offset],(int*)&read_ptr[offset+1],
                       (int*)&read_ptr[offset+3],read_ptr[offset+7],
                       (int)read_ptr[offset+5]); 
    // Don't bother adding new kernels unseen by the current processor.
    // Kernels un-updated will use per-process count as an approximation
    if (comm_kernel_map.find(id) != comm_kernel_map.end()){
      active_kernels[comm_kernel_map[id].val_index].num_local_schedules = read_ptr[offset+8];
    }
  }
}

void path::exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm){
  // Save and accumulate the computation time between last communication routine as both execution-time and computation time
  //   into both the execution-time critical path data structures and the per-process data structures.
  auto save_comp_time = MPI_Wtime() - computation_timer;
  if (mode==1){
    cp_costs[num_cp_measures-1] += save_comp_time;
    cp_costs[num_cp_measures-3] += save_comp_time;
    vol_costs[num_vol_measures-1] += save_comp_time;
    vol_costs[num_vol_measures-3] += save_comp_time;
  }

  generate_aggregate_channels(oldcomm,newcomm);
  PMPI_Barrier(oldcomm);
  if (mode==1){
    computation_timer = MPI_Wtime();
  }
}

bool path::initiate_comp(size_t id, volatile double curtime, float flop_count, int param1, int param2, int param3, int param4, int param5){
  // Save and accumulate the computation time between last communication routine as both execution-time and computation time
  //   into both the execution-time critical path data structures and the per-process data structures.
  auto save_comp_time = curtime - computation_timer;
  cp_costs[num_cp_measures-1] += save_comp_time;
  cp_costs[num_cp_measures-3] += save_comp_time;
  vol_costs[num_vol_measures-1] += save_comp_time;
  vol_costs[num_vol_measures-3] += save_comp_time;
  // Special exit if no kernels are to be scheduled -- the goal is to track the total overhead time (no comp/comm kernels), which should
  //   be attained with timers outside of critter.
  if (schedule_kernels==0){ return false; }
  volatile auto overhead_start_time = MPI_Wtime();

  bool schedule_decision = true;
  if (!(tuning_delta==0 || tuning_delta==2)){// these tuning_deltas (0, 2) signify unconditional scheduling of computation kernels
    comp_kernel_key key(-1,id,flop_count,param1,param2,param3,param4,param5);// '-1' argument is arbitrary, does not influence overloaded operators
    // Below, the idea is that key doesn't exist in comp_kernel_map iff the key hasn't been seen before. If the key has been seen, we automatically
    //   create an entry in comp_kernel_key, although it will be empty.
    if (comp_kernel_map.find(key) != comp_kernel_map.end()){
      schedule_decision = should_schedule(comp_kernel_map[key])==1;
    }
  }
  intercept_overhead[0] += MPI_Wtime() - overhead_start_time;

  // start compunication timer for computation routine
  global_schedule_decision = schedule_decision;// 'global' as in global variable
  comp_start_time = MPI_Wtime();
  return schedule_decision;
}

void path::complete_comp(double errtime, size_t id, float flop_count, int param1, int param2, int param3, int param4, int param5){
  volatile auto comp_time = MPI_Wtime() - comp_start_time - errtime;	// complete computation time

  // Special exit if no kernels are to be scheduled -- the goal is to track the total overhead time (no comp/comm kernels), which should
  //   be attained with timers outside of critter.
  if (schedule_kernels==0){ return; }
  volatile auto overhead_start_time = MPI_Wtime();

  if (!(tuning_delta==0 || tuning_delta==2)){// these deltas (0, 2) signify unconditional scheduling of computation kernels
    comp_kernel_key key(active_kernels.size(),id,flop_count,param1,param2,param3,param4,param5);// 'active_kernels.size()' argument is arbitrary, does not influence overloaded operators
    // Below, the idea is that key doesn't exist in comp_kernel_map iff the key hasn't been seen before. If the key has been seen, we automatically
    //   create an entry in comp_kernel_key, although it will be empty.
    if (comp_kernel_map.find(key) == comp_kernel_map.end()){
      active_comp_kernel_keys.push_back(key);
      active_kernels.emplace_back();
      comp_kernel_map[key] = kernel_key_id(true,active_comp_kernel_keys.size()-1,active_kernels.size()-1,false);
    } else{
      if (global_schedule_decision) comp_kernel_map[key].is_active = true;
    }
    // For debugging sanity, if the autotuner shut off execution of this kernel, then have debugger use only existing mean
    if (autotuning_debug == 1 && comp_kernel_save_map.find(key) != comp_kernel_save_map.end() && should_schedule(comp_kernel_save_map[key]) == 0){
      comp_time = get_estimate(comp_kernel_map[key],comp_analysis_param,flop_count);
    } else{
      update_kernel_stats(comp_kernel_map[key],comp_analysis_param,comp_time,flop_count);
      // Note: 'get_estimate' must be called before setting the updated kernel state. If kernel was not scheduled, comp_time set below overwrites 'comp_time'
      if (should_schedule(comp_kernel_map[key]) == 0){
        comp_time = get_estimate(comp_kernel_map[key],comp_analysis_param,flop_count);
      } else{
        bool _is_steady = steady_test(key,comp_kernel_map[key],comp_analysis_param);
        set_kernel_state(comp_kernel_map[key],!_is_steady);
        if (comp_state_aggregation_mode==0 || sample_constraint_mode==-1) set_kernel_state_global(comp_kernel_map[key],!_is_steady);
      }
    }
  }

  cp_costs[num_cp_measures-1] += comp_time;
  cp_costs[num_cp_measures-2] += comp_time;
  cp_costs[num_cp_measures-3] += comp_time;
  vol_costs[num_vol_measures-1] += comp_time;
  vol_costs[num_vol_measures-2] += comp_time;
  vol_costs[num_vol_measures-3] += comp_time;

  intercept_overhead[0] += MPI_Wtime() - overhead_start_time;
  computation_timer = MPI_Wtime();
}

bool path::initiate_comm(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                         bool is_sender, int partner1, int partner2){
  assert(partner1 != MPI_ANY_SOURCE); if ((tracker.tag == 13) || (tracker.tag == 14)){ assert(partner2 != MPI_ANY_SOURCE); }
  // Save and accumulate the computation time between last communication routine as both execution-time and computation time
  //   into both the execution-time critical path data structures and the per-process data structures.
  tracker.comp_time = curtime - computation_timer;
  volatile auto overhead_start_time = MPI_Wtime();
  cp_costs[num_cp_measures-1] += tracker.comp_time;
  cp_costs[num_cp_measures-3] += tracker.comp_time;
  vol_costs[num_vol_measures-1] += tracker.comp_time;
  vol_costs[num_vol_measures-3] += tracker.comp_time;
  // Special exit if no kernels are to be scheduled -- the goal is to track the total overhead time (no comp/comm kernels), which should
  //   be attained with timers outside of critter.
  if (schedule_kernels==0){ return false; }

  // At this point, 'cp_costs' tells me the process's time up until now. A barrier won't suffice when kernels are conditionally scheduled.
  int rank; MPI_Comm_rank(comm, &rank);

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
  tracker.should_propagate = false;
  tracker.aggregate_comp_kernels=false;
  tracker.aggregate_comm_kernels=false;

  // The process that enters barrier last is not necessarily the critical path root. The
  //   critical path root is decided based on a reduction using 'cp_costs'. Therefore, no explicit barriers are invoked, instead relying on Allreduce
  bool schedule_decision = true;
  // Assume that the communicator of either collective/p2p is registered via comm_split, and that its described using a max of 3 dimension tuples.
  assert(comm_channel_map.find(tracker.comm) != comm_channel_map.end());
  int comm_sizes[2]={0,0}; int comm_strides[2]={0,0};
  for (auto i=0; i<comm_channel_map[tracker.comm]->id.size(); i++){
    comm_sizes[i]=comm_channel_map[tracker.comm]->id[i].first;
    comm_strides[i]=comm_channel_map[tracker.comm]->id[i].second;
  }
  // Below, the idea is that key doesn't exist in comm_kernel_map iff the key hasn't been seen before. If the key has been seen, we automatically
  //   create an entry in comm_kernel_key, although it will be empty.
  comm_kernel_key key(rank,-1,tracker.tag,comm_sizes,comm_strides,tracker.nbytes,tracker.partner1);
  memset(&cp_costs[4],0,14*sizeof(float)); cp_costs[6] = 1;
  if (tuning_delta > 1){// tuning_delta of 1 signifies that communication kernel scheduling is to be unconditional
    if (comm_kernel_map.find(key) != comm_kernel_map.end()){
      schedule_decision = should_schedule(comm_kernel_map[key])==1;
      cp_costs[5] = (!schedule_decision ? 1 : 0);// set to 1 if kernel is globally steady.
      cp_costs[6] = (schedule_decision ? 1 : 0);// set to 0 if kernel is globally steady.
      if (!schedule_decision){
        // If this particular kernel is globally steady, meaning it has exhausted its state aggregation channels,
        //   then we can overwrite the '-1' with the sample mean of the globally-steady kernel
        cp_costs[4] = active_kernels[comm_kernel_map[key].val_index].steady_state;
        cp_costs[9] = active_kernels[comm_kernel_map[key].val_index].hash_id;
        cp_costs[10] = active_kernels[comm_kernel_map[key].val_index].num_schedules;
        cp_costs[11] = active_kernels[comm_kernel_map[key].val_index].num_local_schedules;
        cp_costs[12] = active_kernels[comm_kernel_map[key].val_index].num_scheduled_units;
        cp_costs[13] = active_kernels[comm_kernel_map[key].val_index].num_local_scheduled_units;
        cp_costs[14] = active_kernels[comm_kernel_map[key].val_index].M1;
        cp_costs[15] = active_kernels[comm_kernel_map[key].val_index].M2;
        cp_costs[16] = active_kernels[comm_kernel_map[key].val_index].total_exec_time;
        cp_costs[17] = active_kernels[comm_kernel_map[key].val_index].total_local_exec_time;
      }
    }
  }

  // Use pathsets, not batches, to check if kernel can leverage an aggregation. Such a kernel must be locally steady (i.e.
  //   from its own schedules, its steady), and must be able to aggregate across the channel associated with 'tracker.comm'
  if ((tracker.partner1 == -1) && (comp_state_aggregation_mode>0)){
    for (auto& it : comp_kernel_map){
      if (!((active_kernels[it.second.val_index].steady_state==1) && (should_schedule(it.second)==1))) continue;
      // Any global communicator can fast-track a communication kernel to being in global steady state. No need to match up hashes (that would only be necessary for sample aggregation)
      if (aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->is_final){
        tracker.save_comp_key.push_back(it.first);
        cp_costs[7]++;
      }
      else{
        if (active_kernels[it.second.val_index].registered_channels.find(comm_channel_map[tracker.comm]) != active_kernels[it.second.val_index].registered_channels.end()) continue;
        // TODO: Not exactly sure whether to use global_hash_id below or local_hash_id
        if (aggregate_channel_map.find(active_kernels[it.second.val_index].hash_id ^ aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag) == aggregate_channel_map.end()) continue;
        tracker.save_comp_key.push_back(it.first);
        cp_costs[7]++;
      }
    }
  }
  if ((tracker.partner1 == -1) && (comm_state_aggregation_mode>0)){
    for (auto& it : comm_kernel_map){
      if (!((active_kernels[it.second.val_index].steady_state==1) && (should_schedule(it.second)==1))) continue;
      if (aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->is_final){
        tracker.save_comm_key.push_back(it.first);
        cp_costs[8]++;
      }
      else{
        if (active_kernels[it.second.val_index].registered_channels.find(comm_channel_map[tracker.comm]) != active_kernels[it.second.val_index].registered_channels.end()) continue;
        // TODO: Not exactly sure whether to use global_hash_id below or local_hash_id
        if (aggregate_channel_map.find(active_kernels[it.second.val_index].hash_id ^ aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag) == aggregate_channel_map.end()) continue;
        tracker.save_comm_key.push_back(it.first);
        cp_costs[8]++;
      }
    }
  }

  if (sample_constraint_mode == 2){
    // Fill in -1 first because the number of distinct kernels might be less than 'comm_kernel_select_count',
    //   just to avoid confusion. A -1 tag clearly means that the kernel is void
    memset(&cp_costs[18],-1,sizeof(float)*(cp_costs.size()-18));
    // Iterate over first 'comp_kernel_select_count' keys
    int count=0;
    for (auto it : comp_kernel_map){
      if (comp_kernel_select_count==0) break;
      if (active_kernels[it.second.val_index].num_local_schedules == 0) break;
      auto offset = 18+9*count;
      cp_costs[offset] = it.first.tag;
      cp_costs[offset+1] = it.first.param1;
      cp_costs[offset+2] = it.first.param2;
      cp_costs[offset+3] = it.first.param3;
      cp_costs[offset+4] = it.first.param4;
      cp_costs[offset+5] = it.first.param5;
      cp_costs[offset+6] = it.first.kernel_index;
      cp_costs[offset+7] = it.first.flops;
      cp_costs[offset+8] = active_kernels[it.second.val_index].num_local_schedules;
      count++; if (count==comp_kernel_select_count) break;
    }
    count=0;
    for (auto it : comm_kernel_map){
      if (comm_kernel_select_count==0) break;
      if (active_kernels[it.second.val_index].num_local_schedules == 0) break;
      auto offset = 18+9*comp_kernel_select_count+count*9;
      cp_costs[offset] = it.first.tag;
      cp_costs[offset+1] = it.first.dim_sizes[0];
      cp_costs[offset+2] = it.first.dim_sizes[1];
      cp_costs[offset+3] = it.first.dim_strides[0];
      cp_costs[offset+4] = it.first.dim_strides[1];
      cp_costs[offset+5] = it.first.partner_offset;
      cp_costs[offset+6] = it.first.kernel_index;
      cp_costs[offset+7] = it.first.msg_size;
      cp_costs[offset+8] = active_kernels[it.second.val_index].num_local_schedules;
      count++; if (count==comm_kernel_select_count) break;
    }
  }

  if (partner1 == -1){
    MPI_Op op;
    MPI_Op_create((MPI_User_function*) propagate_cp_op,0,&op);
    PMPI_Allreduce(MPI_IN_PLACE, &cp_costs[0], cp_costs.size(), MPI_FLOAT, op, tracker.comm);
    MPI_Op_free(&op);
    if (collective_state_protocol) schedule_decision = (cp_costs[6] == 0 ? false : true);
    else schedule_decision = (cp_costs[5] == 0 ? true : false);
    tracker.aggregate_comp_kernels = cp_costs[7]>0;
    tracker.aggregate_comm_kernels = cp_costs[8]>0;
    tracker.should_propagate = cp_costs[7]>0 || cp_costs[8]>0;
    if (comm_kernel_map.find(key) != comm_kernel_map.end()){
      if (!schedule_decision){
        if (should_schedule(comm_kernel_map[key])){
          if (comm_kernel_transfer_id == 1){
            // Enter here if a process's local comm kernel is not globally steady, yet its found that at least one of the processors in its communicator is.
            // Completely swap out its kernel statistics for the elements reduced. Note that I am avoiding an extra broadcast, thus the reduced members might
            //   each be from different processors. Likely though, just one is globally steady.
            active_kernels[comm_kernel_map[key].val_index].hash_id = cp_costs[9];
            active_kernels[comm_kernel_map[key].val_index].num_schedules = cp_costs[10];
            active_kernels[comm_kernel_map[key].val_index].num_local_schedules = cp_costs[11];
            active_kernels[comm_kernel_map[key].val_index].num_scheduled_units = cp_costs[12];
            active_kernels[comm_kernel_map[key].val_index].num_local_scheduled_units = cp_costs[13];
            active_kernels[comm_kernel_map[key].val_index].M1 = cp_costs[14];
            active_kernels[comm_kernel_map[key].val_index].M2 = cp_costs[15];
            active_kernels[comm_kernel_map[key].val_index].total_exec_time = cp_costs[16];
            active_kernels[comm_kernel_map[key].val_index].total_local_exec_time  = cp_costs[17];
          }
        }
        set_kernel_state(comm_kernel_map[key],false);
        set_kernel_state_global(comm_kernel_map[key],false);
      }
    }
  }
  else{
    bool has_received=false;
    if ((is_sender) && (rank != partner1)){
      MPI_Buffer_attach(&eager_pad[0],eager_pad.size());
      PMPI_Bsend(&cp_costs[0], cp_costs.size(), MPI_FLOAT, partner1, internal_tag2, tracker.comm);
      void* temp_buf; int temp_size;
      MPI_Buffer_detach(&temp_buf,&temp_size);
    }
    if ((!is_sender) && (rank != partner1)){
      has_received=true;
      PMPI_Recv(&cp_costs_foreign[0], cp_costs_foreign.size(), MPI_FLOAT, partner1, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
    }
    if ((partner2 != -1) && (rank != partner2)){
      has_received=true;
      PMPI_Recv(&cp_costs_foreign[0], cp_costs_foreign.size(), MPI_FLOAT, partner2, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
    }
    if (has_received){
      schedule_decision = (cp_costs_foreign[5] == 0);
      update_frequency(&cp_costs_foreign[0],&cp_costs[0],cp_costs_size);
      if (tracker.tag >= 13 && tracker.tag <= 14) schedule_decision = (cp_costs[5] == 0);
      if (comm_kernel_map.find(key) != comm_kernel_map.end()){
        if (!schedule_decision){
          set_kernel_state(comm_kernel_map[key],false);
          set_kernel_state_global(comm_kernel_map[key],false);
        }
      }
    }
  }
  // Only Receivers always update their kernel maps
  if (sample_constraint_mode == 2){
    if (!(tracker.partner1 != -1 && tracker.is_sender && tracker.partner2 == tracker.partner1)) kernel_update(&cp_costs_foreign[18]);
  }

  global_schedule_decision = schedule_decision;// 'global' as in global variable
  intercept_overhead[1] += MPI_Wtime() - overhead_start_time;
  // start communication timer for communication routine
  tracker.start_time = MPI_Wtime();
  return schedule_decision;
}

// Used only for p2p communication. All blocking collectives use sychronous protocol
void path::complete_comm(blocking& tracker, int recv_source){
  volatile auto comm_time = MPI_Wtime() - tracker.start_time;	// complete communication time
  // Special exit if no kernels are to be scheduled -- the goal is to track the total overhead time (no comp/comm kernels), which should
  //   be attained with timers outside of critter.
  if (schedule_kernels==0){ return; }
  volatile auto overhead_start_time = MPI_Wtime();

  int rank; MPI_Comm_rank(tracker.comm,&rank);
  if (tuning_delta > 1){// tuning_delta of 1 signifies that communication kernel scheduling is to be unconditional
    int comm_sizes[2]={0,0}; int comm_strides[2]={0,0};
    for (auto i=0; i<comm_channel_map[tracker.comm]->id.size(); i++){
      comm_sizes[i]=comm_channel_map[tracker.comm]->id[i].first;
      comm_strides[i]=comm_channel_map[tracker.comm]->id[i].second;
    }
    // Below, the idea is that key doesn't exist in comm_kernel_map iff the key hasn't been seen before. If the key has been seen, we automatically
    //   create an entry in comm_kernel_key, although it will be empty.
    comm_kernel_key key(rank,active_kernels.size(),tracker.tag,comm_sizes,comm_strides,tracker.nbytes,tracker.partner1);
    if (comm_kernel_map.find(key) == comm_kernel_map.end()){
      active_comm_kernel_keys.push_back(key);
      if (tracker.partner1 != -1){
        auto world_partner_rank = channel::translate_rank(tracker.comm,tracker.partner1);
        active_kernels.emplace_back();
      } else{
        assert(comm_channel_map.find(tracker.comm) != comm_channel_map.end());
        active_kernels.emplace_back(comm_channel_map[tracker.comm]);
      }
      comm_kernel_map[key] = kernel_key_id(true,active_comm_kernel_keys.size()-1,active_kernels.size()-1,false);
    } else{
      if (global_schedule_decision) comm_kernel_map[key].is_active = true;
    }
    // For debugging sanity, if the autotuner shut off execution of this kernel, then have debugger use only existing mean
    if (autotuning_debug == 1 && comm_kernel_save_map.find(key) != comm_kernel_save_map.end() && should_schedule(comm_kernel_save_map[key]) == 0){
      comm_time = get_estimate(comm_kernel_map[key],comm_analysis_param,tracker.nbytes);
    } else{
      update_kernel_stats(comm_kernel_map[key],comm_analysis_param,comm_time,tracker.nbytes);
      if (should_schedule(comm_kernel_map[key])==0){
        // Note if this is true, the corresponding entry in the batch map must be cleared. However, I think I delete the entire map in aggregation mode 1, so asserting
        //   on this is difficult.
        // Note: I think branching on aggregation mode is not needed. The pathset should contain all batch samples and the batch should be cleared.
        comm_time = get_estimate(comm_kernel_map[key],comm_analysis_param,tracker.nbytes);
      } else{
        if (tracker.partner1 == -1){
          bool _is_steady = steady_test(key,comm_kernel_map[key],comm_analysis_param);
          set_kernel_state(comm_kernel_map[key],!_is_steady);
          if (sample_constraint_mode == -1) set_kernel_state_global(comm_kernel_map[key],!_is_steady);// Force global state to steady.
          if (comm_state_aggregation_mode == 0) { set_kernel_state_global(comm_kernel_map[key],!_is_steady); }
        } else{
          // If p2p (most notably senders) must delay update, we must check whether the kernel has already been set to local steady state
          if (delay_state_update){
            if (active_kernels[comm_kernel_map[key].val_index].steady_state){
              set_kernel_state_global(comm_kernel_map[key],false);
            } else{
              bool _is_steady = steady_test(key,comm_kernel_map[key],comm_analysis_param);
              set_kernel_state(comm_kernel_map[key],!_is_steady);
            }
          } else{
            bool _is_steady = steady_test(key,comm_kernel_map[key],comm_analysis_param);
            set_kernel_state(comm_kernel_map[key],!_is_steady);
            if (sample_constraint_mode == -1) set_kernel_state_global(comm_kernel_map[key],!_is_steady);// Force global state to steady.
            set_kernel_state_global(comm_kernel_map[key],!_is_steady);
          }
        }
      }
    }
  }

  cp_costs[num_cp_measures-1] += comm_time;
  cp_costs[num_cp_measures-4] += comm_time;
  vol_costs[num_vol_measures-1] += comm_time;
  vol_costs[num_vol_measures-4] += comm_time;

  //vol_costs[num_vol_measures-4] -= std::max(0.,vol_costs[num_vol_measures-1]-cp_costs[num_cp_measures-1]);
  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
  vol_costs[num_vol_measures-4] = vol_costs[num_vol_measures-4] > cp_costs[num_cp_measures-4]
                                          ? cp_costs[num_cp_measures-4] : vol_costs[num_vol_measures-4];
  vol_costs[num_vol_measures-3] = vol_costs[num_vol_measures-3] > cp_costs[num_cp_measures-3]
                                          ? cp_costs[num_cp_measures-3] : vol_costs[num_vol_measures-3];
  vol_costs[num_vol_measures-2] = vol_costs[num_vol_measures-2] > cp_costs[num_cp_measures-2]
                                          ? cp_costs[num_cp_measures-2] : vol_costs[num_vol_measures-2];
  vol_costs[num_vol_measures-1] = vol_costs[num_vol_measures-1] > cp_costs[num_cp_measures-1]
                                          ? cp_costs[num_cp_measures-1] : vol_costs[num_vol_measures-1];
  // Propogate critical paths for all processes in communicator based on what each process has seen up until now (not including this communication)
  if (tracker.should_propagate && tracker.partner1 == -1){
    bool is_world_communication = (tracker.comm == MPI_COMM_WORLD) && (tracker.partner1 == -1);
    if ((rank == tracker.partner1) && (rank == tracker.partner2)) { ; }
    else{
      if (tracker.aggregate_comm_kernels) comm_state_aggregation(tracker);
      if (tracker.aggregate_comp_kernels) comp_state_aggregation(tracker);
    }
  }

  tracker.should_propagate = false;
  tracker.aggregate_comp_kernels = false;
  tracker.aggregate_comm_kernels = false;
  intercept_overhead[2] += MPI_Wtime() - overhead_start_time;
  // Prepare to leave interception and re-enter user code by restarting computation timers.
  tracker.start_time = MPI_Wtime();
  computation_timer = tracker.start_time;
}

// Called by both nonblocking p2p and nonblocking collectives
bool path::initiate_comm(nonblocking& tracker, volatile double curtime, int64_t nelem,
                         MPI_Datatype t, MPI_Comm comm, int user_tag, bool is_sender, int partner){
  // Save and accumulate the computation time between last communication routine as both execution-time and computation time
  //   into both the execution-time critical path data structures and the per-process data structures.
  tracker.comp_time = curtime - computation_timer;
  assert(partner != MPI_ANY_SOURCE);
  cp_costs[num_cp_measures-1] += tracker.comp_time;
  cp_costs[num_cp_measures-3] += tracker.comp_time;
  vol_costs[num_vol_measures-1] += tracker.comp_time;
  vol_costs[num_vol_measures-3] += tracker.comp_time;
  // Special exit if no kernels are to be scheduled -- the goal is to track the total overhead time (no comp/comm kernels), which should
  //   be attained with timers outside of critter.
  if (schedule_kernels==0){ return false; }

  // At this point, 'cp_costs' tells me the process's time up until now. A barrier won't suffice when kernels are conditionally scheduled.
  int rank; MPI_Comm_rank(comm, &rank);
  volatile auto overhead_start_time = MPI_Wtime();

  // Save caller communication attributes into reference object for use in corresponding static method 'complete_comm'
  int word_size,np; MPI_Type_size(t, &word_size);
  int64_t nbytes = word_size * nelem;
  MPI_Comm_size(comm, &np);
  tracker.nbytes = nbytes;
  tracker.comm = comm;
  tracker.comm_size = np;
  tracker.is_sender = is_sender;
  tracker.partner1 = partner;
  tracker.partner2 = partner;
  tracker.synch_time = 0.;// might get updated below
  tracker.barrier_time = 0.;// might get updated below
  tracker.should_propagate = false;
  tracker.aggregate_comp_kernels=false;
  tracker.aggregate_comm_kernels=false;

  bool schedule_decision = true;
  // Assume that the communicator of either collective/p2p is registered via comm_split, and that its described using a max of 3 dimension tuples.
  assert(comm_channel_map.find(tracker.comm) != comm_channel_map.end());
  int comm_sizes[2]={0,0}; int comm_strides[2]={0,0};
  for (auto i=0; i<comm_channel_map[tracker.comm]->id.size(); i++){
    comm_sizes[i]=comm_channel_map[tracker.comm]->id[i].first;
    comm_strides[i]=comm_channel_map[tracker.comm]->id[i].second;
  }
  // Below, the idea is that key doesn't exist in comm_kernel_map iff the key hasn't been seen before. If the key has been seen, we automatically
  //   create an entry in comm_kernel_key, although it will be empty.
  comm_kernel_key key(rank,-1,tracker.tag,comm_sizes,comm_strides,tracker.nbytes,tracker.partner1);
  memset(&cp_costs[4],0,14*sizeof(float)); cp_costs[6] = 1;
  if (tuning_delta > 1){// tuning_delta of 1 signifies that communication kernel scheduling is to be unconditional
    if (comm_kernel_map.find(key) != comm_kernel_map.end()){
      schedule_decision = should_schedule(comm_kernel_map[key])==1;
      cp_costs[5] = (!schedule_decision ? 1 : 0);	// This logic must match that in 'initiate_comm(blocking&,...)'
      cp_costs[6] = (schedule_decision ? 1 : 0);	// This logic must match that in 'initiate_comm(blocking&,...)'
      if (!schedule_decision){
        // If this particular kernel is globally steady, meaning it has exhausted its state aggregation channels,
        //   then we can overwrite the '-1' with the sample mean of the globally-steady kernel
        cp_costs[4] = active_kernels[comm_kernel_map[key].val_index].steady_state;
        cp_costs[9] = active_kernels[comm_kernel_map[key].val_index].hash_id;
        cp_costs[10] = active_kernels[comm_kernel_map[key].val_index].num_schedules;
        cp_costs[11] = active_kernels[comm_kernel_map[key].val_index].num_local_schedules;
        cp_costs[12] = active_kernels[comm_kernel_map[key].val_index].num_scheduled_units;
        cp_costs[13] = active_kernels[comm_kernel_map[key].val_index].num_local_scheduled_units;
        cp_costs[14] = active_kernels[comm_kernel_map[key].val_index].M1;
        cp_costs[15] = active_kernels[comm_kernel_map[key].val_index].M2;
        cp_costs[16] = active_kernels[comm_kernel_map[key].val_index].total_exec_time;
        cp_costs[17] = active_kernels[comm_kernel_map[key].val_index].total_local_exec_time;
      }
    }
  }

  // Note: I will not write special case for rank==partner
  if (sample_constraint_mode == 2){
    // Fill in -1 first because the number of distinct kernels might be less than 'comm_kernel_select_count',
    //   just to avoid confusion. A -1 tag clearly means that the kernel is void
    memset(&cp_costs[18],-1,sizeof(float)*(cp_costs.size()-18));
    // Iterate over first 'comp_kernel_select_count' keys
    int count=0;
    for (auto it : comp_kernel_map){
      if (comp_kernel_select_count==0) break;
      if (active_kernels[it.second.val_index].num_local_schedules == 0) break;
      auto offset = 18+9*count;
      cp_costs[offset] = it.first.tag;
      cp_costs[offset+1] = it.first.param1;
      cp_costs[offset+2] = it.first.param2;
      cp_costs[offset+3] = it.first.param3;
      cp_costs[offset+4] = it.first.param4;
      cp_costs[offset+5] = it.first.param5;
      cp_costs[offset+6] = it.first.kernel_index;
      cp_costs[offset+7] = it.first.flops;
      cp_costs[offset+8] = active_kernels[it.second.val_index].num_local_schedules;
      count++; if (count==comp_kernel_select_count) break;
    }
    count=0;
    for (auto it : comm_kernel_map){
      if (comm_kernel_select_count==0) break;
      if (active_kernels[it.second.val_index].num_local_schedules == 0) break;
      auto offset = 18+9*comp_kernel_select_count+count*9;
      cp_costs[offset] = it.first.tag;
      cp_costs[offset+1] = it.first.dim_sizes[0];
      cp_costs[offset+2] = it.first.dim_sizes[1];
      cp_costs[offset+3] = it.first.dim_strides[0];
      cp_costs[offset+4] = it.first.dim_strides[1];
      cp_costs[offset+5] = it.first.partner_offset;
      cp_costs[offset+6] = it.first.kernel_index;
      cp_costs[offset+7] = it.first.msg_size;
      cp_costs[offset+8] = active_kernels[it.second.val_index].num_local_schedules;
      count++; if (count==comm_kernel_select_count) break;
    }
  }
  save_path_data = nullptr; save_prop_req = MPI_REQUEST_NULL;
  if (tracker.partner1 == -1){
    assert(delay_state_update);
    save_path_data = (float*)malloc(cp_costs_size*sizeof(float));
    std::memcpy(save_path_data,&cp_costs[0],cp_costs_size*sizeof(float));
    MPI_Op op;
    MPI_Op_create((MPI_User_function*) propagate_cp_op,0,&op);
    PMPI_Iallreduce(MPI_IN_PLACE, save_path_data, cp_costs_size, MPI_FLOAT, op, tracker.comm, &save_prop_req);
    //MPI_Op_free(&op);
  }
  else{
    if (is_sender){
      MPI_Buffer_attach(&eager_pad[0],eager_pad.size());
      PMPI_Bsend(&cp_costs[0],cp_costs_size,MPI_FLOAT,tracker.partner1,internal_tag2,tracker.comm);
      void* temp_buf; int temp_size;
      MPI_Buffer_detach(&temp_buf,&temp_size);
    } else{
      assert(delay_state_update);
      save_path_data = (float*)malloc(cp_costs_size*sizeof(float));
      PMPI_Irecv(save_path_data, cp_costs_size, MPI_FLOAT, tracker.partner1, internal_tag2, tracker.comm, &save_prop_req);
    }
  }
  if (!schedule_decision){
    while (1){
      if (request_id == INT_MAX) request_id = 100;// reset to avoid overflow. rare case.
      if ((nonblocking_internal_info.find(request_id) == nonblocking_internal_info.end()) && (request_id != MPI_REQUEST_NULL)){
        nonblocking_info msg_info(save_path_data,save_prop_req,false,is_sender,partner,comm,(float)nbytes,(float)np,user_tag,&tracker);
        nonblocking_internal_info[request_id] = msg_info;
        break;
      }
      request_id++;
    }
  }

  global_schedule_decision = schedule_decision;// 'global' as in global variable
  intercept_overhead[1] += MPI_Wtime() - overhead_start_time;
  computation_timer = MPI_Wtime();
  return schedule_decision;
}

// Called by both nonblocking p2p and nonblocking collectives
void path::initiate_comm(nonblocking& tracker, volatile double itime, int64_t nelem,
                         MPI_Datatype t, MPI_Comm comm, MPI_Request* request, int user_tag, bool is_sender, int partner){
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  tracker.comp_time = itime;
  cp_costs[num_cp_measures-3] += tracker.comp_time;
  cp_costs[num_cp_measures-1] += tracker.comp_time;
  vol_costs[num_vol_measures-3] += tracker.comp_time;
  vol_costs[num_vol_measures-1] += tracker.comp_time;

  int el_size,p;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(comm, &p);

  // Note: I will not write special case for rank==partner
  // These asserts are to prevent the situation in which my synthetic request_id and that of the MPI implementation collide
  assert(nonblocking_internal_info.find(*request) == nonblocking_internal_info.end());
  nonblocking_info msg_info(save_path_data,save_prop_req,true,is_sender,partner,comm,(float)nbytes,(float)p,user_tag,&tracker);
  nonblocking_internal_info[*request] = msg_info;
  computation_timer = MPI_Wtime();
}

void path::complete_comm(nonblocking& tracker, float* path_data, MPI_Request* request, double comp_time, double comm_time){
  auto info_it = nonblocking_internal_info.find(*request);
  assert(info_it != nonblocking_internal_info.end());

  tracker.is_sender = info_it->second.is_sender;
  tracker.comm = info_it->second.comm;
  tracker.partner1 = info_it->second.partner;
  tracker.partner2 = -1;
  tracker.nbytes = info_it->second.nbytes;
  tracker.comm_size = info_it->second.comm_size;
  tracker.synch_time=0;

  int rank; MPI_Comm_rank(tracker.comm,&rank);
  // Right now, I need to check schedule_decision and replace 'comm_time' with the predicted time, if necessary.
  if (tuning_delta > 1){// tuning_delta of 1 signifies that communication kernel scheduling is to be unconditional
    int comm_sizes[2]={0,0}; int comm_strides[2]={0,0};
    for (auto i=0; i<comm_channel_map[tracker.comm]->id.size(); i++){
      comm_sizes[i]=comm_channel_map[tracker.comm]->id[i].first;
      comm_strides[i]=comm_channel_map[tracker.comm]->id[i].second;
    }
    // Below, the idea is that key doesn't exist in comm_kernel_map iff the key hasn't been seen before. If the key has been seen, we automatically
    //   create an entry in comm_kernel_key, although it will be empty.
    comm_kernel_key key(rank,active_kernels.size(),tracker.tag,comm_sizes,comm_strides,tracker.nbytes,tracker.partner1);
    if (comm_kernel_map.find(key) == comm_kernel_map.end()){
      active_comm_kernel_keys.push_back(key);
      if (tracker.partner1 != -1){
        auto world_partner_rank = channel::translate_rank(tracker.comm,tracker.partner1);
        active_kernels.emplace_back();
      } else{
        assert(comm_channel_map.find(tracker.comm) != comm_channel_map.end());
        active_kernels.emplace_back(comm_channel_map[tracker.comm]);
      }
      comm_kernel_map[key] = kernel_key_id(true,active_comm_kernel_keys.size()-1,active_kernels.size()-1,false);
    } else{
      if (info_it->second.is_active) comm_kernel_map[key].is_active = true;
    }
    // For debugging sanity, if the autotuner shut off execution of this kernel, then have debugger use only existing mean
    if (autotuning_debug == 1 && comm_kernel_save_map.find(key) != comm_kernel_save_map.end() && should_schedule(comm_kernel_save_map[key]) == 0){
      comm_time = get_estimate(comm_kernel_map[key],comm_analysis_param,tracker.nbytes);
    } else{
      update_kernel_stats(comm_kernel_map[key],comm_analysis_param,comm_time,tracker.nbytes);
      if (should_schedule(comm_kernel_map[key])==0){
        // Note if this is true, the corresponding entry in the batch map must be cleared. However, I think I delete the entire map in aggregation mode 1, so asserting
        //   on this is difficult.
        // Note: I think branching on aggregation mode is not needed. The pathset should contain all batch samples and the batch should be cleared.
        comm_time = get_estimate(comm_kernel_map[key],comm_analysis_param,tracker.nbytes);
      } else{
        if (delay_state_update){
          // Receivers of any kind must transfer the sender's local state first. This then allows a quick jump from active to globally steady
          if (tracker.partner1 == -1 || !tracker.is_sender) set_kernel_state(comm_kernel_map[key],path_data[4]==0);
          if (active_kernels[comm_kernel_map[key].val_index].steady_state){
            set_kernel_state_global(comm_kernel_map[key],false);
          } else{
            bool is_steady = steady_test(key,comm_kernel_map[key],comm_analysis_param);
            set_kernel_state(comm_kernel_map[key],!is_steady);
          }
        } else{
          bool is_steady = steady_test(key,comm_kernel_map[key],comm_analysis_param);
          set_kernel_state(comm_kernel_map[key],!is_steady);
          if (sample_constraint_mode == -1) set_kernel_state_global(comm_kernel_map[key],!is_steady);// Force global state to steady.
          set_kernel_state_global(comm_kernel_map[key],!is_steady);
        }
      }
    }
  }

  cp_costs[num_cp_measures-4] += comm_time;			// update critical path communication time
  cp_costs[num_cp_measures-3] += comp_time;			// update critical path runtime
  cp_costs[num_cp_measures-1] += comp_time+comm_time;		// update critical path runtime

  vol_costs[num_vol_measures-4] += comm_time;				// update local communication time
  vol_costs[num_vol_measures-3] += comp_time;				// update local runtime
  vol_costs[num_vol_measures-1] += comp_time+comm_time;			// update local runtime
  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
  vol_costs[num_vol_measures-4] = vol_costs[num_vol_measures-4] > cp_costs[num_cp_measures-4]
                                          ? cp_costs[num_cp_measures-4] : vol_costs[num_vol_measures-4];
  vol_costs[num_vol_measures-3] = vol_costs[num_vol_measures-3] > cp_costs[num_cp_measures-3]
                                          ? cp_costs[num_cp_measures-3] : vol_costs[num_vol_measures-3];
  vol_costs[num_vol_measures-2] = vol_costs[num_vol_measures-2] > cp_costs[num_cp_measures-2]
                                          ? cp_costs[num_cp_measures-2] : vol_costs[num_vol_measures-2];
  vol_costs[num_vol_measures-1] = vol_costs[num_vol_measures-1] > cp_costs[num_cp_measures-1]
                                          ? cp_costs[num_cp_measures-1] : vol_costs[num_vol_measures-1];
  if (!tracker.is_sender || tracker.partner1==-1) free(info_it->second.path_data);
  nonblocking_internal_info.erase(*request);
}

int path::complete_comm(double curtime, MPI_Request* request, MPI_Status* status){
  auto comp_time = curtime - computation_timer;
  int ret = MPI_SUCCESS;
  assert(nonblocking_internal_info.find(*request) != nonblocking_internal_info.end());
  auto info_it = nonblocking_internal_info.find(*request);
  auto save_r = info_it->first;
  int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
  if (info_it->second.is_active == 1){
    volatile auto last_start_time = MPI_Wtime();
    ret = PMPI_Wait(request, status);
    auto save_comm_time = MPI_Wtime() - last_start_time;
    auto overhead_start_time = MPI_Wtime();
    // If receiver or collective, complete the barrier and the path data propagation
    if (!info_it->second.is_sender || info_it->second.partner==-1){
      assert(info_it->second.path_data != nullptr);
      PMPI_Wait(&info_it->second.prop_req, MPI_STATUS_IGNORE);
      if (info_it->second.partner != -1) update_frequency(info_it->second.path_data,&cp_costs[0],cp_costs_size);
      if (sample_constraint_mode == 2) kernel_update(&info_it->second.path_data[18]);
    }
    complete_comm(*info_it->second.track, info_it->second.path_data, &save_r, comp_time, save_comm_time);
    intercept_overhead[2] += MPI_Wtime() - overhead_start_time;
  } else{
    auto overhead_start_time = MPI_Wtime();
    // If receiver or collective, complete the barrier and the path data propagation
    if (!info_it->second.is_sender || info_it->second.partner==-1){
      assert(info_it->second.path_data != nullptr);
      PMPI_Wait(&info_it->second.prop_req, MPI_STATUS_IGNORE);
      if (info_it->second.partner != -1) update_frequency(info_it->second.path_data,&cp_costs[0],cp_costs_size);
      if (sample_constraint_mode == 2) kernel_update(&info_it->second.path_data[18]);
      if (status != MPI_STATUS_IGNORE){
        status->MPI_SOURCE = info_it->second.partner;
        status->MPI_TAG = info_it->second.tag;
      }
    }
    complete_comm(*info_it->second.track, info_it->second.path_data, &save_r, comp_time, 1000000.);
    *request = MPI_REQUEST_NULL;
    intercept_overhead[2] += MPI_Wtime() - overhead_start_time;
  }
  computation_timer = MPI_Wtime();
  return ret;
}

int path::complete_comm(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){
  auto comp_time = curtime - computation_timer;
  auto overhead_start_time = MPI_Wtime();
  int ret = MPI_SUCCESS;
  std::vector<MPI_Request> pt(count);
  int num_skips=0; int last_skip;
  for (int i=0;i<count;i++){
    assert(nonblocking_internal_info.find((array_of_requests)[i]) == nonblocking_internal_info.end());
    if (nonblocking_internal_info[(array_of_requests)[i]].is_active){
      pt[i] = (array_of_requests)[i];
    } else{
      pt[i] = MPI_REQUEST_NULL; num_skips++; last_skip = i;
    }
  }
  intercept_overhead[2] += MPI_Wtime() - overhead_start_time;
  if (num_skips < pt.size()){
    volatile auto last_start_time = MPI_Wtime();
    ret = PMPI_Waitany(count,array_of_requests,indx,status);
    auto save_comm_time = MPI_Wtime() - last_start_time;
    overhead_start_time = MPI_Wtime();
    auto info_it = nonblocking_internal_info.find((array_of_requests)[*indx]);
    assert(info_it != nonblocking_internal_info.end());
    int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
    // If receiver or collective, complete the barrier and the path data propagation
    if (!info_it->second.is_sender || info_it->second.partner==-1){
      assert(info_it->second.path_data != nullptr);
      PMPI_Wait(&info_it->second.prop_req, MPI_STATUS_IGNORE);
      if (info_it->second.partner != -1) update_frequency(info_it->second.path_data,&cp_costs[0],cp_costs_size);
      if (sample_constraint_mode == 2) kernel_update(&info_it->second.path_data[18]);
    }
    complete_comm(*info_it->second.track, info_it->second.path_data, &(array_of_requests)[*indx], comp_time, save_comm_time);
  } else{
    overhead_start_time = MPI_Wtime();
    auto info_it = nonblocking_internal_info.find((array_of_requests)[last_skip]);
    assert(info_it != nonblocking_internal_info.end());
    int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
    // If receiver or collective, complete the barrier and the path data propagation
    if (!info_it->second.is_sender || info_it->second.partner==-1){
      assert(info_it->second.path_data != nullptr);
      PMPI_Wait(&info_it->second.prop_req, MPI_STATUS_IGNORE);
      if (info_it->second.partner != -1) update_frequency(info_it->second.path_data,&cp_costs[0],cp_costs_size);
      if (sample_constraint_mode == 2) kernel_update(&info_it->second.path_data[18]);
      if (status != MPI_STATUS_IGNORE){
        status->MPI_SOURCE = info_it->second.partner;
        status->MPI_TAG = info_it->second.tag;
      }
    }
    complete_comm(*info_it->second.track, info_it->second.path_data, &(array_of_requests)[last_skip], comp_time, 1000000.);
    (array_of_requests)[last_skip] = MPI_REQUEST_NULL;
  }
  intercept_overhead[2] += MPI_Wtime() - overhead_start_time;
  computation_timer = MPI_Wtime();
  return ret;
}

int path::complete_comm(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[],
                        MPI_Status array_of_statuses[]){
  int indx; MPI_Status stat;
  int ret = complete_comm(curtime,incount,array_of_requests,&indx,&stat);
  if (array_of_statuses != MPI_STATUSES_IGNORE) array_of_statuses[indx] = stat;
  array_of_indices[0] = indx;
  *outcount=1;
  return ret;
}

int path::complete_comm(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
  auto comp_time = curtime - computation_timer;
  auto overhead_start_time = MPI_Wtime();
  int ret = MPI_SUCCESS;
  int true_count = count;
  std::vector<MPI_Request> pt(count); for (int i=0;i<count;i++){pt[i]=(array_of_requests)[i];}
  // Scan over the requests to identify those that are 'fake'
  for (int i=0; i<count; i++){
    MPI_Request request = array_of_requests[i];
    // 1. check if this request is fake. If so, update its status if a receiver and set its request to MPI_REQUEST_NULL and decrement a count
    auto info_it = nonblocking_internal_info.find(request);
    assert(info_it != nonblocking_internal_info.end());
    int schedule_decision = info_it->second.is_active;
    if (schedule_decision == 0){
      true_count--;
      int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
      if (rank != info_it->second.partner){
        // If receiver or collective, complete the barrier and the path data propagation
        if (!info_it->second.is_sender || info_it->second.partner==-1){
          assert(info_it->second.path_data != nullptr);
          PMPI_Wait(&info_it->second.prop_req, MPI_STATUS_IGNORE);
          update_frequency(info_it->second.path_data,&cp_costs[0],cp_costs_size);
          if (sample_constraint_mode == 2) kernel_update(&info_it->second.path_data[18]);
          if (array_of_statuses != MPI_STATUSES_IGNORE){
            array_of_statuses[i].MPI_SOURCE = info_it->second.partner;
            array_of_statuses[i].MPI_TAG = info_it->second.tag;
          }
        }
      }
      complete_comm(*info_it->second.track, info_it->second.path_data, &pt[i], comp_time, (float)0.);
      comp_time=0;
      array_of_requests[i] = MPI_REQUEST_NULL;
    }
  }
  intercept_overhead[2] += MPI_Wtime() - overhead_start_time;
  // If no requests are fake, issue the user communication the safe way.
  if (true_count == count){
    volatile auto last_start_time = MPI_Wtime();
    ret = PMPI_Waitall(count,array_of_requests,array_of_statuses);
    auto waitall_comm_time = MPI_Wtime() - last_start_time;
    overhead_start_time = MPI_Wtime();
    for (int i=0; i<count; i++){
      MPI_Request request = pt[i];
      auto info_it = nonblocking_internal_info.find(request);
      assert(info_it != nonblocking_internal_info.end());
      int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
      // If receiver or collective, complete the barrier and the path data propagation
      if (!info_it->second.is_sender || info_it->second.partner==-1){
        assert(info_it->second.path_data != nullptr);
        PMPI_Wait(&info_it->second.prop_req, MPI_STATUS_IGNORE);
        if (info_it->second.partner != -1) update_frequency(info_it->second.path_data,&cp_costs[0],cp_costs_size);
        if (sample_constraint_mode == 2) kernel_update(&info_it->second.path_data[18]);
      }
      complete_comm(*info_it->second.track, info_it->second.path_data, &pt[i], comp_time, waitall_comm_time);
      // Although we have to exchange the path data for each request, we do not want to float-count the computation time nor the communicaion time
      comp_time=0; waitall_comm_time=0;
    }
    intercept_overhead[2] += MPI_Wtime() - overhead_start_time;
  }
  else{
    while (true_count>0){
      int indx; MPI_Status status;
      volatile auto start_comm_time = MPI_Wtime();
      int _ret = PMPI_Waitany(count,array_of_requests,&indx,&status);
      auto comm_time = MPI_Wtime() - start_comm_time;
      assert(_ret == MPI_SUCCESS);
      overhead_start_time = MPI_Wtime();
      auto info_it = nonblocking_internal_info.find(pt[indx]);
      assert(info_it != nonblocking_internal_info.end());
      int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
      // If receiver or collective, complete the barrier and the path data propagation
      if (!info_it->second.is_sender || info_it->second.partner==-1){
        assert(info_it->second.path_data != nullptr);
        PMPI_Wait(&info_it->second.prop_req, MPI_STATUS_IGNORE);
        update_frequency(info_it->second.path_data,&cp_costs[0],cp_costs_size);
        if (sample_constraint_mode == 2) kernel_update(&info_it->second.path_data[18]);
        if (array_of_statuses != MPI_STATUSES_IGNORE){
          array_of_statuses[indx].MPI_SOURCE = info_it->second.partner;
          array_of_statuses[indx].MPI_TAG = info_it->second.tag;
        }
      }
      auto save_r = info_it->first;
      complete_comm(*info_it->second.track, info_it->second.path_data, &save_r, comp_time, comm_time);
      comp_time=0;
      array_of_requests[indx] = MPI_REQUEST_NULL;
      true_count--;
      intercept_overhead[2] += MPI_Wtime() - overhead_start_time;
    }
  }
  computation_timer = MPI_Wtime();
  return ret;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void path::comp_state_aggregation(blocking& tracker){
  int size; MPI_Comm_size(tracker.comm,&size);
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  std::vector<kernel_propagate> foreign_active_kernels;
  std::map<comp_kernel_key,kernel_propagate> save_comp_kernels;

  // First save the kernels we want to contribute to the aggregation (because they are steady)
  for (auto& it : tracker.save_comp_key){
    save_comp_kernels[it] = active_kernels[comp_kernel_map[it].val_index];
  }

  size_t active_size = size;
  size_t active_rank = rank;
  size_t active_mult = 1;
  while (active_size>1){
    if (active_rank % 2 == 1){
      // Fill-in the associated comp kernel
      tracker.save_comp_key.clear();
      foreign_active_kernels.clear();
      for (auto& it : save_comp_kernels){
        tracker.save_comp_key.push_back(it.first);
        foreign_active_kernels.push_back(it.second);
      }

      int partner = (active_rank-1)*active_mult;
      int size_array[2] = {tracker.save_comp_key.size(),tracker.save_comp_key.size()};
      // Send sizes before true message so that receiver can be aware of the array sizes for subsequent communication
      PMPI_Send(&size_array[0],2,MPI_INT,partner,internal_tag,tracker.comm);
      // Send active kernels with keys
      PMPI_Send(&tracker.save_comp_key[0],size_array[0],comp_kernel_key_type,partner,internal_tag2,tracker.comm);
      PMPI_Send(&foreign_active_kernels[0],size_array[1],kernel_type,partner,internal_tag2,tracker.comm);
      break;// Incredibely important. Senders must not update {active_size,active_rank,active_mult}
    }
    else if ((active_rank % 2 == 0) && (active_rank < (active_size-1))){
      int partner = (active_rank+1)*active_mult;
      int size_array[2] = {0,0};
      // Recv sizes of arrays to create buffers for subsequent communication
      PMPI_Recv(&size_array[0],2,MPI_INT,partner,internal_tag,tracker.comm,MPI_STATUS_IGNORE);
      // Recv partner's active kernels with keys
      tracker.save_comp_key.resize(size_array[0]);
      foreign_active_kernels.resize(size_array[1]);
      PMPI_Recv(&tracker.save_comp_key[0],size_array[0],comp_kernel_key_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      PMPI_Recv(&foreign_active_kernels[0],size_array[1],kernel_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      // Iterate over all active kernels and simply perform an AND operation on whether a kernel is in steady state.
      //   If just one is active across the world communicator, the kernel must remain active.
      //   If kernel does not exist among the sent kernels, it does not count as active. The logical operation is a trivial (AND 1)
      for (auto i=0; i<tracker.save_comp_key.size(); i++){
        auto& key = tracker.save_comp_key[i];
        if (save_comp_kernels.find(key) != save_comp_kernels.end()){
          auto ci_local = get_error_estimate(key,save_comp_kernels[key],comp_analysis_param);
          auto ci_foreign = get_error_estimate(key,foreign_active_kernels[i],comp_analysis_param);
          if (ci_foreign < ci_local){
            save_comp_kernels[key] = foreign_active_kernels[i];
          }
        } else{
          save_comp_kernels[key] = foreign_active_kernels[i];
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
    foreign_active_kernels.clear();
    for (auto& it : save_comp_kernels){
      tracker.save_comp_key.push_back(it.first);
      foreign_active_kernels.push_back(it.second);
    }
  }
  int size_array[2] = {rank==0 ? tracker.save_comp_key.size() : 0,
                       rank==0 ? tracker.save_comp_key.size() : 0};
  PMPI_Bcast(&size_array[0],2,MPI_INT,0,tracker.comm);
  if (rank != 0){
    tracker.save_comp_key.resize(size_array[0]);
    foreign_active_kernels.resize(size_array[1]);
  }
  PMPI_Bcast(&tracker.save_comp_key[0],size_array[0],comp_kernel_key_type,0,tracker.comm);
  PMPI_Bcast(&foreign_active_kernels[0],size_array[1],kernel_type,0,tracker.comm);
  for (auto i=0; i<tracker.save_comp_key.size(); i++){
    auto& key = tracker.save_comp_key[i];
    if (comp_kernel_map.find(key) != comp_kernel_map.end()){
      if (comp_kernel_transfer_id==0){
        active_kernels[comp_kernel_map[key].val_index].hash_id = foreign_active_kernels[i].hash_id;
      }
      else if (comp_kernel_transfer_id==1){
        auto save_kernel = active_kernels[comp_kernel_map[key].val_index];
        active_kernels[comp_kernel_map[key].val_index] = foreign_active_kernels[i];
        active_kernels[comp_kernel_map[key].val_index].num_non_schedules = save_kernel.num_non_schedules;
        active_kernels[comp_kernel_map[key].val_index].num_non_scheduled_units = save_kernel.num_non_scheduled_units;
      }
    } else{
      // Add new entry.
      active_comp_kernel_keys.push_back(key);
      active_kernels.emplace_back(foreign_active_kernels[i]);
      comp_kernel_map[key] = kernel_key_id(true,active_comp_kernel_keys.size()-1,active_kernels.size()-1,false);
    }
    set_kernel_state(comp_kernel_map[key],false);
    if (comp_state_aggregation_mode==1){
      set_kernel_state_global(comp_kernel_map[key],false);
    } else if (comp_state_aggregation_mode==2){
      if (aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->is_final){
        active_kernels[comp_kernel_map[key].val_index].hash_id = aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag;
        active_kernels[comp_kernel_map[key].val_index].registered_channels.clear();
        active_kernels[comp_kernel_map[key].val_index].registered_channels.insert(comm_channel_map[tracker.comm]);// Add the solo channel, not the aggregate
        assert(aggregate_channel_map.find(active_kernels[comp_kernel_map[key].val_index].hash_id) != aggregate_channel_map.end());
        set_kernel_state_global(comp_kernel_map[key],false);
      }
      else{
        active_kernels[comp_kernel_map[key].val_index].hash_id ^= aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag;
        active_kernels[comp_kernel_map[key].val_index].registered_channels.insert(comm_channel_map[tracker.comm]);// Add the solo channel, not the aggregate
        assert(aggregate_channel_map.find(active_kernels[comp_kernel_map[key].val_index].hash_id) != aggregate_channel_map.end());
        if (aggregate_channel_map[active_kernels[comp_kernel_map[key].val_index].hash_id]->is_final){
          set_kernel_state_global(comp_kernel_map[key],false);
        }
      }
    }
  }
  tracker.save_comp_key.clear();
}

void path::comm_state_aggregation(blocking& tracker){
  int size; MPI_Comm_size(tracker.comm,&size);
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  std::vector<kernel_propagate> foreign_active_kernels;
  std::map<comm_kernel_key,kernel_propagate> save_comm_kernels;

  // First save the kernels we want to contribute to the aggregation (because they are steady)
  for (auto& it : tracker.save_comm_key){
    save_comm_kernels[it] = active_kernels[comm_kernel_map[it].val_index];
  }

  size_t active_size = size;
  size_t active_rank = rank;
  size_t active_mult = 1;
  while (active_size>1){
    if (active_rank % 2 == 1){
      // Fill-in the associated comm kernel
      tracker.save_comm_key.clear();
      foreign_active_kernels.clear();
      for (auto& it : save_comm_kernels){
        tracker.save_comm_key.push_back(it.first);
        foreign_active_kernels.push_back(it.second);
      }

      int partner = (active_rank-1)*active_mult;
      int size_array[2] = {tracker.save_comm_key.size(),tracker.save_comm_key.size()};
      // Send sizes before true message so that receiver can be aware of the array sizes for subsequent communication
      PMPI_Send(&size_array[0],2,MPI_INT,partner,internal_tag,tracker.comm);
      // Send active kernels with keys
      PMPI_Send(&tracker.save_comm_key[0],size_array[0],comm_kernel_key_type,partner,internal_tag2,tracker.comm);
      PMPI_Send(&foreign_active_kernels[0],size_array[1],kernel_type,partner,internal_tag2,tracker.comm);
      break;// Incredibely important. Senders must not update {active_size,active_rank,active_mult}
    }
    else if ((active_rank % 2 == 0) && (active_rank < (active_size-1))){
      int partner = (active_rank+1)*active_mult;
      int size_array[2] = {0,0};
      // Recv sizes of arrays to create buffers for subsequent communication
      PMPI_Recv(&size_array[0],2,MPI_INT,partner,internal_tag,tracker.comm,MPI_STATUS_IGNORE);
      // Recv partner's active kernels with keys
      tracker.save_comm_key.resize(size_array[0]);
      foreign_active_kernels.resize(size_array[1]);
      PMPI_Recv(&tracker.save_comm_key[0],size_array[0],comm_kernel_key_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      PMPI_Recv(&foreign_active_kernels[0],size_array[1],kernel_type,partner,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      // Iterate over all active kernels and simply perform an AND operation on whether a kernel is in steady state.
      //   If just one is active across the world communicator, the kernel must remain active.
      //   If kernel does not exist among the sent kernels, it does not count as active. The logical operation is a trivial (AND 1)
      for (auto i=0; i<tracker.save_comm_key.size(); i++){
        auto& key = tracker.save_comm_key[i];
        if (save_comm_kernels.find(key) != save_comm_kernels.end()){
          auto ci_local = get_error_estimate(key,save_comm_kernels[key],comm_analysis_param);
          auto ci_foreign = get_error_estimate(key,foreign_active_kernels[i],comm_analysis_param);
          if (ci_foreign < ci_local){
            save_comm_kernels[key] = foreign_active_kernels[i];
          }
        } else{
          save_comm_kernels[key] = foreign_active_kernels[i];
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
    foreign_active_kernels.clear();
    for (auto& it : save_comm_kernels){
      tracker.save_comm_key.push_back(it.first);
      foreign_active_kernels.push_back(it.second);
    }
  }
  int size_array[2] = {rank==0 ? tracker.save_comm_key.size() : 0,
                       rank==0 ? tracker.save_comm_key.size() : 0};
  PMPI_Bcast(&size_array[0],2,MPI_INT,0,tracker.comm);
  if (rank != 0){
    tracker.save_comm_key.resize(size_array[0]);
    foreign_active_kernels.resize(size_array[1]);
  }
  PMPI_Bcast(&tracker.save_comm_key[0],size_array[0],comm_kernel_key_type,0,tracker.comm);
  PMPI_Bcast(&foreign_active_kernels[0],size_array[1],kernel_type,0,tracker.comm);
  for (auto i=0; i<tracker.save_comm_key.size(); i++){
    auto& key = tracker.save_comm_key[i];
    if (comm_kernel_map.find(key) != comm_kernel_map.end()){
      if (comm_kernel_transfer_id==0){
        active_kernels[comm_kernel_map[key].val_index].hash_id = foreign_active_kernels[i].hash_id;
      }
      else if (comm_kernel_transfer_id==1){
        auto save_kernel = active_kernels[comm_kernel_map[key].val_index];
        active_kernels[comm_kernel_map[key].val_index] = foreign_active_kernels[i];
        active_kernels[comm_kernel_map[key].val_index].num_non_schedules = save_kernel.num_non_schedules;
        active_kernels[comm_kernel_map[key].val_index].num_non_scheduled_units = save_kernel.num_non_scheduled_units;
      }
    } else{
      // Add new entry.
      active_comm_kernel_keys.push_back(key);
      active_kernels.emplace_back(foreign_active_kernels[i]);
      comm_kernel_map[key] = kernel_key_id(true,active_comm_kernel_keys.size()-1,active_kernels.size()-1,false);
    }
    set_kernel_state(comm_kernel_map[key],false);
    if (comm_state_aggregation_mode==1){
      set_kernel_state_global(comm_kernel_map[key],false);
    } else if (comm_state_aggregation_mode==2){
      if (aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->is_final){
        active_kernels[comm_kernel_map[key].val_index].hash_id = aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag;
        active_kernels[comm_kernel_map[key].val_index].registered_channels.clear();
        active_kernels[comm_kernel_map[key].val_index].registered_channels.insert(comm_channel_map[tracker.comm]);// Add the solo channel, not the aggregate
        assert(aggregate_channel_map.find(active_kernels[comm_kernel_map[key].val_index].hash_id) != aggregate_channel_map.end());
        set_kernel_state_global(comm_kernel_map[key],false);
      }
      else{
        active_kernels[comm_kernel_map[key].val_index].hash_id ^= aggregate_channel_map[comm_channel_map[tracker.comm]->global_hash_tag]->global_hash_tag;
        active_kernels[comm_kernel_map[key].val_index].registered_channels.insert(comm_channel_map[tracker.comm]);// Add the solo channel, not the aggregate
        assert(aggregate_channel_map.find(active_kernels[comm_kernel_map[key].val_index].hash_id) != aggregate_channel_map.end());
        if (aggregate_channel_map[active_kernels[comm_kernel_map[key].val_index].hash_id]->is_final){
          set_kernel_state_global(comm_kernel_map[key],false);
        }
      }
    }
  }
  tracker.save_comm_key.clear();
}

}
}
}
