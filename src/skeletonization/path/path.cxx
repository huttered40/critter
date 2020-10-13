#include "path.h"
#include "../container/symbol_tracker.h"
#include "../util/util.h"
#include "../../util/util.h"

namespace critter{
namespace internal{
namespace skeletonization{

void path::exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm){}

bool path::initiate_comp(size_t id, volatile double curtime, double flop_count, int param1, int param2, int param3, int param4, int param5){
  // Always skip, and update statistics in 'complete_comp'
  return false;
}

void path::complete_comp(size_t id, double flop_count, int param1, int param2, int param3, int param4, int param5){
/*
  comp_pattern_key key(-1,id,flop_count,param1,param2,param3,param4,param5);// '-1' argument is arbitrary, does not influence overloaded operators
  if (comp_pattern_map.find(key) == comp_pattern_map.end()){
    comp_pattern_map[key] = pattern_key_id(true,active_comp_pattern_keys.size()-1,active_patterns.size()-1,false);
  }
  critical_path_costs[num_critical_path_measures-1] += flop_count;	// execution cost
  volume_costs[num_volume_measures-1] += flop_count;			// execution cost
*/
}

bool path::initiate_comm(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                         bool is_sender, int partner1, int partner2){
/*
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
*/
  return false;
}

// Used only for p2p communication. All blocking collectives use sychronous protocol
void path::complete_comm(blocking& tracker, int recv_source){
/*
  volatile double overhead_start_time = MPI_Wtime();
  std::pair<double,double> cost_bsp    = tracker.cost_func_bsp(tracker.nbytes, tracker.comm_size);
  std::pair<double,double> cost_alphabeta = tracker.cost_func_alphabeta(tracker.nbytes, tracker.comm_size);

  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  int comm_sizes[2]={0,0}; int comm_strides[2]={0,0};
  for (auto i=0; i<comm_channel_map[tracker.comm]->id.size(); i++){
    comm_sizes[i]=comm_channel_map[tracker.comm]->id[i].first;
    comm_strides[i]=comm_channel_map[tracker.comm]->id[i].second;
  }
  // Below, the idea is that key doesn't exist in comm_pattern_map iff the key hasn't been seen before. If the key has been seen, we automatically
  //   create an entry in comm_pattern_key, although it will be empty.
  comm_pattern_key key(rank,-1,tracker.tag,comm_sizes,comm_strides,tracker.nbytes,tracker.partner1);
  //if (world_rank == 8) std::cout << "******" << key.tag << " " << key.dim_sizes[0] << " " << key.dim_strides[0] << " " << key.msg_size << " " << key.partner_offset << ")\n";
  if (comm_pattern_map.find(key) == comm_pattern_map.end()){
    comm_pattern_map[key] = pattern_key_id(true,active_comm_pattern_keys.size()-1,active_patterns.size()-1,false);
  }

  critical_path_costs[num_critical_path_measures-1] += cost_alpha_beta.first*1e6 + cost_alpha_beta.second*1e3;
  volume_costs[num_volume_measures-1] += cost_alpha_beta.first*1e6 + cost_alpha_beta.second*1e3;

  // Attain the min and max costs to determine whether or not to propogate the root's kernel counts to the rest of the processors
  PMPI_Allreduce();
  comm_intercept_overhead_stage2 += MPI_Wtime() - overhead_start_time;
*/
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
/*
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
*/

}
}
}
