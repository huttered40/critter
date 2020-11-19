#include "path.h"
#include "../container/symbol_tracker.h"
#include "../util/util.h"
#include "../../util/util.h"

namespace critter{
namespace internal{
namespace skeletonization{

void path::exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm){}

bool path::initiate_comp(size_t id, double flop_count, int param1, int param2, int param3, int param4, int param5){
  // Always skip, and update statistics in 'complete_comp'
  return false;
}

void path::complete_comp(size_t id, double flop_count, int param1, int param2, int param3, int param4, int param5){
  comp_kernel_key key(active_kernels.size(),id,flop_count,param1,param2,param3,param4,param5);// '-1' argument is arbitrary, does not influence overloaded operators
  if (comp_kernel_map.find(key) == comp_kernel_map.end()){
    active_comp_kernel_keys.push_back(key);
    active_kernels.push_back(1);
    comp_kernel_map[key] = kernel_key_id(true,active_comp_kernel_keys.size()-1,active_kernels.size()-1,false);
  }
  else{
    active_kernels[comp_kernel_map[key].val_index]++;
  }

  critical_path_costs[num_critical_path_measures-1] += flop_count;	// execution cost
  volume_costs[num_volume_measures-1] += flop_count;			// execution cost
}

bool path::initiate_comm(blocking& tracker, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                         bool is_sender, int partner1, int partner2){
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
  return false;
}

void path::complete_comm(blocking& tracker, int recv_source){
  volatile double overhead_start_time = MPI_Wtime();
  std::pair<double,double> cost_bsp    = tracker.cost_func_bsp(tracker.nbytes, tracker.comm_size);
  std::pair<double,double> cost_alphabeta = tracker.cost_func_alphabeta(tracker.nbytes, tracker.comm_size);

  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  int comm_sizes[2]={0,0}; int comm_strides[2]={0,0};
  assert(comm_channel_map.find(tracker.comm) != comm_channel_map.end());
  for (auto i=0; i<comm_channel_map[tracker.comm]->id.size(); i++){
    comm_sizes[i]=comm_channel_map[tracker.comm]->id[i].first;
    comm_strides[i]=comm_channel_map[tracker.comm]->id[i].second;
  }
  // Below, the idea is that key doesn't exist in comm_kernel_map iff the key hasn't been seen before. If the key has been seen, we automatically
  //   create an entry in comm_kernel_key, although it will be empty.
  comm_kernel_key key(rank,active_kernels.size(),tracker.tag,comm_sizes,comm_strides,tracker.nbytes,tracker.partner1);
  //if (world_rank == 8) std::cout << "******" << key.tag << " " << key.dim_sizes[0] << " " << key.dim_strides[0] << " " << key.msg_size << " " << key.partner_offset << ")\n";
  if (comm_kernel_map.find(key) == comm_kernel_map.end()){
    active_comm_kernel_keys.push_back(key);
    active_kernels.push_back(1);
    comm_kernel_map[key] = kernel_key_id(true,active_comm_kernel_keys.size()-1,active_kernels.size()-1,false);
  }
  else{
    active_kernels[comm_kernel_map[key].val_index]++;
  }

  critical_path_costs[num_critical_path_measures-1] += cost_alphabeta.first*1e6 + cost_alphabeta.second*1e3;
  volume_costs[num_volume_measures-1] += cost_alphabeta.first*1e6 + cost_alphabeta.second*1e3;

  // Temporary
  if (recv_kernel_override==1){
    propagate_kernels(tracker);
  }
  comm_intercept_overhead_stage2 += MPI_Wtime() - overhead_start_time;
}

// Called by both nonblocking p2p and nonblocking collectives
bool path::initiate_comm(nonblocking& tracker, int64_t nelem, MPI_Datatype t, MPI_Comm comm, bool is_sender, int partner){
  volatile double overhead_start_time = MPI_Wtime();
  // Save caller communication attributes into reference object for use in corresponding static method 'complete_comm'
  int word_size,np; MPI_Type_size(t, &word_size);
  int64_t nbytes = word_size * nelem;
  MPI_Comm_size(comm, &np);
  tracker.nbytes = nbytes;
  tracker.comm = comm;
  tracker.comm_size = np;
  tracker.is_sender = is_sender;
  tracker.partner1 = partner;
  tracker.partner2 = -1;

  std::pair<double,double> cost_bsp    = tracker.cost_func_bsp(tracker.nbytes, tracker.comm_size);
  std::pair<double,double> cost_alphabeta = tracker.cost_func_alphabeta(tracker.nbytes, tracker.comm_size);

  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  int comm_sizes[2]={0,0}; int comm_strides[2]={0,0};
  assert(comm_channel_map.find(tracker.comm) != comm_channel_map.end());
  for (auto i=0; i<comm_channel_map[tracker.comm]->id.size(); i++){
    comm_sizes[i]=comm_channel_map[tracker.comm]->id[i].first;
    comm_strides[i]=comm_channel_map[tracker.comm]->id[i].second;
  }
  // Below, the idea is that key doesn't exist in comm_kernel_map iff the key hasn't been seen before. If the key has been seen, we automatically
  //   create an entry in comm_kernel_key, although it will be empty.
  comm_kernel_key key(rank,active_kernels.size(),tracker.tag,comm_sizes,comm_strides,tracker.nbytes,tracker.partner1);
  //if (world_rank == 8) std::cout << "******" << key.tag << " " << key.dim_sizes[0] << " " << key.dim_strides[0] << " " << key.msg_size << " " << key.partner_offset << ")\n";
  if (comm_kernel_map.find(key) == comm_kernel_map.end()){
    active_comm_kernel_keys.push_back(key);
    active_kernels.push_back(1);
    comm_kernel_map[key] = kernel_key_id(true,active_comm_kernel_keys.size()-1,active_kernels.size()-1,false);
  }
  else{
    active_kernels[comm_kernel_map[key].val_index]++;
  }

  critical_path_costs[num_critical_path_measures-1] += cost_alphabeta.first*1e6 + cost_alphabeta.second*1e3;
  volume_costs[num_volume_measures-1] += cost_alphabeta.first*1e6 + cost_alphabeta.second*1e3;

  //propagate_kernels(tracker);
  comm_intercept_overhead_stage2 += MPI_Wtime() - overhead_start_time;
  return false;
}

// Called by both nonblocking p2p and nonblocking collectives
void path::initiate_comm(nonblocking& tracker, int64_t nelem,
                         MPI_Datatype t, MPI_Comm comm, MPI_Request* request, bool is_sender, int partner){
  // No-op - this function should never be invoked.
  assert(0);
}

void path::complete_comm(nonblocking& tracker, MPI_Request* request){}

int path::complete_comm(MPI_Request* request, MPI_Status* status){
  return MPI_SUCCESS;
}

int path::complete_comm(int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){
  return MPI_SUCCESS;
}

int path::complete_comm(int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[],
                        MPI_Status array_of_statuses[]){
  return MPI_SUCCESS;
}

int path::complete_comm(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
  return MPI_SUCCESS;
}

void path::propagate_kernels(blocking& tracker){
  // Use info_receiver[last].second when deciding who to issue 3 broadcasts from
  // First need to broadcast the size of each of the 3 broadcasts so that the receiving buffers can prepare the size of their receiving buffers
  // Only the active kernels need propagating. Steady-state are treated differently depending on the communicator.
  int rank; MPI_Comm_rank(tracker.comm,&rank);

  // Attain the min and max costs to determine whether or not to propogate the root's kernel counts to the rest of the processors
  int max_proc,max_val;
  int min_proc,min_val;
  info_sender[0].first = critical_path_costs[0];
  info_sender[0].second = rank;
  if (tracker.partner1 == -1){
    PMPI_Allreduce(&info_sender[0].first, &info_receiver[0].first, 1, MPI_DOUBLE_INT, MPI_MAXLOC, tracker.comm);
    max_proc = info_receiver[0].second;
    max_val = info_receiver[0].first;
    PMPI_Allreduce(&info_sender[0].first, &info_receiver[0].first, 1, MPI_DOUBLE_INT, MPI_MINLOC, tracker.comm);
    min_proc = info_receiver[0].second;
    min_val = info_receiver[0].first;
    // This only works for blocking. This won't even work for blocking+nonblocking (isend+recv)
    if (max_val == min_val){ return; }
  }
  else{
    if ((recv_kernel_override) || ((tracker.tag >= 13) && (tracker.tag <= 14))){
      PMPI_Sendrecv(&info_sender[0].first, 1, MPI_DOUBLE_INT, tracker.partner1, internal_tag,
                    &info_receiver[0].first, 1, MPI_DOUBLE_INT, tracker.partner2, internal_tag, tracker.comm, MPI_STATUS_IGNORE);
      if (info_sender[0].first>info_receiver[0].first){max_proc = rank;}
      else if (info_sender[0].first==info_receiver[0].first){ max_proc = std::min(rank,tracker.partner1); }
      max_val = std::max(info_sender[0].first, info_receiver[0].first);
      if (info_sender[0].first<info_receiver[0].first){min_proc = rank;}
      else if (info_sender[0].first==info_receiver[0].first){ min_proc = std::min(rank,tracker.partner1); }
      min_val = std::min(info_sender[0].first, info_receiver[0].first);
      // This only works for blocking. This won't even work for blocking+nonblocking (isend+recv)
      if (max_val == min_val){ return; }
    }
    else{
      // Here, knowledge of max/min proc and val is not important. We know the sender will receive additional information from the sender, and vice versa.
      if ((tracker.is_sender) && (rank != tracker.partner1)){
        if (0/*true_eager_p2p*/) { /*PMPI_Bsend(&sbuf, 17, MPI_DOUBLE, partner1, internal_tag3, comm);*/ assert(0); }
        else                { PMPI_Ssend(&info_sender[0].first, 1, MPI_DOUBLE_INT, tracker.partner1, internal_tag2, tracker.comm); }
      }
      if ((!tracker.is_sender) && (rank != tracker.partner1)){
        PMPI_Recv(&info_receiver[0].first, 1, MPI_DOUBLE_INT, tracker.partner1, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
      }
    }
  }

  int size_array[3] = {0,0,0};
  if (rank == max_proc){
    size_array[0] = active_comm_kernel_keys.size();
    size_array[1] = active_comp_kernel_keys.size();
    size_array[2] = active_kernels.size();
  }
  if (tracker.partner1 == -1){
    PMPI_Bcast(&size_array[0],3,MPI_INT,max_proc,tracker.comm);
    if (rank != max_proc){
        active_comm_kernel_keys.resize(size_array[0]);
        active_comp_kernel_keys.resize(size_array[1]);
        active_kernels.resize(size_array[2]);
    }
    PMPI_Bcast(&active_comm_kernel_keys[0],size_array[0],comm_kernel_key_type,max_proc,tracker.comm);
    PMPI_Bcast(&active_comp_kernel_keys[0],size_array[1],comp_kernel_key_type,max_proc,tracker.comm);
    PMPI_Bcast(&active_kernels[0],size_array[2],MPI_INT,max_proc,tracker.comm);
    // receivers update their maps and data structures. Yes, this is costly, but I guess this is what has to happen.
    // basically everything is set via the broadcasts themselves, except for the two maps.
  }
  else{
    // We leverage the fact that we know the path-defining process rank
    //TODO: Note this will not work if process sendrecvs to itself. It only works for CholInv because the cp cost is the same and so it returns before
    //        reaching this part.
    if (rank == max_proc){
      PMPI_Send(&size_array[0],3,MPI_INT,tracker.partner1,internal_tag2,tracker.comm);
      PMPI_Send(&active_comm_kernel_keys[0],size_array[0],comm_kernel_key_type,tracker.partner1,internal_tag2,tracker.comm);
      PMPI_Send(&active_comp_kernel_keys[0],size_array[1],comp_kernel_key_type,tracker.partner1,internal_tag2,tracker.comm);
      PMPI_Send(&active_kernels[0],size_array[2],MPI_INT,tracker.partner1,internal_tag2,tracker.comm);
    } else{
      //TODO: I want to keep both variants. The first is easy to rewrite, just take from above.
      PMPI_Recv(&size_array[0],3,MPI_INT,tracker.partner1,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      active_comm_kernel_keys.resize(size_array[0]);
      active_comp_kernel_keys.resize(size_array[1]);
      active_kernels.resize(size_array[2]);
      PMPI_Recv(&active_comm_kernel_keys[0],size_array[0],comm_kernel_key_type,tracker.partner1,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      PMPI_Recv(&active_comp_kernel_keys[0],size_array[1],comp_kernel_key_type,tracker.partner1,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
      PMPI_Recv(&active_kernels[0],size_array[2],MPI_INT,tracker.partner1,internal_tag2,tracker.comm,MPI_STATUS_IGNORE);
    }
  }

  comm_kernel_map.clear();
  comp_kernel_map.clear();

  // Lets have all processes update, even the root, so that they leave this routine (and subsequently leave the interception) at approximately the same time.
  for (auto i=0; i<active_comm_kernel_keys.size(); i++){
    comm_kernel_key id(active_comm_kernel_keys[i].kernel_index,active_comm_kernel_keys[i].tag,active_comm_kernel_keys[i].dim_sizes,
                        active_comm_kernel_keys[i].dim_strides,active_comm_kernel_keys[i].msg_size,active_comm_kernel_keys[i].partner_offset); 
    if (comm_kernel_map.find(id) == comm_kernel_map.end()){
      comm_kernel_map[id] = kernel_key_id(true, i, active_comm_kernel_keys[i].kernel_index, false);
    } else{
      comm_kernel_map[id].is_active = true;	// always assumed true
      comm_kernel_map[id].key_index = i;
      comm_kernel_map[id].val_index = active_comm_kernel_keys[i].kernel_index;
      comm_kernel_map[id].is_updated = false;
    }
  }

  for (auto i=0; i<active_comp_kernel_keys.size(); i++){
    comp_kernel_key id(active_comp_kernel_keys[i].kernel_index,active_comp_kernel_keys[i].tag,active_comp_kernel_keys[i].flops,
                               active_comp_kernel_keys[i].param1,active_comp_kernel_keys[i].param2,active_comp_kernel_keys[i].param3,
                               active_comp_kernel_keys[i].param4,active_comp_kernel_keys[i].param5); 
    if (comp_kernel_map.find(id) == comp_kernel_map.end()){
      comp_kernel_map[id] = kernel_key_id(true, i, active_comp_kernel_keys[i].kernel_index,false);
    } else{
      comp_kernel_map[id].is_active = true;	// always assumed true
      comp_kernel_map[id].key_index = i;
      comp_kernel_map[id].val_index = active_comp_kernel_keys[i].kernel_index;
      comp_kernel_map[id].is_updated = false;
    }
  }
}

}
}
}
