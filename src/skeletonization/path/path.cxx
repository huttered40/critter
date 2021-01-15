#include "path.h"
#include "../container/symbol_tracker.h"
#include "../util/util.h"
#include "../../util/util.h"

namespace critter{
namespace internal{
namespace skeletonization{

static void update_frequency(float* in, float* inout, size_t len){
  assert(len == cp_costs_size);	// this assert prevents user from obtaining wrong output if MPI implementation cuts up the message.
  if (in[0] > inout[0]){
    std::memcpy(inout,in,len*sizeof(float));
  }
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
      active_kernels[comp_kernel_map[id].val_index] = read_ptr[offset+8];
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
      active_kernels[comm_kernel_map[id].val_index] = read_ptr[offset+8];
    }
  }
}

void path::exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm){
  auto save_comp_time = MPI_Wtime() - computation_timer;
  if (mode==1 && skeleton_type==0){
    cp_costs[0] += save_comp_time;
    vol_costs[0] += save_comp_time;
  }

  generate_aggregate_channels(oldcomm,newcomm);
  PMPI_Barrier(oldcomm);
  if (mode==1){
    computation_timer = MPI_Wtime();
  }
}

bool path::initiate_comp(size_t id, volatile double curtime, float flop_count, int param1, int param2, int param3, int param4, int param5){
  // Always skip, and update statistics in 'complete_comp'
  auto save_comp_time = curtime - computation_timer;
  if (skeleton_type==0){
    cp_costs[0] += save_comp_time;
    vol_costs[0] += save_comp_time;
  }
  comp_start_time = MPI_Wtime();
  if (skeleton_type){ return false; }
  else{ return true; }
}

void path::complete_comp(double errtime, size_t id, float flop_count, int param1, int param2, int param3, int param4, int param5){
  volatile auto comp_time = MPI_Wtime() - comp_start_time - errtime;	// complete computation time
  comp_kernel_key key(active_kernels.size(),id,flop_count,param1,param2,param3,param4,param5);// '-1' argument is arbitrary, does not influence overloaded operators
  if (comp_kernel_map.find(key) == comp_kernel_map.end()){
    active_comp_kernel_keys.push_back(key);
    active_kernels.push_back(1);
    comp_kernel_map[key] = kernel_key_id(true,active_comp_kernel_keys.size()-1,active_kernels.size()-1,false);
  }
  else{
    active_kernels[comp_kernel_map[key].val_index]++;
  }

  if (skeleton_type){
    cp_costs[0] += flop_count;
    vol_costs[0] += flop_count;
  } else{
    cp_costs[0] += comp_time;
    vol_costs[0] += comp_time;
  }
  computation_timer = MPI_Wtime();
}

bool path::initiate_comm(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                         bool is_sender, int partner1, int partner2){
  assert(partner1 != MPI_ANY_SOURCE); if ((tracker.tag == 13) || (tracker.tag == 14)){ assert(partner2 != MPI_ANY_SOURCE); }
  tracker.comp_time = curtime - computation_timer;
  if (skeleton_type==0){
    cp_costs[0] += tracker.comp_time;
    vol_costs[0] += tracker.comp_time;
  }

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

  // Must post these barriers to identify the rank that enters last
  volatile auto init_time = MPI_Wtime();
  if (partner1 == -1){
    PMPI_Barrier(tracker.comm);
  }
  else {
    MPI_Request barrier_reqs[3]; int barrier_count=0;
    char sbuf='H'; char rbuf='H';
    if ((is_sender) && (rank != partner1)){
      MPI_Buffer_attach(&eager_pad[0],eager_pad.size());
      PMPI_Bsend(&sbuf, 1, MPI_CHAR, partner1, internal_tag1, tracker.comm);
      void* temp_buf; int temp_size;
      MPI_Buffer_detach(&temp_buf,&temp_size);
    }
    if ((!is_sender) && (rank != partner1)){
      PMPI_Irecv(&rbuf, 1, MPI_CHAR, partner1, internal_tag1, tracker.comm, &barrier_reqs[barrier_count]); barrier_count++;
    }
    if ((partner2 != -1) && (rank != partner2)){
      PMPI_Irecv(&rbuf, 1, MPI_CHAR, partner2, internal_tag1, tracker.comm, &barrier_reqs[barrier_count]); barrier_count++;
    }
    PMPI_Waitall(barrier_count,&barrier_reqs[0],MPI_STATUSES_IGNORE);
  }

  // No reason to post barrier because we are not interested in decomposing the execution time
  tracker.start_time = MPI_Wtime();
  if (skeleton_type){ return false; }
  else{ return true; }
}

void path::complete_comm(blocking& tracker, int recv_source){
  volatile auto comm_time = MPI_Wtime() - tracker.start_time;	// complete communication time
  std::pair<float,float> cost_alphabeta = tracker.cost_func_alphabeta(tracker.nbytes, tracker.comm_size);

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
  if (comm_kernel_map.find(key) == comm_kernel_map.end()){
    active_comm_kernel_keys.push_back(key);
    active_kernels.push_back(1);
    comm_kernel_map[key] = kernel_key_id(true,active_comm_kernel_keys.size()-1,active_kernels.size()-1,false);
  }
  else{ active_kernels[comm_kernel_map[key].val_index]++; }

  if (skeleton_type==0){
    cp_costs[0] += comm_time;
    vol_costs[0] += comm_time;
    vol_costs[0] = vol_costs[0] > cp_costs[0] ? cp_costs[0] : vol_costs[0];
  } else{
    cp_costs[0] += cost_alphabeta.first*1e6 + cost_alphabeta.second*1e3;
    vol_costs[0] += cost_alphabeta.first*1e6 + cost_alphabeta.second*1e3;
  }

  propagate_kernels(tracker);
  computation_timer = MPI_Wtime();
}

// Called by both nonblocking p2p and nonblocking collectives
bool path::initiate_comm(nonblocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm, int user_tag, bool is_sender, int partner){

  tracker.comp_time = curtime - computation_timer;
  assert(partner != MPI_ANY_SOURCE);
  // Save caller communication attributes into reference object for use in corresponding static method 'complete_comm'
  int rank; MPI_Comm_rank(comm,&rank); 
  int word_size,np; MPI_Type_size(t, &word_size);
  int64_t nbytes = word_size * nelem;
  MPI_Comm_size(comm, &np);
  tracker.nbytes = nbytes;
  tracker.comm = comm;
  tracker.comm_size = np;
  tracker.is_sender = is_sender;
  tracker.partner1 = partner;
  tracker.partner2 = -1;

  if (skeleton_type==0){
    cp_costs[0] += tracker.comp_time;
    vol_costs[0] += tracker.comp_time;
  } else if (rank != partner){
    MPI_Request barrier_req = MPI_REQUEST_NULL;// Only necessary for nonblocking receives
    if (partner!=-1){// Branch protects against nonblocking collectives
      if (is_sender){
        MPI_Buffer_attach(&eager_pad[0],eager_pad.size());
        PMPI_Bsend(&barrier_pad_send, 1, MPI_CHAR, partner, internal_tag1, comm);
        void* temp_buf; int temp_size;
        MPI_Buffer_detach(&temp_buf,&temp_size);
      }
      else{
        PMPI_Irecv(&barrier_pad_recv, 1, MPI_CHAR, partner, internal_tag1, comm, &barrier_req);
      }
    } else{
      PMPI_Iallreduce(&barrier_pad_send, &barrier_pad_recv, 1, MPI_CHAR, MPI_MAX, comm, &barrier_req);
    }
    MPI_Request prop_req = MPI_REQUEST_NULL;
    float* path_data = nullptr;
    propagate_kernels(tracker,path_data,&prop_req);
    while (1){
      if (request_id == INT_MAX) request_id = 100;// reset to avoid overflow. rare case.
      if ((nonblocking_internal_info.find(request_id) == nonblocking_internal_info.end()) && (request_id != MPI_REQUEST_NULL)){
        nonblocking_info msg_info(path_data,barrier_req,prop_req,false,is_sender,partner,comm,(float)nbytes,(float)np,user_tag,&tracker);
        nonblocking_internal_info[request_id] = msg_info;
        break;
      }
      request_id++;
    }
  }

  computation_timer = MPI_Wtime();
  if (skeleton_type){ return false; }
  else{ return true; }
}

// Called by both nonblocking p2p and nonblocking collectives
void path::initiate_comm(nonblocking& tracker, volatile double itime, int64_t nelem,
                         MPI_Datatype t, MPI_Comm comm, MPI_Request* request, int user_tag, bool is_sender, int partner){
  // Note: this function is invoked only when skeleton_type==0
  assert(skeleton_type == 0);
  MPI_Request barrier_req = MPI_REQUEST_NULL;// Only necessary for nonblocking receives
  cp_costs[0] += itime;
  vol_costs[0] += itime;
  int rank; MPI_Comm_rank(comm,&rank); 
  int word_size,np; MPI_Type_size(t, &word_size);
  int64_t nbytes = word_size * nelem;
  MPI_Comm_size(comm, &np);
  if (rank == partner){
    computation_timer = MPI_Wtime();
    return;
  }
  // Issue the barrier call, regardless of msg size
  // Note that this is only necessary due to blocking+nonblocking p2p communication
  // Therefore, nonblocking collectives need not participate in the barrier call
  if (partner!=-1 && rank != partner){// Branch protects against nonblocking collectives
    if (is_sender){
      MPI_Buffer_attach(&eager_pad[0],eager_pad.size());
      PMPI_Bsend(&barrier_pad_send, 1, MPI_CHAR, partner, internal_tag1, comm);
      void* temp_buf; int temp_size;
      MPI_Buffer_detach(&temp_buf,&temp_size);
    }
    else{
      PMPI_Irecv(&barrier_pad_recv, 1, MPI_CHAR, partner, internal_tag1, comm, &barrier_req);
    }
  } else{
    PMPI_Iallreduce(&barrier_pad_send, &barrier_pad_recv, 1, MPI_CHAR, MPI_MAX, comm, &barrier_req);
  }
  MPI_Request prop_req = MPI_REQUEST_NULL;
  float* path_data = nullptr;
  propagate_kernels(tracker,path_data,&prop_req);
  nonblocking_info msg_info(path_data,barrier_req,prop_req,true,is_sender,partner,comm,(float)nbytes,(float)np,user_tag,&tracker);
  nonblocking_internal_info[*request] = msg_info;
  computation_timer = MPI_Wtime();
}

void path::complete_comm(nonblocking& tracker, MPI_Request* request, double comp_time, double comm_time){
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
  std::pair<float,float> cost_alphabeta = tracker.cost_func_alphabeta(tracker.nbytes, tracker.comm_size);

  int comm_sizes[2]={0,0}; int comm_strides[2]={0,0};
  assert(comm_channel_map.find(tracker.comm) != comm_channel_map.end());
  for (auto i=0; i<comm_channel_map[tracker.comm]->id.size(); i++){
    comm_sizes[i]=comm_channel_map[tracker.comm]->id[i].first;
    comm_strides[i]=comm_channel_map[tracker.comm]->id[i].second;
  }
  // Below, the idea is that key doesn't exist in comm_kernel_map iff the key hasn't been seen before. If the key has been seen, we automatically
  //   create an entry in comm_kernel_key, although it will be empty.
  comm_kernel_key key(rank,active_kernels.size(),tracker.tag,comm_sizes,comm_strides,tracker.nbytes,tracker.partner1);
  if (comm_kernel_map.find(key) == comm_kernel_map.end()){
    active_comm_kernel_keys.push_back(key);
    active_kernels.push_back(1);
    comm_kernel_map[key] = kernel_key_id(true,active_comm_kernel_keys.size()-1,active_kernels.size()-1,false);
  }
  else{
    active_kernels[comm_kernel_map[key].val_index]++;
  }

  if (skeleton_type==0){
    cp_costs[0] += comp_time+comm_time;
    vol_costs[0] += comp_time+comm_time;
    vol_costs[0] = vol_costs[0] > cp_costs[0] ? cp_costs[0] : vol_costs[0];
  } else{
    cp_costs[0] += cost_alphabeta.first*1e6 + cost_alphabeta.second*1e3;
    vol_costs[0] += cost_alphabeta.first*1e6 + cost_alphabeta.second*1e3;
  }
  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
  nonblocking_internal_info.erase(*request);
}

int path::complete_comm(double curtime, MPI_Request* request, MPI_Status* status){
  auto comp_time = curtime - computation_timer;
  int ret = MPI_SUCCESS;
  assert(nonblocking_internal_info.find(*request) != nonblocking_internal_info.end());
  auto info_it = nonblocking_internal_info.find(*request);
  auto save_r = info_it->first;
  // If the request is active, that must mean that skeleton_type==0
  if (info_it->second.is_active == 1){
    volatile auto last_start_time = MPI_Wtime();
    ret = PMPI_Wait(request, status);
    auto save_comm_time = MPI_Wtime() - last_start_time;
    int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
    if (rank != info_it->second.partner){
      // If receiver or collective, complete the barrier and the path data propagation
      if (!info_it->second.is_sender || info_it->second.partner==-1){
        assert(info_it->second.path_data != nullptr);
        MPI_Request req_array[] = {info_it->second.barrier_req, info_it->second.prop_req};
        PMPI_Waitall(2, &req_array[0], MPI_STATUSES_IGNORE);
        if (info_it->second.partner != -1) update_frequency(info_it->second.path_data,&cp_costs[0],cp_costs_size);
        kernel_update(&info_it->second.path_data[1]);
        free(info_it->second.path_data);
      }
      complete_comm(*info_it->second.track, &save_r, comp_time, save_comm_time);
    }
  } else{
    int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
    if (rank != info_it->second.partner){
      // If receiver or collective, complete the barrier and the path data propagation
      if (!info_it->second.is_sender || info_it->second.partner==-1){
        assert(info_it->second.path_data != nullptr);
        MPI_Request req_array[] = {info_it->second.barrier_req, info_it->second.prop_req};
        PMPI_Waitall(2, &req_array[0], MPI_STATUSES_IGNORE);
        if (info_it->second.partner != -1) update_frequency(info_it->second.path_data,&cp_costs[0],cp_costs_size);
        kernel_update(&info_it->second.path_data[1]);
        free(info_it->second.path_data);
        if (status != MPI_STATUS_IGNORE){
          status->MPI_SOURCE = info_it->second.partner;
          status->MPI_TAG = info_it->second.tag;
        }
      }
      complete_comm(*info_it->second.track, &save_r, comp_time, (float)1000000.);
    }
    *request = MPI_REQUEST_NULL;
  }
  computation_timer = MPI_Wtime();
  return ret;
}

int path::complete_comm(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){
  auto comp_time = curtime - computation_timer;
  int ret = MPI_SUCCESS;
  if (skeleton_type == 0){
    std::vector<MPI_Request> pt(count); for (int i=0;i<count;i++){pt[i]=(array_of_requests)[i];}
    volatile auto last_start_time = MPI_Wtime();
    ret = PMPI_Waitany(count,array_of_requests,indx,status);
    auto waitany_comm_time = MPI_Wtime() - last_start_time;
    MPI_Request request = pt[*indx];
    auto info_it = nonblocking_internal_info.find(request);
    assert(info_it != nonblocking_internal_info.end());
    int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
    if (rank != info_it->second.partner){
      // If receiver or collective, complete the barrier and the path data propagation
      if (!info_it->second.is_sender || info_it->second.partner==-1){
        assert(info_it->second.path_data != nullptr);
        MPI_Request req_array[] = {info_it->second.barrier_req, info_it->second.prop_req};
        PMPI_Waitall(2, &req_array[0], MPI_STATUSES_IGNORE);
        if (info_it->second.partner != -1) update_frequency(info_it->second.path_data,&cp_costs[0],cp_costs_size);
        kernel_update(&info_it->second.path_data[1]);
        free(info_it->second.path_data);
      }
      complete_comm(*info_it->second.track, &request, comp_time, waitany_comm_time);
    }
  } else{
    int save_request_index=-1;
    // Find the last open fake request
    for (auto i=0; i<count; i++) if (array_of_requests[i] != MPI_REQUEST_NULL) { save_request_index = i; }
    assert(save_request_index >= 0);
    auto info_it = nonblocking_internal_info.find(array_of_requests[save_request_index]);
    assert(info_it != nonblocking_internal_info.end());
    int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
    if (rank != info_it->second.partner){
      // If receiver or collective, complete the barrier and the path data propagation
      if (!info_it->second.is_sender || info_it->second.partner==-1){
        assert(info_it->second.path_data != nullptr);
        MPI_Request req_array[] = {info_it->second.barrier_req, info_it->second.prop_req};
        PMPI_Waitall(2, &req_array[0], MPI_STATUSES_IGNORE);
        if (info_it->second.partner != -1) update_frequency(info_it->second.path_data,&cp_costs[0],cp_costs_size);
        kernel_update(&info_it->second.path_data[1]);
        free(info_it->second.path_data);
        if (status != MPI_STATUS_IGNORE){
          status->MPI_SOURCE = info_it->second.partner;
          status->MPI_TAG = info_it->second.tag;
        }
      }
      complete_comm(*info_it->second.track, &array_of_requests[save_request_index], comp_time, (float)1000000.);
    }
    array_of_requests[save_request_index] = MPI_REQUEST_NULL;
  }
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
  computation_timer = MPI_Wtime();
  return ret;
}

int path::complete_comm(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
  auto comp_time = curtime - computation_timer;
  int ret = MPI_SUCCESS;
  is_first_request=true;
  if (skeleton_type==0){
    std::vector<MPI_Request> pt(count); for (int i=0;i<count;i++){pt[i]=(array_of_requests)[i];}
    volatile auto last_start_time = MPI_Wtime();
    ret = PMPI_Waitall(count,array_of_requests,array_of_statuses);
    auto waitall_comm_time = MPI_Wtime() - last_start_time;
    for (int i=0; i<count; i++){
      auto info_it = nonblocking_internal_info.find(pt[i]);
      assert(info_it != nonblocking_internal_info.end());
      int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
      if (rank != info_it->second.partner){
        // If receiver or collective, complete the barrier and the path data propagation
        if (!info_it->second.is_sender || info_it->second.partner==-1){
          assert(info_it->second.path_data != nullptr);
          MPI_Request req_array[] = {info_it->second.barrier_req, info_it->second.prop_req};
          PMPI_Waitall(2, &req_array[0], MPI_STATUSES_IGNORE);
          if (info_it->second.partner != -1) update_frequency(info_it->second.path_data,&cp_costs[0],cp_costs_size);
          kernel_update(&info_it->second.path_data[1]);
          free(info_it->second.path_data);
        }
        complete_comm(*info_it->second.track, &pt[i], comp_time, waitall_comm_time);
        // Although we have to exchange the path data for each request, we do not want to float-count the computation time nor the communicaion time
        comp_time=0; waitall_comm_time=0; is_first_request=false;
      }
    }
  } else{
    for (int i=0; i<count; i++){
      auto info_it = nonblocking_internal_info.find(array_of_requests[i]);
      assert(info_it != nonblocking_internal_info.end());
      int rank; MPI_Comm_rank(info_it->second.track->comm,&rank);
      if (rank != info_it->second.partner){
        // If receiver or collective, complete the barrier and the path data propagation
        if (!info_it->second.is_sender || info_it->second.partner==-1){
          assert(info_it->second.path_data != nullptr);
          MPI_Request req_array[] = {info_it->second.barrier_req, info_it->second.prop_req};
          PMPI_Waitall(2, &req_array[0], MPI_STATUSES_IGNORE);
          update_frequency(info_it->second.path_data,&cp_costs[0],cp_costs_size);
          kernel_update(&info_it->second.path_data[1]);
          free(info_it->second.path_data);
          if (array_of_statuses != MPI_STATUSES_IGNORE){
            array_of_statuses[i].MPI_SOURCE = info_it->second.partner;
            array_of_statuses[i].MPI_TAG = info_it->second.tag;
          }
        }
        complete_comm(*info_it->second.track, &array_of_requests[i], comp_time, (float)100000.);
        comp_time=0; is_first_request=false;
      }
      array_of_requests[i] = MPI_REQUEST_NULL;
    }
  }
  computation_timer = MPI_Wtime();
  return ret;
}

void path::propagate_kernels(blocking& tracker){
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  if ((rank == tracker.partner1) && (rank == tracker.partner2)) { return; } 

  // Fill in -1 first because the number of distinct kernels might be less than 'comm_kernel_select_count',
  //   just to avoid confusion. A -1 tag clearly means that the kernel is void
  memset(&cp_costs[1],-1,sizeof(float)*(cp_costs.size()-1));
  // Iterate over first 'comp_kernel_select_count' keys
  int count=0;
  for (auto it : comp_kernel_map){
    if (comp_kernel_select_count==0) break;
    auto offset = 1+9*count;
    cp_costs[offset] = it.first.tag;
    cp_costs[offset+1] = it.first.param1;
    cp_costs[offset+2] = it.first.param2;
    cp_costs[offset+3] = it.first.param3;
    cp_costs[offset+4] = it.first.param4;
    cp_costs[offset+5] = it.first.param5;
    cp_costs[offset+6] = it.first.kernel_index;
    cp_costs[offset+7] = it.first.flops;
    cp_costs[offset+8] = active_kernels[it.second.val_index];
    count++; if (count==comp_kernel_select_count) break;
  }
  count=0;
  for (auto it : comm_kernel_map){
    if (comm_kernel_select_count==0) break;
    auto offset = 1+9*comp_kernel_select_count+count*9;
    cp_costs[offset] = it.first.tag;
    cp_costs[offset+1] = it.first.dim_sizes[0];
    cp_costs[offset+2] = it.first.dim_sizes[1];
    cp_costs[offset+3] = it.first.dim_strides[0];
    cp_costs[offset+4] = it.first.dim_strides[1];
    cp_costs[offset+5] = it.first.partner_offset;
    cp_costs[offset+6] = it.first.kernel_index;
    cp_costs[offset+7] = it.first.msg_size;
    cp_costs[offset+8] = active_kernels[it.second.val_index];
    count++; if (count==comm_kernel_select_count) break;
  }

  // Exchange the tracked routine critical path data
  if (tracker.partner1 == -1){
    MPI_Op op;
    MPI_Op_create((MPI_User_function*) propagate_cp_op,0,&op);
    PMPI_Allreduce(&cp_costs[0], &cp_costs_foreign[0], cp_costs.size(), MPI_FLOAT, op, tracker.comm);
    MPI_Op_free(&op);
  } else{
    if (tracker.is_sender){
      MPI_Buffer_attach(&eager_pad[0],eager_pad.size());
      PMPI_Bsend(&cp_costs[0], cp_costs.size(), MPI_FLOAT, tracker.partner1, internal_tag2, tracker.comm);
      void* temp_buf; int temp_size;
      MPI_Buffer_detach(&temp_buf,&temp_size);
    } else{
      PMPI_Recv(&cp_costs_foreign[0], cp_costs_foreign.size(), MPI_FLOAT, tracker.partner1, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
      update_frequency(&cp_costs_foreign[0],&cp_costs[0],cp_costs_size);
    }
    if (tracker.partner2 != tracker.partner1 || tracker.tag==13 || tracker.tag==14){
      PMPI_Recv(&cp_costs_foreign[0], cp_costs_foreign.size(), MPI_FLOAT, tracker.partner2, internal_tag2, tracker.comm, MPI_STATUS_IGNORE);
      update_frequency(&cp_costs_foreign[0],&cp_costs[0],cp_costs_size);
    }
  }

  // Senders can leave early
  if (tracker.partner1 != -1 && tracker.is_sender && tracker.partner2 == tracker.partner1) return;
  // Receivers always update their kernel maps
  kernel_update(&cp_costs_foreign[1]);
}

void path::propagate_kernels(nonblocking& tracker, float*& path_data, MPI_Request* prop_req){
  int rank; MPI_Comm_rank(tracker.comm,&rank);
  if (rank == tracker.partner1) { return; } 

  // Fill in -1 first because the number of distinct kernels might be less than 'com*_kernel_select_count',
  //   just to avoid confusion. A -1 tag clearly means that the kernel is void
  memset(&cp_costs[1],-1,sizeof(float)*(cp_costs.size()-1));
  // Iterate over map
  int count=0;
  for (auto it : comp_kernel_map){
    if (comp_kernel_select_count==0) break;
    auto offset = 1+9*count;
    cp_costs[offset] = it.first.tag;
    cp_costs[offset+1] = it.first.param1;
    cp_costs[offset+2] = it.first.param2;
    cp_costs[offset+3] = it.first.param3;
    cp_costs[offset+4] = it.first.param4;
    cp_costs[offset+5] = it.first.param5;
    cp_costs[offset+6] = it.first.kernel_index;
    cp_costs[offset+7] = it.first.flops;
    cp_costs[offset+8] = active_kernels[it.second.val_index];
    count++; if (count==comp_kernel_select_count) break;
  }
  count=0;
  for (auto it : comm_kernel_map){
    if (comm_kernel_select_count==0) break;
    auto offset = 1+9*comp_kernel_select_count+9*count;
    cp_costs[offset] = it.first.tag;
    cp_costs[offset+1] = it.first.dim_sizes[0];
    cp_costs[offset+2] = it.first.dim_sizes[1];
    cp_costs[offset+3] = it.first.dim_strides[0];
    cp_costs[offset+4] = it.first.dim_strides[1];
    cp_costs[offset+5] = it.first.partner_offset;
    cp_costs[offset+6] = it.first.kernel_index;
    cp_costs[offset+7] = it.first.msg_size;
    cp_costs[offset+8] = active_kernels[it.second.val_index];
    count++; if (count==comm_kernel_select_count) break;
  }

  // Exchange the tracked routine critical path data
  if (tracker.partner1 == -1){
    MPI_Op op;
    MPI_Op_create((MPI_User_function*) propagate_cp_op,0,&op);
    path_data = (float*)malloc(cp_costs.size()*sizeof(float));
    std::memcpy(path_data, &cp_costs[0], cp_costs.size()*sizeof(float));
    PMPI_Iallreduce(MPI_IN_PLACE, path_data, cp_costs.size(), MPI_FLOAT, op, tracker.comm, prop_req);
    //MPI_Op_free(&op);
  } else{
    if (tracker.is_sender){
      MPI_Buffer_attach(&eager_pad[0],eager_pad.size());
      PMPI_Bsend(&cp_costs[0], cp_costs.size(), MPI_FLOAT, tracker.partner1, internal_tag2, tracker.comm);
      void* temp_buf; int temp_size;
      MPI_Buffer_detach(&temp_buf,&temp_size);
    } else{
      path_data = (float*)malloc(cp_costs.size()*sizeof(float));
      PMPI_Irecv(path_data, cp_costs.size(), MPI_FLOAT, tracker.partner1, internal_tag2, tracker.comm, prop_req);
    }
  }
}

}
}
}
