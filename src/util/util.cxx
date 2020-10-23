#include <cstring>
#include <stdint.h>
#include <limits.h>

#include "util.h"

namespace critter{
namespace internal{

// ****************************************************************************************************************************************************
static int gcd(int a, int b){
  if (a==b) return a;
  if (a>b) return gcd(a-b,b);
  else return gcd(a,b-a);
}
static int lcm(int a, int b){
  return a*b / gcd(a,b);
}

static double truncate(double val, int unit_param){
  if (unit_param==0) return val;
  // returns the next highest power of 2 as a double
  if (val==0) return 0;
  // I'm counting on the int64_t to hold sufficiently large number, especially to store the flop_count exactly.
  int64_t v = val;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v |= v >> 32;
  v++;
  if (unit_param==1){
    return (double)v;
  }
  else{
    int pos=0; int64_t temp = v;
    while (temp){
      pos++; temp>>=1;
    }
    pos--;
    if (pos%unit_param==0) { return (double)v; }
    for (int z = pos%unit_param; z<unit_param; z++){
      v<<=1;
    }
    return (double)v;
  }
}
static double truncate(int v, int unit_param){
  if (unit_param==0) return v;
  // returns the next highest power of 2 as a double
  if (v==0) return 0;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  if (unit_param==1){
    return v;
  }
  else{
    int pos=0; int64_t temp = v;
    while (temp){
      pos++; temp>>=1;
    }
    pos--;
    if (pos%unit_param==0) { return (double)v; }
    for (int z = pos%unit_param; z<unit_param; z++){
      v<<=1;
    }
    return (double)v;
  }
}

// ****************************************************************************************************************************************************
int comm_envelope_param;
int comp_envelope_param;
int comm_unit_param;
int comp_unit_param;
int comm_analysis_param;
int comp_analysis_param;
// ****************************************************************************************************************************************************
comm_kernel_key::comm_kernel_key(int _rank, int _kernel_index, int _tag, int _dim_sizes[2], int _dim_strides[2], double _msg_size, int _partner){
  this->kernel_index = _kernel_index;
  // Envelope (non-message-size) parameterization specification
  this->tag = _tag;
  if (_dim_sizes != nullptr) std::memcpy(&this->dim_sizes[0],&_dim_sizes[0],2*sizeof(int));
  if (_dim_strides != nullptr) std::memcpy(&this->dim_strides[0],&_dim_strides[0],2*sizeof(int));
  if (comm_envelope_param == 0){
    this->partner_offset = _rank-_partner;
  }
  else if (comm_envelope_param == 1){
    this->partner_offset = abs(_rank-_partner);
  }
  else if (comm_envelope_param == 2){
    this->partner_offset = _partner;
  }
  else if (comm_envelope_param == 3){// Note: this parameter specifies that partner rank will not factor into parameterization (but we still need to give the member a value)
    this->partner_offset = -1;
  }
  // Unit (message-size) parameterization specification
  this->msg_size = truncate(_msg_size,comm_unit_param);
  // Regardless of the specified envelope parameterization, non-p2p communication requires all processes to set partner_offset <- INT_MIN (-1 is not sufficient)
  if (_partner == -1){ this->partner_offset = INT_MIN; }
}

// This constructor is used when transferring ownership of kernels following path propagation.
comm_kernel_key::comm_kernel_key(int _kernel_index, int _tag, int _dim_sizes[2], int _dim_strides[2], double _msg_size, int _partner_offset){
  this->kernel_index = _kernel_index;
  this->tag = _tag;
  if (_dim_sizes != nullptr) std::memcpy(&this->dim_sizes[0],&_dim_sizes[0],2*sizeof(int));
  if (_dim_strides != nullptr) std::memcpy(&this->dim_strides[0],&_dim_strides[0],2*sizeof(int));
  this->partner_offset = _partner_offset;
  this->msg_size = _msg_size;
}

comm_kernel_key::comm_kernel_key(const comm_kernel_key& _copy){
  this->kernel_index = _copy.kernel_index;
  this->tag = _copy.tag;
  if (_copy.dim_sizes != nullptr) std::memcpy(&this->dim_sizes[0],&_copy.dim_sizes[0],2*sizeof(int));
  if (_copy.dim_strides != nullptr) std::memcpy(&this->dim_strides[0],&_copy.dim_strides[0],2*sizeof(int));
  this->msg_size = _copy.msg_size;
  this->partner_offset = _copy.partner_offset;
}

comm_kernel_key& comm_kernel_key::operator=(const comm_kernel_key& _copy){
  this->kernel_index = _copy.kernel_index;
  this->tag = _copy.tag;
  if (_copy.dim_sizes != nullptr) std::memcpy(&this->dim_sizes[0],&_copy.dim_sizes[0],2*sizeof(int));
  if (_copy.dim_strides != nullptr) std::memcpy(&this->dim_strides[0],&_copy.dim_strides[0],2*sizeof(int));
  this->msg_size = _copy.msg_size;
  this->partner_offset = _copy.partner_offset;
  return *this;
}

bool operator==(const comm_kernel_key& ref1, const comm_kernel_key& ref2){
  // Note that because of how we set the member variables in the constructor based on envlope, unit, and analysis parameterizations, no branching is required here.
  if ((ref1.tag==ref2.tag) &&
      (ref1.dim_sizes[0] == ref2.dim_sizes[0]) && (ref1.dim_sizes[1] == ref2.dim_sizes[1]) && (ref1.dim_sizes[2] == ref2.dim_sizes[2]) &&
      (ref1.dim_strides[0] == ref2.dim_strides[0]) && (ref1.dim_strides[1] == ref2.dim_strides[1]) && (ref1.dim_strides[2] == ref2.dim_strides[2]) &&
      (ref1.msg_size == ref2.msg_size) && (ref1.partner_offset == ref2.partner_offset)) return true;
  else return false;
}

bool operator<(const comm_kernel_key& ref1, const comm_kernel_key& ref2){
  if (ref1.tag < ref2.tag) return true;
  else if (ref1.tag > ref2.tag) return false;
  if (ref1.dim_sizes[0] < ref2.dim_sizes[0]) return true;
  else if (ref1.dim_sizes[0] > ref2.dim_sizes[0]) return false;
  if (ref1.dim_sizes[1] < ref2.dim_sizes[1]) return true;
  else if (ref1.dim_sizes[1] > ref2.dim_sizes[1]) return false;
  if (ref1.dim_sizes[2] < ref2.dim_sizes[2]) return true;
  else if (ref1.dim_sizes[2] > ref2.dim_sizes[2]) return false;
  if (ref1.dim_strides[0] < ref2.dim_strides[0]) return true;
  else if (ref1.dim_strides[0] > ref2.dim_strides[0]) return false;
  if (ref1.dim_strides[1] < ref2.dim_strides[1]) return true;
  else if (ref1.dim_strides[1] > ref2.dim_strides[1]) return false;
  if (ref1.dim_strides[2] < ref2.dim_strides[2]) return true;
  else if (ref1.dim_strides[2] > ref2.dim_strides[2]) return false;
  if (ref1.msg_size < ref2.msg_size) return true;
  else if (ref1.msg_size > ref2.msg_size) return false;
  if (ref1.partner_offset < ref2.partner_offset) return true;
  else if (ref1.partner_offset > ref2.partner_offset) return false;
  return false;
}

// ****************************************************************************************************************************************************
comp_kernel_key::comp_kernel_key(int _kernel_index, int _tag, double _flops, int _param1, int _param2, int _param3, int _param4, int _param5){
  this->kernel_index = _kernel_index;
  this->tag = _tag;
  this->param1 = truncate(_param1,comp_unit_param);
  this->param2 = truncate(_param2,comp_unit_param);
  this->param3 = truncate(_param3,comp_unit_param);
  this->param4 = truncate(_param4,comp_unit_param);
  this->param5 = truncate(_param5,comp_unit_param);
  this->flops = truncate(_flops,comp_unit_param);
}

comp_kernel_key::comp_kernel_key(const comp_kernel_key& _copy){
  this->kernel_index = _copy.kernel_index;
  this->tag = _copy.tag;
  this->flops = _copy.flops;
  this->param1 = _copy.param1;
  this->param2 = _copy.param2;
  this->param3 = _copy.param3;
  this->param4 = _copy.param4;
  this->param5 = _copy.param5;
}

comp_kernel_key& comp_kernel_key::operator=(const comp_kernel_key& _copy){
  this->kernel_index = _copy.kernel_index;
  this->tag = _copy.tag;
  this->flops = _copy.flops;
  this->param1 = _copy.param1;
  this->param2 = _copy.param2;
  this->param3 = _copy.param3;
  this->param4 = _copy.param4;
  this->param5 = _copy.param5;
  return *this;
}

bool operator==(const comp_kernel_key& ref1, const comp_kernel_key& ref2){
  if ((ref1.tag==ref2.tag) && (ref1.param1 == ref2.param1) && (ref1.param2 == ref2.param2) && (ref1.param3 == ref2.param3) && (ref1.param4 == ref2.param4) && (ref1.param5 == ref2.param5)) return true;
  else return false;
}

bool operator<(const comp_kernel_key& ref1, const comp_kernel_key& ref2){
  if (ref1.tag < ref2.tag) return true;
  else if (ref1.tag > ref2.tag) return false;
  if (ref1.param1 < ref2.param1) return true;
  else if (ref1.param1 > ref2.param1) return false;
  if (ref1.param2 < ref2.param2) return true;
  else if (ref1.param2 > ref2.param2) return false;
  if (ref1.param3 < ref2.param3) return true;
  else if (ref1.param3 > ref2.param3) return false;
  if (ref1.param4 < ref2.param4) return true;
  else if (ref1.param4 > ref2.param4) return false;
  if (ref1.param5 < ref2.param5) return true;
  else if (ref1.param5 > ref2.param5) return false;
  return false;
}

// ****************************************************************************************************************************************************
kernel_key_id::kernel_key_id(bool _is_active, int _key_index, int _val_index, bool _is_updated){
  this->is_active = _is_active;
  this->key_index = _key_index;
  this->val_index = _val_index;
  this->is_updated = _is_updated;
}

kernel_key_id::kernel_key_id(const kernel_key_id& _copy){
  this->is_active = _copy.is_active;
  this->key_index = _copy.key_index;
  this->val_index = _copy.val_index;
  this->is_updated = _copy.is_updated;
}

kernel_key_id& kernel_key_id::operator=(const kernel_key_id& _copy){
  this->is_active = _copy.is_active;
  this->key_index = _copy.key_index;
  this->val_index = _copy.val_index;
  this->is_updated = _copy.is_updated;
  return *this;
}


int bsp_counter;
int reset_counter;
int clear_counter;
int communicator_count;
std::string schedule_tag;
volatile double computation_timer;
std::vector<double> wall_timer;
double _wall_time;
size_t auto_capture;
bool is_world_root;
size_t mechanism,mode,stack_id;
double scratch_pad;
size_t track_blas;
size_t track_lapack;
size_t track_collective;
size_t track_p2p;
size_t track_p2p_idle;
size_t eager_p2p;
size_t delete_comm;
size_t
         _MPI_Send__id,
         _MPI_Ssend__id,
         _MPI_Bsend__id,
         _MPI_Recv__id,
         _MPI_Sendrecv__id,
         _MPI_Sendrecv_replace__id,
         _MPI_Barrier__id,
         _MPI_Bcast__id,
         _MPI_Reduce__id,
         _MPI_Allreduce__id,
         _MPI_Gather__id,
         _MPI_Allgather__id,
         _MPI_Scatter__id,
         _MPI_Reduce_scatter__id,
         _MPI_Alltoall__id,
         _MPI_Gatherv__id,
         _MPI_Allgatherv__id,
         _MPI_Scatterv__id,
         _MPI_Alltoallv__id,
         _MPI_Isend__id,
         _MPI_Irecv__id,
         _MPI_Ibcast__id,
         _MPI_Iallreduce__id,
         _MPI_Ireduce__id,
         _MPI_Igather__id,
         _MPI_Igatherv__id,
         _MPI_Iallgather__id,
         _MPI_Iallgatherv__id,
         _MPI_Iscatter__id,
         _MPI_Iscatterv__id,
         _MPI_Ireduce_scatter__id,
         _MPI_Ialltoall__id,
         _MPI_Ialltoallv__id;
size_t
	_BLAS_axpy__id,
	_BLAS_scal__id,
	_BLAS_ger__id,
	_BLAS_gemm__id,
	_BLAS_trmm__id,
	_BLAS_trsm__id,
	_BLAS_syrk__id;
size_t
	_LAPACK_getrf__id,
	_LAPACK_potrf__id,
	_LAPACK_trtri__id,
	_LAPACK_geqrf__id,
	_LAPACK_orgqr__id,
	_LAPACK_ormqr__id,
	_LAPACK_getri__id,
	_LAPACK_tpqrt__id,
	_LAPACK_tpmqrt__id;
size_t
	_CAPITAL_blktocyc__id;
//std::map<std::pair<std::string,size_t>,bool> schedule_map;

std::map<MPI_Comm,solo_channel*> comm_channel_map;
std::map<int,solo_channel*> p2p_channel_map;
std::map<int,aggregate_channel*> aggregate_channel_map;

// ****************************************************************************************************************************************************
channel::channel(){
  this->offset = INT_MIN;	// will be updated later
}

std::vector<std::pair<int,int>> channel::generate_tuple(std::vector<int>& ranks, int new_comm_size){
  std::vector<std::pair<int,int>> tuple_list;
  if (new_comm_size<=1){
    tuple_list.emplace_back(new_comm_size,1);
  }
  else{
    int stride = ranks[1]-ranks[0];
    int count = 0;
    int jump=1;
    int extra=0;
    int i=0;
    while (i < ranks.size()-1){
      if ((ranks[i+jump]-ranks[i]) != stride){
        tuple_list.emplace_back(count+extra+1,stride);// node->id.push_back(std::make_pair(count+extra+1,stride));
        stride = ranks[i+1]-ranks[0];
        i += jump;
        if (tuple_list.size()==1){
          jump=count+extra+1;
        }
        else{
          jump = (count+extra+1)*tuple_list[tuple_list.size()-2].first;
        }
        extra=1;
        count = 0;
      } else{
        count++;
        i += jump;
      }
    }
    if (count != 0){
      tuple_list.emplace_back(count+extra+1,stride);//node->id.push_back(std::make_pair(count+extra+1,stride));
    }
  }
  assert(tuple_list.size() >= 1);
  assert(tuple_list.size() <= 2);
  return tuple_list;
}
void channel::contract_tuple(std::vector<std::pair<int,int>>& tuple_list){
  int index=0;
  for (int i=1; i<tuple_list.size(); i++){
    if (tuple_list[index].first*tuple_list[index].second == tuple_list[i].second){
      tuple_list[index].first *= tuple_list[i].first;
      tuple_list[index].second = std::min(tuple_list[index].second,tuple_list[i].second);
    }
    else{
      index++;
    }
  }
  int remove_len = tuple_list.size() - index - 1;
  for (auto i=0; i<remove_len; i++){ tuple_list.pop_back(); }
}
int channel::enumerate_tuple(channel* node, std::vector<int>& process_list){
  int count = 0;
  if (node->id.size()==1){
    int offset = node->offset;
    for (auto i=0; i<node->id[0].first; i++){
      process_list.push_back(offset + i*node->id[0].second);
      count++;
      if (node->id[0].second==0) break;// this might occur if a p2p send/receives with itself
    }
  } else{
    assert(node->id.size()==2);
    bool op = (node->id[0].first*node->id[0].second) < node->id[0].second;
    int max_process;
    if (op){
      max_process = node->offset + channel::span(node->id[0]) + channel::span(node->id[1]);
    } else{
      max_process = node->offset + std::max(channel::span(node->id[0]), channel::span(node->id[1]));
    }
    int offset = node->offset;
    for (auto i=0; i<node->id[1].first; i++){
      for (auto j=0; j<node->id[0].first; j++){
        if (offset + i*node->id[0].second <= max_process){ process_list.push_back(offset + i*node->id[0].second); count++; }
      }
      offset += node->id[1].second;
    }
  }
  return count;
}
int channel::duplicate_process_count(std::vector<int>& process_list){
  int count=0;
  int save = process_list[0];
  for (int i=0; i<process_list.size(); i++){
    if (process_list[i] == save){
      count++;
    } else{
      save = process_list[i];// restart duplicate tracking
    }
  }
  return count;
}
int channel::translate_rank(MPI_Comm comm, int rank){
  // Returns new rank relative to world communicator
  auto node = comm_channel_map[comm];
  int new_rank = node->offset;
  for (auto i=0; i<node->id.size(); i++){
    new_rank += node->id[i].second*(rank%node->id[i].first);
    rank /= node->id[i].first;
  }
  return new_rank;
}
std::string channel::generate_tuple_string(channel* comm){
  std::string str1 = "{ offset = " + std::to_string(comm->offset) + ", ";
  for (auto it : comm->id){
    str1 += " (" + std::to_string(it.first) + "," + std::to_string(it.second) + ")";
  }
  str1+=" }";
  return str1;
}

bool channel::verify_ancestor_relation(channel* comm1, channel* comm2){
/*
  // First check that the parent 'comm2' is not a p2p, regardless of whether 'tree_node' is a subcomm or p2p
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  if (comm2->tag >= ((-1)*world_size)) return false;
  //TODO: Commented out is old support for checking if a p2p is a child of a channel. Not sure if still correct.
  // p2p nodes need special machinery for detecting whether its a child of 'node'
  int rank = tree_node->offset - node->offset;//TODO; Fix if tree_node describes a aggregated p2p (rare case, we haven't needed it yet)
  if (rank<0) return false;// corner case
  for (int i=node->id.size()-1; i>=0; i--){
    if (rank >= (node->id[i].first*node->id[i].second)){ return false; }
    if (i>0) { rank %= node->id[i].second; }
  }
  if (rank%node->id[0].second != 0) return false;
  return true;
*/
  if ((comm1->id.size() == 1) && (comm2->id.size() == 1)){
    return ((comm1->offset >= comm2->offset) &&
            (comm1->offset+channel::span(comm1->id[0]) <= comm2->offset+channel::span(comm2->id[0])) &&
            ((comm1->offset - comm2->offset) % comm2->id[0].second == 0) &&
            ((comm1->id[0].second % comm2->id[0].second) == 0));
  }
  // Spill to process list, sort, identify if |comm1| duplicates.
  std::vector<int> process_list;
  int count1 = channel::enumerate_tuple(comm1,process_list);
  int count2 = channel::enumerate_tuple(comm2,process_list);
  std::sort(process_list.begin(),process_list.end());
  assert(process_list.size()>0);
  int count = channel::duplicate_process_count(process_list);
  return (count==count1);
}
bool channel::verify_sibling_relation(channel* comm1, channel* comm2){
  if ((comm1->id.size() == 1) && (comm2->id.size() == 1)){
    int min1 = comm1->offset + span(comm1->id[0]);
    int min2 = comm2->offset + span(comm2->id[0]);
    int _min_ = std::min(min1,min2);
    int max1 = comm1->offset;
    int max2 = comm2->offset;
    int _max_ = std::max(max1,max2);
     // Two special cases -- if stride of one channel is 0, that means that a single intersection point is guaranteed.
    if (comm2->id[0].second == 0 || comm1->id[0].second == 0){ return true; }
    int _lcm_ = lcm(comm1->id[0].second,comm2->id[0].second);
    return _lcm_ > (_min_ - _max_);
  }
  // Spill to process list, sort, identify if 1 duplicate.
  // This is in place of a more specialized routine that iterates over the tuples to identify the pairs that identify as siblings.
  std::vector<int> process_list;
  int count1 = channel::enumerate_tuple(comm1,process_list);
  int count2 = channel::enumerate_tuple(comm2,process_list);
  std::sort(process_list.begin(),process_list.end());
  assert(process_list.size()>0);
  int count = channel::duplicate_process_count(process_list);
  return (count==1);
}
int channel::span(std::pair<int,int>& id){
  return (id.first-1)*id.second;
}

solo_channel::solo_channel(){
  this->tag = 0;	// This member MUST be overwritten immediately after instantiation
  this->frequency=0;
  this->children.push_back(std::vector<solo_channel*>());
}
bool solo_channel::verify_sibling_relation(solo_channel* node, int subtree_idx, std::vector<int>& skip_indices){
  // Perform recursive permutation generation to identify if a permutation of tuples among siblings is valid
  // Return true if node's children are valid siblings
  if (node->children[subtree_idx].size()<=1) return true;
  // Check if all are p2p, all are subcomms, or if there is a mixture.
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int p2p_count=0; int subcomm_count=0;
  for (auto i=0; i<node->children[subtree_idx].size(); i++){
    if (node->children[subtree_idx][i]->tag >= ((-1)*world_size)){ p2p_count++; }
    else { subcomm_count++; }
  }
  if ((subcomm_count>0) && (p2p_count>0)) return false;
  if (subcomm_count==0) return true;// all p2p siblings are fine
/*
  std::vector<std::pair<int,int>> static_info;
  int skip_index=0;
  for (auto i=0; i<node->children[subtree_idx].size(); i++){
    if ((skip_index<skip_indices.size()) && (i==skip_indices[skip_index])){
      skip_index++;
      continue;
    }
    if (node->children[subtree_idx][i]->tag >= ((-1)*world_size)){
      continue;
    }
    for (auto j=0; j<node->children[subtree_idx][i]->id.size(); j++){
      static_info.push_back(node->children[subtree_idx][i]->id[j]);
    }
  }
  std::vector<std::pair<int,int>> gen_info;
  std::vector<std::pair<int,int>> save_info;
  bool valid_siblings=false;
  generate_sibling_perm(static_info,gen_info,save_info,0,valid_siblings);
  if (valid_siblings){
    // Now we must run through the p2p, if any. I don't think we need to check if a p2p is also skipped via 'skip_indices', because if its a child of one of the communicators,
    //   it'll b a child of the generated span
    if (p2p_count>0){
      auto* span_node = new solo_channel;
      generate_span(span_node,save_info);
      for (auto i=0; i<node->children[subtree_idx].size(); i++){
        if (node->children[subtree_idx][i]->tag >= ((-1)*world_size)){
          valid_siblings = channel::verify_ancestor_relation(node->children[subtree_idx][i],span_node);
          if (!valid_siblings){
            delete span_node;
            return false;
          }
        }
      }
      delete span_node;
    }
  }
  return valid_siblings;
*/
  // Compare last node in 'node->children[subtree_idx]' against those that are not skipped
  auto& potential_sibling = node->children[subtree_idx][node->children[subtree_idx].size()-1];
  assert(potential_sibling->id.size() == 1);
  int skip_index=0;
  for (auto i=0; i<node->children[subtree_idx].size()-1; i++){
    if ((skip_index < skip_indices.size()) && (i == skip_indices[skip_index])){
      skip_index++;
      continue;
    }
    assert(node->children[subtree_idx][i]->id.size() == 1);
    int min1 = potential_sibling->offset + span(potential_sibling->id[0]);
    int min2 = node->children[subtree_idx][i]->offset + span(node->children[subtree_idx][i]->id[0]);
    int _min_ = std::min(min1,min2);
    int max1 = potential_sibling->offset;
    int max2 = node->children[subtree_idx][i]->offset;
    int _max_ = std::max(max1,max2);
    int _lcm_ = lcm(potential_sibling->id[0].second,node->children[subtree_idx][i]->id[0].second);
    if (_lcm_ < (_min_-_max_)) return false;
    else{// corner case check
      int count=0;
      if (_max_ == max1){
        if ((_max_ - max2) % node->children[subtree_idx][i]->id[0].second == 0) count++;
      } else{
        if ((_max_ - max1) % potential_sibling->id[0].second == 0) count++;
      }
      if (_min_ == min1){
        if ((min2 - _min_) % node->children[subtree_idx][i]->id[0].second == 0) count++;
      } else{
        if ((min1 - _min_) % potential_sibling->id[0].second == 0) count++;
      }
      // If count == 2, then we have a situation in which the endpoints match, and we must have the lcm be >= _min_-_max_, but > _min_-_max_
      if (count==2){
        if (_lcm_ == (_min_-_max_)) return false;
      }
    }
  }
  return true;
}
aggregate_channel::aggregate_channel(std::vector<std::pair<int,int>>& tuple_list, int local_hash, int global_hash, int offset, int channel_size){
  this->local_hash_tag = local_hash;
  this->global_hash_tag = global_hash;
  this->is_final = false;
  this->num_channels=channel_size;
  this->offset = offset;
  this->id = tuple_list;
}
std::string aggregate_channel::generate_hash_history(aggregate_channel* comm){
  std::string str1 = "{ hashes = ";
  int count=0;
  for (auto it : comm->channels){
    if (count>0) str1 += ",";
    str1 += std::to_string(it);
    count++;
  }
  str1+=" }";
  return str1;
}


void generate_aggregate_channels(MPI_Comm oldcomm, MPI_Comm newcomm){
  if (comm_channel_map.find(newcomm) != comm_channel_map.end()){ return; }

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
  //spf.insert_node(node);// This call will just fill in SPT via node's parent/children members, and the members of related channels
  comm_channel_map[newcomm] = node;

  // Recursively build up other legal aggregate channels that include 'node'
  std::vector<int> local_hash_array;
  std::vector<aggregate_channel*> new_aggregate_channels;
  int max_sibling_node_size=0;
  std::vector<int> save_max_indices;
  std::set<int> previous_channels_set;
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
        previous_channels_set.insert(it_2);
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
    // assert(aggregate_channel_map.find(new_aggregate_channels[i]->local_hash_tag) == aggregate_channel_map.end());
    if (aggregate_channel_map.find(new_aggregate_channels[i]->local_hash_tag) == aggregate_channel_map.end()){
      aggregate_channel_map[new_aggregate_channels[i]->local_hash_tag] = new_aggregate_channels[i];
      if (world_comm_rank == 8){
        auto str1 = channel::generate_tuple_string(new_aggregate_channels[i]);
        auto str2 = aggregate_channel::generate_hash_history(new_aggregate_channels[i]);
        std::cout << "Process " << world_comm_rank << " has aggregate " << str1 << " " << str2 << " with hashes (" << new_aggregate_channels[i]->local_hash_tag << " " << new_aggregate_channels[i]->global_hash_tag << "), num_channels - " << new_aggregate_channels[i]->num_channels << std::endl;
      }
    }
  }

  // Verify that the aggregates are build with the same hashes
  int local_sibling_size = local_hash_array.size();
  // Always treat 1-communicator channels as trivial aggregate channels.
  aggregate_channel* agg_node = new aggregate_channel(node->id,node->local_hash_tag,node->global_hash_tag,node->offset,1);
  agg_node->channels.insert(node->local_hash_tag);
  // assert(aggregate_channel_map.find(node->local_hash_tag) == aggregate_channel_map.end());
  if (aggregate_channel_map.find(node->local_hash_tag) == aggregate_channel_map.end()){
    aggregate_channel_map[node->local_hash_tag] = agg_node;
    if (world_comm_rank == 8){
      auto str1 = channel::generate_tuple_string(agg_node);
      auto str2 = aggregate_channel::generate_hash_history(agg_node);
      std::cout << "Process " << world_comm_rank << " has aggregate " << str1 << " " << str2 << " with hashes (" << agg_node->local_hash_tag << " " << agg_node->global_hash_tag << "), num_channels - " << agg_node->num_channels << std::endl;
    }
    if (local_sibling_size==0){// Only if 'node' exists as the smallest trivial aggregate should it be considered final. Think of 'node==world' of the very first registered channel
      aggregate_channel_map[node->local_hash_tag]->is_final=true;
    } else{
    }
  }
}

}
}
