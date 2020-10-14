#include <cstring>
#include <stdint.h>
#include <limits.h>

#include "util.h"

namespace critter{
namespace internal{

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
comm_pattern_key::comm_pattern_key(int _rank, int _pattern_index, int _tag, int _dim_sizes[2], int _dim_strides[2], double _msg_size, int _partner){
  this->pattern_index = _pattern_index;
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
comm_pattern_key::comm_pattern_key(int _pattern_index, int _tag, int _dim_sizes[2], int _dim_strides[2], double _msg_size, int _partner_offset){
  this->pattern_index = _pattern_index;
  this->tag = _tag;
  if (_dim_sizes != nullptr) std::memcpy(&this->dim_sizes[0],&_dim_sizes[0],2*sizeof(int));
  if (_dim_strides != nullptr) std::memcpy(&this->dim_strides[0],&_dim_strides[0],2*sizeof(int));
  this->partner_offset = _partner_offset;
  this->msg_size = _msg_size;
}

comm_pattern_key::comm_pattern_key(const comm_pattern_key& _copy){
  this->pattern_index = _copy.pattern_index;
  this->tag = _copy.tag;
  if (_copy.dim_sizes != nullptr) std::memcpy(&this->dim_sizes[0],&_copy.dim_sizes[0],2*sizeof(int));
  if (_copy.dim_strides != nullptr) std::memcpy(&this->dim_strides[0],&_copy.dim_strides[0],2*sizeof(int));
  this->msg_size = _copy.msg_size;
  this->partner_offset = _copy.partner_offset;
}

comm_pattern_key& comm_pattern_key::operator=(const comm_pattern_key& _copy){
  this->pattern_index = _copy.pattern_index;
  this->tag = _copy.tag;
  if (_copy.dim_sizes != nullptr) std::memcpy(&this->dim_sizes[0],&_copy.dim_sizes[0],2*sizeof(int));
  if (_copy.dim_strides != nullptr) std::memcpy(&this->dim_strides[0],&_copy.dim_strides[0],2*sizeof(int));
  this->msg_size = _copy.msg_size;
  this->partner_offset = _copy.partner_offset;
  return *this;
}

bool operator==(const comm_pattern_key& ref1, const comm_pattern_key& ref2){
  // Note that because of how we set the member variables in the constructor based on envlope, unit, and analysis parameterizations, no branching is required here.
  if ((ref1.tag==ref2.tag) &&
      (ref1.dim_sizes[0] == ref2.dim_sizes[0]) && (ref1.dim_sizes[1] == ref2.dim_sizes[1]) && (ref1.dim_sizes[2] == ref2.dim_sizes[2]) &&
      (ref1.dim_strides[0] == ref2.dim_strides[0]) && (ref1.dim_strides[1] == ref2.dim_strides[1]) && (ref1.dim_strides[2] == ref2.dim_strides[2]) &&
      (ref1.msg_size == ref2.msg_size) && (ref1.partner_offset == ref2.partner_offset)) return true;
  else return false;
}

bool operator<(const comm_pattern_key& ref1, const comm_pattern_key& ref2){
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
comp_pattern_key::comp_pattern_key(int _pattern_index, int _tag, double _flops, int _param1, int _param2, int _param3, int _param4, int _param5){
  this->pattern_index = _pattern_index;
  this->tag = _tag;
  this->param1 = truncate(_param1,comp_unit_param);
  this->param2 = truncate(_param2,comp_unit_param);
  this->param3 = truncate(_param3,comp_unit_param);
  this->param4 = truncate(_param4,comp_unit_param);
  this->param5 = truncate(_param5,comp_unit_param);
  this->flops = truncate(_flops,comp_unit_param);
}

comp_pattern_key::comp_pattern_key(const comp_pattern_key& _copy){
  this->pattern_index = _copy.pattern_index;
  this->tag = _copy.tag;
  this->flops = _copy.flops;
  this->param1 = _copy.param1;
  this->param2 = _copy.param2;
  this->param3 = _copy.param3;
  this->param4 = _copy.param4;
  this->param5 = _copy.param5;
}

comp_pattern_key& comp_pattern_key::operator=(const comp_pattern_key& _copy){
  this->pattern_index = _copy.pattern_index;
  this->tag = _copy.tag;
  this->flops = _copy.flops;
  this->param1 = _copy.param1;
  this->param2 = _copy.param2;
  this->param3 = _copy.param3;
  this->param4 = _copy.param4;
  this->param5 = _copy.param5;
  return *this;
}

bool operator==(const comp_pattern_key& ref1, const comp_pattern_key& ref2){
  if ((ref1.tag==ref2.tag) && (ref1.param1 == ref2.param1) && (ref1.param2 == ref2.param2) && (ref1.param3 == ref2.param3) && (ref1.param4 == ref2.param4) && (ref1.param5 == ref2.param5)) return true;
  else return false;
}

bool operator<(const comp_pattern_key& ref1, const comp_pattern_key& ref2){
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
pattern_key_id::pattern_key_id(bool _is_active, int _key_index, int _val_index, bool _is_updated){
  this->is_active = _is_active;
  this->key_index = _key_index;
  this->val_index = _val_index;
  this->is_updated = _is_updated;
}

pattern_key_id::pattern_key_id(const pattern_key_id& _copy){
  this->is_active = _copy.is_active;
  this->key_index = _copy.key_index;
  this->val_index = _copy.val_index;
  this->is_updated = _copy.is_updated;
}

pattern_key_id& pattern_key_id::operator=(const pattern_key_id& _copy){
  this->is_active = _copy.is_active;
  this->key_index = _copy.key_index;
  this->val_index = _copy.val_index;
  this->is_updated = _copy.is_updated;
  return *this;
}


int bsp_counter;
volatile double computation_timer;
std::vector<double> wall_timer;
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
}
}
