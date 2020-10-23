#ifndef CRITTER__UTIL__UTIL_H_
#define CRITTER__UTIL__UTIL_H_

#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <vector>
#include <stack>
#include <stdint.h>
#include <functional>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_map>
#include <cmath>
#include <assert.h>

namespace critter{
namespace internal{

// ****************************************************************************************************************************************************
extern int comm_envelope_param;
extern int comp_envelope_param;
extern int comm_unit_param;
extern int comp_unit_param;
extern int comm_analysis_param;
extern int comp_analysis_param;

// ****************************************************************************************************************************************************
struct kernel_key{};

// ****************************************************************************************************************************************************
struct comm_kernel_key : public kernel_key{

  comm_kernel_key(int _rank=-1, int _kernel_index=-1, int _tag=-1, int _dim_sizes[2]=nullptr, int _dim_strides[2]=nullptr, double _msg_size=-1, int _partner=-1);
  comm_kernel_key(int _kernel_index, int _tag, int _dim_sizes[2], int _dim_strides[2], double _msg_size, int _partner_offset);
  comm_kernel_key(const comm_kernel_key& _copy);
  comm_kernel_key& operator=(const comm_kernel_key& _copy);
  friend bool operator==(const comm_kernel_key& ref1, const comm_kernel_key& ref2);
  friend bool operator<(const comm_kernel_key& ref1, const comm_kernel_key& ref2);

  int tag;
  int dim_sizes[2];// Allow up to 2 dimensions
  int dim_strides[2];// Allow up to 2 dimensions
  int partner_offset;
  int kernel_index;
  double msg_size;
};


// ****************************************************************************************************************************************************
struct comp_kernel_key : public kernel_key{

  comp_kernel_key(int _kernel_index=-1, int _tag=-1, double _flops=-1, int =-1, int _param2=-1, int _param3=-1, int _param4=-1, int _param5=-1);
  comp_kernel_key(const comp_kernel_key& _copy);
  comp_kernel_key& operator=(const comp_kernel_key& _copy);
  friend bool operator==(const comp_kernel_key& ref1, const comp_kernel_key& ref2);
  friend bool operator<(const comp_kernel_key& ref1, const comp_kernel_key& ref2);

  int tag;
  int param1,param2,param3,param4,param5;
  int kernel_index;
  double flops;
};


// ****************************************************************************************************************************************************
struct double_int{
  double_int(){first=0; second=0;}
  double_int(double one, int two){first=one; second=two;}
  double first; int second;
};
struct int_int_double{
  int_int_double(){first=0; second=0; third=0;}
  int_int_double(int one, int two, double three){first=one; second=two; third=three;}
  int first; int second; double third;
};

// ****************************************************************************************************************************************************
struct kernel_key_id{

  kernel_key_id(bool _is_active=false, int _key_index=0, int _val_index=0, bool _is_updated=false);
  kernel_key_id(const kernel_key_id& _copy);
  kernel_key_id& operator=(const kernel_key_id& _copy);

  // Active just means its still being propogated. It acts as a switch betweeh steady_state arrays and active arrays
  bool is_active;
  bool is_updated;
  int key_index;
  int val_index;
};

// ****************************************************************************************************************************************************
extern int bsp_counter;
extern int reset_counter;
extern int clear_counter;
extern int communicator_count;
extern std::string schedule_tag;
extern volatile double computation_timer;
extern std::vector<double> wall_timer;
extern double _wall_time;
extern size_t auto_capture;
extern bool is_world_root;
extern size_t mechanism,mode,stack_id;
extern double scratch_pad;
extern size_t track_blas;
extern size_t track_lapack;
extern size_t track_collective;
extern size_t track_p2p;
extern size_t track_p2p_idle;
extern size_t eager_p2p;
extern size_t delete_comm;
extern size_t
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
extern size_t
	_BLAS_axpy__id,
	_BLAS_scal__id,
	_BLAS_ger__id,
	_BLAS_gemm__id,
	_BLAS_trmm__id,
	_BLAS_trsm__id,
	_BLAS_syrk__id;
extern size_t
	_LAPACK_getrf__id,
	_LAPACK_potrf__id,
	_LAPACK_trtri__id,
	_LAPACK_geqrf__id,
	_LAPACK_orgqr__id,
	_LAPACK_ormqr__id,
	_LAPACK_getri__id,
	_LAPACK_tpqrt__id,
	_LAPACK_tpmqrt__id;
extern size_t
	_CAPITAL_blktocyc__id;

// ****************************************************************************************************************************************************
struct channel{
  channel();
  static std::vector<std::pair<int,int>> generate_tuple(std::vector<int>& ranks, int new_comm_size);
  static void contract_tuple(std::vector<std::pair<int,int>>& tuple_list);
  static int enumerate_tuple(channel* node, std::vector<int>& process_list);
  static int duplicate_process_count(std::vector<int>& process_list);
  static int translate_rank(MPI_Comm comm, int rank);
  static bool verify_ancestor_relation(channel* comm1, channel* comm2);
  static bool verify_sibling_relation(channel* comm1, channel* comm2);
  static int span(std::pair<int,int>& id);
  static std::string generate_tuple_string(channel* comm);
 
  int offset;
  int local_hash_tag;
  int global_hash_tag;
  std::vector<std::pair<int,int>> id;
};

struct aggregate_channel : public channel{
  aggregate_channel(std::vector<std::pair<int,int>>& tuple_list, int local_hash, int global_hash, int offset, int channel_size);
  static std::string generate_hash_history(aggregate_channel* comm);

  bool is_final;
  int num_channels;
  std::set<int> channels;
};

struct solo_channel : public channel{
  solo_channel();
  static bool verify_sibling_relation(solo_channel* node, int subtree_idx, std::vector<int>& skip_indices);

  int tag;
  int frequency;
  solo_channel* parent;
  std::vector<std::vector<solo_channel*>> children;
};

extern std::map<MPI_Comm,solo_channel*> comm_channel_map;
extern std::map<int,solo_channel*> p2p_channel_map;
extern std::map<int,aggregate_channel*> aggregate_channel_map;

void generate_aggregate_channels(MPI_Comm oldcomm, MPI_Comm newcomm);

}
}

#endif /*CRITTER__UTIL__UTIL_H_*/
