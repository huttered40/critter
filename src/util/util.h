#ifndef CRITTER__UTIL__UTIL_H_
#define CRITTER__UTIL__UTIL_H_

#include <mpi.h>
#include <cstring>
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

#include "../discretization/parameterization/parameterization.h"

namespace critter{
namespace internal{

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
struct event{
  event(std::string _kernel, std::vector<double> _measurements){
    tag = -1; measurements = _measurements;
    kernel = _kernel;
  }
  event(std::string _kernel, std::vector<double> _measurements, int _tag, MPI_Comm _comm, int _partner1, int _partner2, bool _is_sender, bool _is_eager){
    measurements = _measurements;
    tag = _tag; comm = _comm; partner1 = _partner1; partner2 = _partner2; is_sender = _is_sender; is_eager = _is_eager;
    kernel = _kernel;
  }
  event(std::string _kernel, std::vector<double> _measurements, int _tag, MPI_Comm _comm, int _partner1, bool _is_sender, bool _is_eager, int _match_id, bool blah){
    measurements = _measurements;
    tag = _tag; comm = _comm; partner1 = _partner1; is_sender = _is_sender; is_eager = _is_eager; match_id = _match_id; is_close = false;;
    kernel = _kernel;
  }
  event(std::string _kernel,std::vector<double> _measurements, std::vector<int> _match_vec){
    measurements = _measurements;
    tag=18; match_vec = _match_vec; match_size = _match_vec.size(); is_close = true;
    kernel = _kernel;
  }
  MPI_Comm comm;
  int tag,partner1,partner2,match_size,match_id;
  bool is_sender,is_eager,is_close;
  std::vector<int> match_vec;
  std::vector<double> measurements;
  std::string kernel;
};

// ****************************************************************************************************************************************************
extern volatile double computation_timer;
extern std::vector<double> wall_timer;
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
//extern std::map<std::pair<std::string,size_t>,bool> schedule_map;

}
}

#endif /*CRITTER__UTIL__UTIL_H_*/
