#ifndef CRITTER__SKELETONIZATION__PATH__PATH_H_
#define CRITTER__SKELETONIZATION__PATH__PATH_H_

#include "../container/comm_tracker.h"

namespace critter{
namespace internal{
namespace skeletonization{

class path{
public:
  static void exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm);
  static bool initiate_comp(size_t id, double flop_count, int param1, int param2, int param3, int param4, int param5);
  static void complete_comp(size_t id, double flop_count, int param1, int param2, int param3, int param4, int param5);
  static bool initiate_comm(blocking& tracker, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                            bool is_sender, int partner1, int partner2);
  static bool initiate_comm(nonblocking& tracker, int64_t nelem, MPI_Datatype t, MPI_Comm comm, bool is_sender, int partner);
  static void initiate_comm(nonblocking& tracker, int64_t nelem,
                       MPI_Datatype t, MPI_Comm comm, MPI_Request* request, bool is_sender, int partner);
  static void complete_comm(blocking& tracker, int recv_source=-1);
  static int complete_comm(MPI_Request* request, MPI_Status* status);
  static int complete_comm(int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status);
  static int complete_comm(int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]);
  static int complete_comm(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);

private:
  static void complete_comm(nonblocking& tracker, MPI_Request* request);
  static void propagate_kernels(blocking& tracker);
};

}
}
}

#endif /*CRITTER__SKELETONIZATION__PATH__PATH_H_*/
