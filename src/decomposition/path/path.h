#ifndef CRITTER__DECOMPOSITION__PATH__PATH_H_
#define CRITTER__DECOMPOSITION__PATH__PATH_H_

#include "../container/comm_tracker.h"

namespace critter{
namespace internal{
namespace decomposition{

class path{
public:
  static void exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm);
  static bool initiate_comp(size_t id, volatile double curtime, double flop_count, int param1, int param2, int param3, int param4, int param5);
  static void complete_comp(size_t id, double flop_count, int param1, int param2, int param3, int param4, int param5);
  static bool initiate_comm(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                       bool is_sender, int partner1, int partner2);
  static bool initiate_comm(nonblocking& tracker, volatile double curtime, volatile double itime, int64_t nelem,
                       MPI_Datatype t, MPI_Comm comm, MPI_Request* request, bool is_sender, int partner);
  static void complete_comm(blocking& tracker, int recv_source=-1);
  static void complete_comm(double curtime, MPI_Request* request, MPI_Status* status);
  static void complete_comm(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status);
  static void complete_comm(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]);
  static void complete_comm(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);
  static void propagate(blocking& tracker);
  static void propagate(nonblocking& tracker);

private:
  static void complete_comm(nonblocking& tracker, MPI_Request* request, double comp_time, double comm_time);
  static void propagate_symbols(blocking& tracker, int rank);
  static void propagate_symbols(nonblocking& tracker, int rank);
};

}
}
}

#endif /*CRITTER__DECOMPOSITION__PATH__PATH_H_*/
