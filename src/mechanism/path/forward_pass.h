#ifndef CRITTER__MECHANISM__PATH__FORWARD_PASS_H_
#define CRITTER__MECHANISM__PATH__FORWARD_PASS_H_

#include "../../container/comm_tracker.h"

namespace critter{
namespace internal{

class forward_pass{
public:
  static void initiate(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
                       bool is_sender, int partner1, int partner2);
  static void initiate(nonblocking& tracker, volatile double curtime, volatile double itime, int64_t nelem,
                       MPI_Datatype t, MPI_Comm cm, MPI_Request* request, bool is_sender, int partner);
  static void complete(blocking& tracker);
  static void wait(double curtime, MPI_Request* request, MPI_Status* status);
  static void wait(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status);
  static void wait(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]);
  static void wait(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);
  static void propagate(MPI_Comm cm, int tag, bool is_sender, int partner1, int partner2);

private:
  static void wait(int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status);
  static void complete(nonblocking& tracker, MPI_Request* request, double comp_time, double comm_time);
  static void initiate_timers(int rank, int tag, MPI_Comm cm, int partner1, int partner2);
  static void propagate_timers(int rank, int tag, MPI_Comm cm, int partner1, int partner2);
};

}
}

#endif /*CRITTER__MECHANISM__PATH__FORWARD_PASS_H_*/
