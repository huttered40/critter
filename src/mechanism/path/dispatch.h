#ifndef CRITTER__MECHANISMS__PATH__DISPATCH_H_
#define CRITTER__MECHANISMS__PATH__DISPATCH_H_

#include "../../container/comm_tracker.h"
#include "decomposition.h"

namespace critter{
namespace internal{

void allocate(MPI_Comm comm);
void initiate(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
              bool is_sender=false, int partner1=-1, int partner2=-1);
void initiate(nonblocking& tracker, volatile double curtime, volatile double itime, int64_t nelem,
              MPI_Datatype t, MPI_Comm cm, MPI_Request* request, bool is_sender=false, int partner=-1);
void complete(blocking& tracker, int recv_source=-1);
void complete(double curtime, MPI_Request* request, MPI_Status* status);
void complete(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status);
void complete(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]);
void complete(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);
void propagate(MPI_Comm comm);
void final_accumulate(double last_time);

}
}

#endif /*CRITTER__MECHANISMS__PATH__DISPATCH_H_*/
