#ifndef CRITTER__DISPATCH__DISPATCH_H_
#define CRITTER__DISPATCH__DISPATCH_H_

#include "../util/util.h"

namespace critter{
namespace internal{

void allocate(MPI_Comm comm);
void reset();

void initiate(size_t id, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
              bool is_sender=false, int partner1=-1, int partner2=-1);
void initiate(size_t id, volatile double curtime, volatile double itime, int64_t nelem,
              MPI_Datatype t, MPI_Comm cm, MPI_Request* request, bool is_sender=false, int partner=-1);
void complete(size_t id, int recv_source=-1);
void complete(double curtime, MPI_Request* request, MPI_Status* status);
void complete(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status);
void complete(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]);
void complete(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);
void propagate(MPI_Comm comm);
void collect(MPI_Comm comm);
void final_accumulate(double last_time);

void open_symbol(const char* symbol, double curtime);
void close_symbol(const char* symbol, double curtime);

void record(std::ofstream& Stream);
void record(std::ostream& Stream);

void clear();

}
}

#endif /*CRITTER__DISPATCH__DISPATCH_H_*/
