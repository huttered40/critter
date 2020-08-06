#ifndef CRITTER__DISPATCH__DISPATCH_H_
#define CRITTER__DISPATCH__DISPATCH_H_

#include "../util/util.h"

namespace critter{
namespace internal{

void allocate(MPI_Comm comm);
void reset();

void exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm);
bool initiate_comp(size_t id, volatile double curtime, double flop_count, int param1=-1, int param2=-1, int param3=-1, int param4=-1, int param5=-1);
void complete_comp(size_t id, double flop_count, int param1=-1, int param2=-1, int param3=-1, int param4=-1, int param5=-1);
bool initiate_comm(size_t id, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
              bool is_sender=false, int partner1=-1, int partner2=-1);
bool initiate_comm(size_t id, volatile double curtime, volatile double itime, int64_t nelem,
              MPI_Datatype t, MPI_Comm cm, MPI_Request* request, bool is_sender=false, int partner=-1);
void complete_comm(size_t id, int recv_source=-1);
void complete_comm(double curtime, MPI_Request* request, MPI_Status* status);
void complete_comm(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status);
void complete_comm(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]);
void complete_comm(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);
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
