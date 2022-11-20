#ifndef CRITTER__PROFILE__PROFILER_H_
#define CRITTER__PROFILE__PROFILER_H_

#include "container/comm_tracker.h"

namespace internal{

class profiler{
public:
  static bool initiate_comp(volatile double curtime, float flop_count);
  static void complete_comp();
  static bool initiate_comm(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                       bool is_sender=false, int partner1=-1, int user_tag1=-1, int partner2=-1, int user_tag2=-1);
  static bool inspect_comm(nonblocking& tracker, volatile double curtime, int64_t nelem,
                       MPI_Datatype t, MPI_Comm comm, bool is_sender=false, int partner=-1, int user_tag=-1);
  static void initiate_comm(nonblocking& tracker, volatile double itime, int64_t nelem,
                       MPI_Datatype t, MPI_Comm comm, MPI_Request* request, bool is_sender=false, int partner=-1, int user_tag=-1);
  static void complete_comm(blocking& tracker);
  static void complete_comm(double comp_time);
  static int complete_comm(double curtime, MPI_Request* request, MPI_Status* status, int is_test=0, int* flag=nullptr);
  static int complete_comm(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status, int is_test=0, int* flag=nullptr);
  static int complete_comm(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[], int is_test=0);
  static int complete_comm(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[], int is_test=0, int* flag=nullptr);
  static void propagate(blocking& tracker);
  static void propagate(nonblocking& tracker, float*& path_data, MPI_Request* prop_req);

private:
  static void complete_comm(nonblocking& tracker, MPI_Request* request, double comp_time, double comm_time);
};

}

#endif /*CRITTER__PROFILE__PROFILER_H_*/
