#ifndef CRITTER__PROFILE__PATH__PATH_H_
#define CRITTER__PROFILE__PATH__PATH_H_

#include "../container/comm_tracker.h"

namespace critter{
namespace internal{
namespace profile{

class path{
public:
  static void exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm);
  static bool initiate_comp(size_t id, volatile double curtime, float flop_count, int param1, int param2, int param3, int param4, int param5);
  static void complete_comp(double errtime, size_t id, float flop_count, int param1, int param2, int param3, int param4, int param5);
  static bool initiate_comm(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                       bool is_sender, int partner1, int user_tag1, int partner2, int user_tag2);
  static bool initiate_comm(nonblocking& tracker, volatile double curtime, int64_t nelem,
                       MPI_Datatype t, MPI_Comm comm, bool is_sender, int partner, int user_tag);
  static void initiate_comm(nonblocking& tracker, volatile double itime, int64_t nelem,
                       MPI_Datatype t, MPI_Comm comm, MPI_Request* request, bool is_sender, int partner, int user_tag);
  static void complete_comm(blocking& tracker);
  static void complete_comm(double comp_time);
  static int complete_comm(double curtime, MPI_Request* request, MPI_Status* status, int is_test, int* flag);
  static int complete_comm(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status, int is_test, int* flag);
  static int complete_comm(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[], int is_test);
  static int complete_comm(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[], int is_test, int* flag);
  static void propagate(blocking& tracker);
  static void propagate(nonblocking& tracker, float*& path_data, MPI_Request* prop_req);

private:
  static void complete_comm(nonblocking& tracker, MPI_Request* request, double comp_time, double comm_time);
};

}
}
}

#endif /*CRITTER__PROFILE__PATH__PATH_H_*/
