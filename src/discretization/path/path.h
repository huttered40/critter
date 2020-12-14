#ifndef CRITTER__DISCRETIZATION__PATH__PATH_H_
#define CRITTER__DISCRETIZATION__PATH__PATH_H_

#include "../container/comm_tracker.h"

namespace critter{
namespace internal{
namespace discretization{

class path{
public:
  static void exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm);
  static bool initiate_comp(size_t id, volatile float curtime, float flop_count, int param1, int param2, int param3, int param4, int param5);
  static void complete_comp(float errtime, size_t id, float flop_count, int param1, int param2, int param3, int param4, int param5);
  static bool initiate_comm(blocking& tracker, volatile float curtime, int64_t nelem, MPI_Datatype t, MPI_Comm comm,
                       bool is_sender, int partner1, int partner2);
  static bool initiate_comm(nonblocking& tracker, volatile float curtime, int64_t nelem,
                       MPI_Datatype t, MPI_Comm comm, int user_tag, bool is_sender, int partner);
  static void initiate_comm(nonblocking& tracker, volatile float itime, int64_t nelem,
                       MPI_Datatype t, MPI_Comm comm, MPI_Request* request, int user_tag, bool is_sender, int partner);
  static void complete_comm(blocking& tracker, int recv_source=-1);
  static int complete_comm(float curtime, MPI_Request* request, MPI_Status* status);
  static int complete_comm(float curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status);
  static int complete_comm(float curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]);
  static int complete_comm(float curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);

  static void comp_state_aggregation(blocking& tracker);
  static void comm_state_aggregation(blocking& tracker);

private:
  static void complete_comm(nonblocking& tracker, MPI_Request* request, float comp_time, float comm_time);
  static void single_stage_sample_aggregation(blocking& tracker);
  static void multi_stage_sample_aggregation(blocking& tracker);
};

}
}
}

#endif /*CRITTER__DISCRETIZATION__PATH__PATH_H_*/
