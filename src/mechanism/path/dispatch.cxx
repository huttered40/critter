#include "dispatch.h"

namespace critter{
namespace internal{

void allocate(MPI_Comm comm){
  switch (mechanism){
    case 0:
      decomposition::allocate(comm);
  }
}

void initiate(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
              bool is_sender, int partner1, int partner2){
  switch (mechanism){
    case 0:
      decomposition::initiate(tracker,curtime,nelem,t,cm,is_sender,partner1,partner2);
      break;
  }
}

void initiate(nonblocking& tracker, volatile double curtime, volatile double itime, int64_t nelem,
              MPI_Datatype t, MPI_Comm cm, MPI_Request* request, bool is_sender, int partner){
  switch (mechanism){
    case 0:
      decomposition::initiate(tracker,curtime,itime,nelem,t,cm,request,is_sender,partner);
      break;
  }
}

void complete(blocking& tracker, int recv_source){
  switch (mechanism){
    case 0:
      decomposition::complete(tracker,recv_source);
      break;
  }
}

void complete(double curtime, MPI_Request* request, MPI_Status* status){
  switch (mechanism){
    case 0:
      decomposition::complete(curtime,request,status);
      break;
  }
}

void complete(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){
  switch (mechanism){
    case 0:
      decomposition::complete(curtime,count,array_of_requests,indx,status);
      break;
  }
}

void complete(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]){
  switch (mechanism){
    case 0:
      decomposition::complete(curtime,incount,array_of_requests,outcount,array_of_indices,array_of_statuses);
      break;
  }
}

void complete(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
  switch (mechanism){
    case 0:
      decomposition::complete(curtime,count,array_of_requests,array_of_statuses);
      break;
  }
}

void propagate(MPI_Comm comm){
  _MPI_Barrier.comm = comm;
  switch (mechanism){
    case 0:
      decomposition::propagate(_MPI_Barrier);
      break;
  }
}

void final_accumulate(double last_time){
  switch (mechanism){
    case 0:
      decomposition::final_accumulate(last_time);
  }
}

}
}
