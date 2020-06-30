#include "dispatch.h"

namespace critter{
namespace internal{

void initiate(blocking& tracker, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
              bool is_sender, int partner1, int partner2){
  switch (mechanism){
    case 0:
      forward_pass::initiate(tracker,curtime,nelem,t,cm,is_sender,partner1,partner2);
      break;
  }
}

void initiate(nonblocking& tracker, volatile double curtime, volatile double itime, int64_t nelem,
              MPI_Datatype t, MPI_Comm cm, MPI_Request* request, bool is_sender, int partner){
  switch (mechanism){
    case 0:
      forward_pass::initiate(tracker,curtime,itime,nelem,t,cm,request,is_sender,partner);
      break;
  }
}

void complete(blocking& tracker){
  switch (mechanism){
    case 0:
      forward_pass::complete(tracker);
      break;
  }
}

void wait(double curtime, MPI_Request* request, MPI_Status* status){
  switch (mechanism){
    case 0:
      forward_pass::wait(curtime,request,status);
      break;
  }
}

void wait(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){
  switch (mechanism){
    case 0:
      forward_pass::wait(curtime,count,array_of_requests,indx,status);
      break;
  }
}

void wait(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]){
  switch (mechanism){
    case 0:
      forward_pass::wait(curtime,incount,array_of_requests,outcount,array_of_indices,array_of_statuses);
      break;
  }
}

void wait(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
  switch (mechanism){
    case 0:
      forward_pass::wait(curtime,count,array_of_requests,array_of_statuses);
      break;
  }
}

}
}
