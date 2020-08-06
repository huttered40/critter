#include "dispatch.h"
#include "../decomposition/container/comm_tracker.h"
#include "../decomposition/util/util.h"
#include "../decomposition/volumetric/volumetric.h"
#include "../decomposition/path/path.h"
#include "../decomposition/record/record.h"
#include "../replay/path/path.h"

namespace critter{
namespace internal{

void allocate(MPI_Comm comm){
  switch (mechanism){
    case 0:
      decomposition::allocate(comm);
  }
}

void reset(){
  switch (mechanism){
    case 0:
      decomposition::reset();
  }
}

void exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm){
  switch (mechanism){
    case 0:
      decomposition::path::exchange_communicators(oldcomm,newcomm);
      break;
  }
}

bool initiate_comp(size_t id, volatile double curtime, double flop_count,  int param1, int param2, int param3, int param4, int param5){
  bool schedule_decision;
  switch (mechanism){
    case 0:
      schedule_decision = decomposition::path::initiate_comp(id,curtime,flop_count,param1,param2,param3,param4,param5);
      break;
  }
  return schedule_decision;
}

void complete_comp(size_t id, double flop_count,  int param1, int param2, int param3, int param4, int param5){
  switch (mechanism){
    case 0:
      decomposition::path::complete_comp(id,flop_count,param1,param2,param3,param4,param5);
      break;
  }
}

bool initiate_comm(size_t id, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
              bool is_sender, int partner1, int partner2){
  bool schedule_decision;
  switch (mechanism){
    case 0:
      schedule_decision = decomposition::path::initiate_comm(*(decomposition::blocking*)decomposition::list[id],curtime,nelem,t,cm,is_sender,partner1,partner2);
      break;
  }
  return schedule_decision;
}

bool initiate_comm(size_t id, volatile double curtime, volatile double itime, int64_t nelem,
              MPI_Datatype t, MPI_Comm cm, MPI_Request* request, bool is_sender, int partner){
  bool schedule_decision;
  switch (mechanism){
    case 0:
      schedule_decision = decomposition::path::initiate_comm(*(decomposition::nonblocking*)decomposition::list[id],curtime,itime,nelem,t,cm,request,is_sender,partner);
      break;
  }
  return schedule_decision;
}

void complete_comm(size_t id, int recv_source){
  switch (mechanism){
    case 0:
      decomposition::path::complete_comm(*(decomposition::blocking*)decomposition::list[id],recv_source);
      break;
  }
}

void complete_comm(double curtime, MPI_Request* request, MPI_Status* status){
  switch (mechanism){
    case 0:
      decomposition::path::complete_comm(curtime,request,status);
      break;
  }
}

void complete_comm(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){
  switch (mechanism){
    case 0:
      decomposition::path::complete_comm(curtime,count,array_of_requests,indx,status);
      break;
  }
}

void complete_comm(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]){
  switch (mechanism){
    case 0:
      decomposition::path::complete_comm(curtime,incount,array_of_requests,outcount,array_of_indices,array_of_statuses);
      break;
  }
}

void complete_comm(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
  switch (mechanism){
    case 0:
      decomposition::path::complete_comm(curtime,count,array_of_requests,array_of_statuses);
      break;
  }
}

void propagate(MPI_Comm comm){
  switch (mechanism){
    case 0:
      decomposition::_MPI_Barrier.comm = comm;
      decomposition::path::propagate(decomposition::_MPI_Barrier);
      break;
  }
}

void collect(MPI_Comm comm){
  switch (mechanism){
    case 0:
      decomposition::volumetric::collect(comm);
      break;
  }
}

void final_accumulate(double last_time){
  switch (mechanism){
    case 0:
      decomposition::final_accumulate(last_time);
      break;
  }
}

void open_symbol(const char* symbol, double curtime){
  switch (mechanism){
    case 0:
      decomposition::open_symbol(symbol,curtime);
      break;
    case 1:
      //discretization::open_symbol(symbol,curtime);
      break;
  }
}

void close_symbol(const char* symbol, double curtime){
  switch (mechanism){
    case 0:
      decomposition::close_symbol(symbol,curtime);
      break;
    case 1:
      //discretization::close_symbol(symbol,curtime);
      break;
  }
}

void clear(){
  switch (mechanism){
    case 0:
      decomposition::clear();
      break;
  }
}

void record(std::ofstream& Stream){
  switch (mechanism){
    case 0:
      decomposition::record::invoke(Stream);
      break;
  }
}

void record(std::ostream& Stream){
  switch (mechanism){
    case 0:
      decomposition::record::invoke(Stream);
      break;
  }
}

}
}
