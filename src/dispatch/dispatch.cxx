#include "dispatch.h"
#include "../decomposition/container/comm_tracker.h"
#include "../decomposition/util/util.h"
#include "../decomposition/volumetric/volumetric.h"
#include "../decomposition/path/path.h"
#include "../decomposition/record/record.h"
#include "../discretization/container/comm_tracker.h"
#include "../discretization/util/util.h"
#include "../discretization/volumetric/volumetric.h"
#include "../discretization/path/path.h"
#include "../discretization/record/record.h"
#include "../replay/path/path.h"

namespace critter{
namespace internal{

void allocate(MPI_Comm comm){
  switch (mechanism){
    case 0:
      decomposition::allocate(comm);
      break;
    case 1:
      discretization::allocate(comm);
      break;
  }
}

void reset(bool schedule_kernels_override, bool force_steady_statistical_data_overide){
  switch (mechanism){
    case 0:
      decomposition::reset();
      break;
    case 1:
      discretization::reset(schedule_kernels_override,force_steady_statistical_data_overide);
      break;
  }
}

void exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm){
  switch (mechanism){
    case 0:
      decomposition::path::exchange_communicators(oldcomm,newcomm);
      break;
    case 1:
      discretization::path::exchange_communicators(oldcomm,newcomm);
      break;
  }
}

bool initiate_comp(size_t id, volatile double curtime, double flop_count,  int param1, int param2, int param3, int param4, int param5){
  bool schedule_decision;
  switch (mechanism){
    case 0:
      schedule_decision = decomposition::path::initiate_comp(id,curtime,flop_count,param1,param2,param3,param4,param5);
      break;
    case 1:
      schedule_decision = discretization::path::initiate_comp(id,curtime,flop_count,param1,param2,param3,param4,param5);
      break;
  }
  return schedule_decision;
}

void complete_comp(size_t id, double flop_count,  int param1, int param2, int param3, int param4, int param5){
  switch (mechanism){
    case 0:
      decomposition::path::complete_comp(id,flop_count,param1,param2,param3,param4,param5);
      break;
    case 1:
      discretization::path::complete_comp(id,flop_count,param1,param2,param3,param4,param5);
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
    case 1:
      schedule_decision = discretization::path::initiate_comm(*(discretization::blocking*)discretization::list[id],curtime,nelem,t,cm,is_sender,partner1,partner2);
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
    case 1:
      schedule_decision = discretization::path::initiate_comm(*(discretization::nonblocking*)discretization::list[id],curtime,itime,nelem,t,cm,request,is_sender,partner);
      break;
  }
  return schedule_decision;
}

void complete_comm(size_t id, int recv_source){
  switch (mechanism){
    case 0:
      decomposition::path::complete_comm(*(decomposition::blocking*)decomposition::list[id],recv_source);
      break;
    case 1:
      discretization::path::complete_comm(*(discretization::blocking*)discretization::list[id],recv_source);
      break;
  }
}

void complete_comm(double curtime, MPI_Request* request, MPI_Status* status){
  switch (mechanism){
    case 0:
      decomposition::path::complete_comm(curtime,request,status);
      break;
    case 1:
      discretization::path::complete_comm(curtime,request,status);
      break;
  }
}

void complete_comm(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){
  switch (mechanism){
    case 0:
      decomposition::path::complete_comm(curtime,count,array_of_requests,indx,status);
      break;
    case 1:
      discretization::path::complete_comm(curtime,count,array_of_requests,indx,status);
      break;
  }
}

void complete_comm(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]){
  switch (mechanism){
    case 0:
      decomposition::path::complete_comm(curtime,incount,array_of_requests,outcount,array_of_indices,array_of_statuses);
      break;
    case 1:
      discretization::path::complete_comm(curtime,incount,array_of_requests,outcount,array_of_indices,array_of_statuses);
      break;
  }
}

void complete_comm(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
  switch (mechanism){
    case 0:
      decomposition::path::complete_comm(curtime,count,array_of_requests,array_of_statuses);
      break;
    case 1:
      discretization::path::complete_comm(curtime,count,array_of_requests,array_of_statuses);
      break;
  }
}

void propagate(MPI_Comm comm){
  switch (mechanism){
    case 0:
      decomposition::_MPI_Barrier.comm = comm;
      decomposition::path::propagate(decomposition::_MPI_Barrier);
      break;
    case 1:
      // Do nothing: 4-double reduction is performed in 'discretization::final_accumulate'
      break;
  }
}

void collect(MPI_Comm comm){
  switch (mechanism){
    case 0:
      decomposition::volumetric::collect(comm);
      break;
    case 1:
      discretization::volumetric::collect(comm);
      break;
  }
}

void final_accumulate(MPI_Comm comm, double last_time){
  switch (mechanism){
    case 0:
      decomposition::final_accumulate(comm,last_time);
      break;
    case 1:
      discretization::final_accumulate(comm,last_time);
      break;
  }
}

void open_symbol(const char* symbol, double curtime){
  switch (mechanism){
    case 0:
      decomposition::open_symbol(symbol,curtime);
      break;
    case 1:
      discretization::open_symbol(symbol,curtime);
      break;
  }
}

void close_symbol(const char* symbol, double curtime){
  switch (mechanism){
    case 0:
      decomposition::close_symbol(symbol,curtime);
      break;
    case 1:
      discretization::close_symbol(symbol,curtime);
      break;
  }
}

void clear(){
  switch (mechanism){
    case 0:
      decomposition::clear();
      break;
    case 1:
      discretization::clear();
      break;
  }
}

void _finalize(){
  switch (mechanism){
    case 0:
      decomposition::finalize();
      break;
    case 1:
      discretization::finalize();
      break;
  }
}

void reset_frequencies(){
  switch (mechanism){
    case 0:
      break;
    case 1:
      discretization::reset_frequencies();
      break;
  }
}

void record(std::ofstream& Stream, double* data, bool print_statistical_data, bool save_statistical_data){
  switch (mechanism){
    case 0:
      decomposition::record::invoke(Stream,data);
      break;
    case 1:
      discretization::record::invoke(Stream,data,print_statistical_data,save_statistical_data);
      break;
  }
}

void record(std::ostream& Stream, double* data, bool print_statistical_data, bool save_statistical_data){
  switch (mechanism){
    case 0:
      decomposition::record::invoke(Stream);
      break;
    case 1:
      discretization::record::invoke(Stream,data,print_statistical_data,save_statistical_data);
      break;
  }
}

}
}
