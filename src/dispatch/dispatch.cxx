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
#include "../skeleton/container/comm_tracker.h"
#include "../skeleton/util/util.h"
#include "../skeleton/volumetric/volumetric.h"
#include "../skeleton/path/path.h"
#include "../skeleton/record/record.h"

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
    case 2:
      skeleton::allocate(comm);
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
    case 2:
      skeleton::reset();
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
    case 2:
      skeleton::path::exchange_communicators(oldcomm,newcomm);
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
    case 2:
      schedule_decision = skeleton::path::initiate_comp(id,curtime,flop_count,param1,param2,param3,param4,param5);
      break;
  }
  return schedule_decision;
}

void complete_comp(double errtime, size_t id, double flop_count,  int param1, int param2, int param3, int param4, int param5){
  switch (mechanism){
    case 0:
      decomposition::path::complete_comp(errtime,id,flop_count,param1,param2,param3,param4,param5);
      break;
    case 1:
      discretization::path::complete_comp(errtime,id,flop_count,param1,param2,param3,param4,param5);
      break;
    case 2:
      skeleton::path::complete_comp(errtime,id,flop_count,param1,param2,param3,param4,param5);
      break;
  }
}

bool initiate_comm(size_t id, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
              bool is_sender, int partner1, int user_tag1, int partner2, int user_tag2){
  bool schedule_decision;
  switch (mechanism){
    case 0:
      schedule_decision = decomposition::path::initiate_comm(*(decomposition::blocking*)decomposition::list[id],curtime,nelem,t,cm,is_sender,partner1,user_tag1,partner2,user_tag2);
      break;
    case 1:
      schedule_decision = discretization::path::initiate_comm(*(discretization::blocking*)discretization::list[id],curtime,nelem,t,cm,is_sender,partner1,user_tag1,partner2,user_tag2);
      break;
    case 2:
      schedule_decision = skeleton::path::initiate_comm(*(skeleton::blocking*)skeleton::list[id],curtime,nelem,t,cm,is_sender,partner1,user_tag1,partner2,user_tag2);
      break;
  }
  return schedule_decision;
}

bool inspect_comm(size_t id, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
              bool is_sender, int partner, int user_tag){
  bool schedule_decision;
  switch (mechanism){
    case 0:
      schedule_decision = decomposition::path::initiate_comm(*(decomposition::nonblocking*)decomposition::list[id],curtime,nelem,t,cm,is_sender,partner,user_tag);
      break;
    case 1:
      schedule_decision = discretization::path::initiate_comm(*(discretization::nonblocking*)discretization::list[id],curtime,nelem,t,cm,is_sender,partner,user_tag);
      break;
    case 2:
      schedule_decision = skeleton::path::initiate_comm(*(skeleton::nonblocking*)skeleton::list[id],curtime,nelem,t,cm,is_sender,partner,user_tag);
      break;
  }
  return schedule_decision;
}

void initiate_comm(size_t id, volatile double itime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, MPI_Request* request, bool is_sender, int partner, int user_tag){
  switch (mechanism){
    case 0:
      decomposition::path::initiate_comm(*(decomposition::nonblocking*)decomposition::list[id],itime,nelem,t,cm,request,is_sender,partner,user_tag);
      break;
    case 1:
      discretization::path::initiate_comm(*(discretization::nonblocking*)discretization::list[id],itime,nelem,t,cm,request,is_sender,partner,user_tag);
      break;
    case 2:
      skeleton::path::initiate_comm(*(skeleton::nonblocking*)skeleton::list[id],itime,nelem,t,cm,request,is_sender,partner,user_tag);
      break;
  }
}

void complete_comm(size_t id){
  switch (mechanism){
    case 0:
      decomposition::path::complete_comm(*(decomposition::blocking*)decomposition::list[id]);
      break;
    case 1:
      discretization::path::complete_comm(*(discretization::blocking*)discretization::list[id]);
      break;
    case 2:
      skeleton::path::complete_comm(*(skeleton::blocking*)skeleton::list[id]);
      break;
  }
}

int complete_comm(double curtime, MPI_Request* request, MPI_Status* status, int is_test, int* flag){
  switch (mechanism){
    case 0:
      return decomposition::path::complete_comm(curtime,request,status,is_test,flag);
    case 1:
      return discretization::path::complete_comm(curtime,request,status);
    case 2:
      return skeleton::path::complete_comm(curtime,request,status);
  }
}

int complete_comm(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status, int is_test, int* flag){
  switch (mechanism){
    case 0:
      return decomposition::path::complete_comm(curtime,count,array_of_requests,indx,status,is_test,flag);
    case 1:
      return discretization::path::complete_comm(curtime,count,array_of_requests,indx,status);
    case 2:
      return skeleton::path::complete_comm(curtime,count,array_of_requests,indx,status);
  }
}

int complete_comm(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[], int is_test){
  switch (mechanism){
    case 0:
      return decomposition::path::complete_comm(curtime,incount,array_of_requests,outcount,array_of_indices,array_of_statuses,is_test);
    case 1:
      return discretization::path::complete_comm(curtime,incount,array_of_requests,outcount,array_of_indices,array_of_statuses);
    case 2:
      return skeleton::path::complete_comm(curtime,incount,array_of_requests,outcount,array_of_indices,array_of_statuses);
  }
}

int complete_comm(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[], int is_test, int* flag){
  switch (mechanism){
    case 0:
      return decomposition::path::complete_comm(curtime,count,array_of_requests,array_of_statuses,is_test,flag);
    case 1:
      return discretization::path::complete_comm(curtime,count,array_of_requests,array_of_statuses);
    case 2:
      return skeleton::path::complete_comm(curtime,count,array_of_requests,array_of_statuses);
  }
}

void propagate(MPI_Comm comm){
  switch (mechanism){
    case 0:
      decomposition::_MPI_Barrier.comm = comm;
      decomposition::_MPI_Barrier.partner1 = -1;
      decomposition::_MPI_Barrier.partner2 = -1;
      decomposition::path::propagate(decomposition::_MPI_Barrier);
      break;
    case 1:
      // Do nothing: final reduction is performed in 'discretization::final_accumulate'
      break;
    case 2:
      if (skeleton::skeleton_type==0){
        skeleton::_MPI_Barrier.comm = comm;
        skeleton::_MPI_Barrier.partner1 = -1;
        skeleton::_MPI_Barrier.partner2 = -1;
        skeleton::path::propagate_kernels(skeleton::_MPI_Barrier);
      }
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
    case 2:
      skeleton::volumetric::collect(comm);
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
      if (discretization::_MPI_Barrier.aggregate_comp_kernels){
        discretization::path::comp_state_aggregation(discretization::_MPI_Barrier);
      }
      if (discretization::_MPI_Barrier.aggregate_comm_kernels){
        discretization::path::comm_state_aggregation(discretization::_MPI_Barrier);
      }
      break;
    case 2:
      skeleton::final_accumulate(comm,last_time);
      skeleton::_MPI_Barrier.comm = comm;
      skeleton::_MPI_Barrier.partner1 = -1;
      skeleton::_MPI_Barrier.partner2 = -1;
      skeleton::path::propagate_kernels(skeleton::_MPI_Barrier);
      break;
  }
}

void set_reference_values(){
  switch (mechanism){
    case 0:
      break;
    case 1:
      discretization::reference_initiate();
      break;
    case 2:
      break;
  }
}

void save_reference_values(){
  switch (mechanism){
    case 0:
      discretization::cp_costs_ref[0] = decomposition::cp_costs[decomposition::num_cp_measures-5];
      discretization::cp_costs_ref[1] = decomposition::cp_costs[decomposition::num_cp_measures-3];
      discretization::cp_costs_ref[2] = decomposition::cp_costs[decomposition::num_cp_measures-2];
      discretization::cp_costs_ref[3] = decomposition::cp_costs[decomposition::num_cp_measures-1];
      discretization::max_pp_costs_ref[0] = decomposition::max_pp_costs[decomposition::num_pp_measures-5];
      discretization::max_pp_costs_ref[1] = decomposition::max_pp_costs[decomposition::num_pp_measures-3];
      discretization::max_pp_costs_ref[2] = decomposition::max_pp_costs[decomposition::num_pp_measures-2];
      discretization::max_pp_costs_ref[3] = decomposition::max_pp_costs[decomposition::num_pp_measures-1];
      discretization::vol_costs_ref[0] = decomposition::vol_costs[decomposition::num_vol_measures-5];
      discretization::vol_costs_ref[1] = decomposition::vol_costs[decomposition::num_vol_measures-3];
      discretization::vol_costs_ref[2] = decomposition::vol_costs[decomposition::num_vol_measures-2];
      discretization::vol_costs_ref[3] = decomposition::vol_costs[decomposition::num_vol_measures-1];
      break;
    case 1:
      discretization::reference_transfer();
      break;
    case 2:
      break;
  }
}

void init_symbol(std::vector<std::string>& symbols){
  switch (mechanism){
    case 0:
      decomposition::init_symbol(symbols);
      break;
    case 1:
      discretization::init_symbol(symbols);
      break;
    case 2:
      skeleton::init_symbol(symbols);
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
    case 2:
      skeleton::open_symbol(symbol,curtime);
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
    case 2:
      skeleton::close_symbol(symbol,curtime);
      break;
  }
}

void clear(int tag_count, int* distribution_tags){
  switch (mechanism){
    case 0:
      decomposition::clear();
      break;
    case 1:
      discretization::clear(tag_count,distribution_tags);
      break;
    case 2:
      skeleton::clear();
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
    case 2:
      skeleton::finalize();
      break;
  }
}

void write_file(int variantID, int print_mode, double overhead_time){
  if (std::getenv("CRITTER_VIZ_FILE") == NULL) return;
  switch (mechanism){
    case 0:
      decomposition::record::write_file(variantID,overhead_time);
      break;
    case 1:
      discretization::record::write_file(variantID,print_mode,overhead_time);
      break;
    case 2:
      skeleton::record::write_file(variantID,print_mode,overhead_time);
      break;
  }
}

void print(int variantID, int print_mode, double overhead_time){
  switch (mechanism){
    case 0:
      decomposition::record::print(variantID,overhead_time);
      break;
    case 1:
      discretization::record::print(variantID,print_mode,overhead_time);
      break;
    case 2:
      skeleton::record::print(variantID,print_mode,overhead_time);
      break;
  }
}

int get_critical_path_costs(){
  switch (mechanism){
    case 0:
      return decomposition::num_cp_measures;
    case 1:
      return discretization::num_cp_measures;
    case 2:
      return skeleton::num_cp_measures;
  }
  assert(-1);
  return -1;
}
void get_critical_path_costs(float* costs){
  switch (mechanism){
    case 0:
      std::memcpy(costs,&decomposition::cp_costs[0],sizeof(float)*decomposition::num_cp_measures);
      break;
    case 1:
      std::memcpy(costs,&discretization::cp_costs[0],sizeof(float)*discretization::num_cp_measures);
      break;
    case 2:
      std::memcpy(costs,&skeleton::cp_costs[0],sizeof(float)*skeleton::num_cp_measures);
      break;
  }
  return;
}
int get_max_per_process_costs(){
  switch (mechanism){
    case 0:
      return decomposition::num_pp_measures;
    case 1:
      return discretization::num_pp_measures;
    case 2:
      return skeleton::num_pp_measures;
  }
  assert(-1);
  return -1;
}
void get_max_per_process_costs(float* costs){
  switch (mechanism){
    case 0:
      std::memcpy(costs,&decomposition::max_pp_costs[0],sizeof(float)*decomposition::num_pp_measures);
      break;
    case 1:
      std::memcpy(costs,&discretization::max_pp_costs[0],sizeof(float)*discretization::num_pp_measures);
      break;
    case 2:
      std::memcpy(costs,&skeleton::max_pp_costs[0],sizeof(float)*skeleton::num_pp_measures);
      break;
  }
  return;
}
int get_volumetric_costs(){
  switch (mechanism){
    case 0:
      return decomposition::num_vol_measures;
    case 1:
      return discretization::num_vol_measures;
    case 2:
      return skeleton::num_vol_measures;
  }
  assert(-1);
  return -1;
}
void get_volumetric_costs(float* costs){
  switch (mechanism){
    case 0:
      std::memcpy(costs,&decomposition::vol_costs[0],sizeof(float)*decomposition::num_vol_measures);
      break;
    case 1:
      std::memcpy(costs,&discretization::vol_costs[0],sizeof(float)*discretization::num_vol_measures);
      break;
    case 2:
      std::memcpy(costs,&skeleton::vol_costs[0],sizeof(float)*skeleton::num_vol_measures);
      break;
  }
  return;
}

}
}
