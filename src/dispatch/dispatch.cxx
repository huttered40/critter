#include "dispatch.h"
#include "../profile/container/comm_tracker.h"
#include "../profile/util/util.h"
#include "../profile/volumetric/volumetric.h"
#include "../profile/path/path.h"
#include "../profile/record/record.h"
#include "../accelerate/container/comm_tracker.h"
#include "../accelerate/util/util.h"
#include "../accelerate/volumetric/volumetric.h"
#include "../accelerate/path/path.h"
#include "../accelerate/record/record.h"
#include "../skeletonize/container/comm_tracker.h"
#include "../skeletonize/util/util.h"
#include "../skeletonize/volumetric/volumetric.h"
#include "../skeletonize/path/path.h"
#include "../skeletonize/record/record.h"

namespace critter{
namespace internal{

void allocate(MPI_Comm comm){
  switch (mechanism){
    case 0:
      profile::allocate(comm);
      break;
    case 1:
      accelerate::allocate(comm);
      break;
    case 2:
      skeletonize::allocate(comm);
      break;
  }
}

void reset(bool schedule_kernels_override, bool force_steady_statistical_data_overide){
  switch (mechanism){
    case 0:
      profile::reset();
      break;
    case 1:
      accelerate::reset(schedule_kernels_override,force_steady_statistical_data_overide);
      break;
    case 2:
      skeletonize::reset();
      break;
  }
}

void exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm){
  switch (mechanism){
    case 0:
      //profile::path::exchange_communicators(oldcomm,newcomm);
      break;
    case 1:
      accelerate::path::exchange_communicators(oldcomm,newcomm);
      break;
    case 2:
      skeletonize::path::exchange_communicators(oldcomm,newcomm);
      break;
  }
}

bool initiate_comp(size_t id, volatile double curtime, double flop_count,  int param1, int param2, int param3, int param4, int param5){
  bool schedule_decision;
  switch (mechanism){
    case 0:
      schedule_decision = profile::path::initiate_comp(id,curtime,flop_count,param1,param2,param3,param4,param5);
      break;
    case 1:
      schedule_decision = accelerate::path::initiate_comp(id,curtime,flop_count,param1,param2,param3,param4,param5);
      break;
    case 2:
      schedule_decision = skeletonize::path::initiate_comp(id,curtime,flop_count,param1,param2,param3,param4,param5);
      break;
  }
  return schedule_decision;
}

void complete_comp(double errtime, size_t id, double flop_count,  int param1, int param2, int param3, int param4, int param5){
  switch (mechanism){
    case 0:
      profile::path::complete_comp(errtime,id,flop_count,param1,param2,param3,param4,param5);
      break;
    case 1:
      accelerate::path::complete_comp(errtime,id,flop_count,param1,param2,param3,param4,param5);
      break;
    case 2:
      skeletonize::path::complete_comp(errtime,id,flop_count,param1,param2,param3,param4,param5);
      break;
  }
}

bool initiate_comm(size_t id, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
              bool is_sender, int partner1, int user_tag1, int partner2, int user_tag2){
  bool schedule_decision;
  switch (mechanism){
    case 0:
      schedule_decision = profile::path::initiate_comm(*(profile::blocking*)profile::list[id],curtime,nelem,t,cm,is_sender,partner1,user_tag1,partner2,user_tag2);
      break;
    case 1:
      schedule_decision = accelerate::path::initiate_comm(*(accelerate::blocking*)accelerate::list[id],curtime,nelem,t,cm,is_sender,partner1,user_tag1,partner2,user_tag2);
      break;
    case 2:
      schedule_decision = skeletonize::path::initiate_comm(*(skeletonize::blocking*)skeletonize::list[id],curtime,nelem,t,cm,is_sender,partner1,user_tag1,partner2,user_tag2);
      break;
  }
  return schedule_decision;
}

bool inspect_comm(size_t id, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
              bool is_sender, int partner, int user_tag){
  bool schedule_decision;
  switch (mechanism){
    case 0:
      schedule_decision = profile::path::initiate_comm(*(profile::nonblocking*)profile::list[id],curtime,nelem,t,cm,is_sender,partner,user_tag);
      break;
    case 1:
      schedule_decision = accelerate::path::initiate_comm(*(accelerate::nonblocking*)accelerate::list[id],curtime,nelem,t,cm,is_sender,partner,user_tag);
      break;
    case 2:
      schedule_decision = skeletonize::path::initiate_comm(*(skeletonize::nonblocking*)skeletonize::list[id],curtime,nelem,t,cm,is_sender,partner,user_tag);
      break;
  }
  return schedule_decision;
}

void initiate_comm(size_t id, volatile double itime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, MPI_Request* request, bool is_sender, int partner, int user_tag){
  switch (mechanism){
    case 0:
      profile::path::initiate_comm(*(profile::nonblocking*)profile::list[id],itime,nelem,t,cm,request,is_sender,partner,user_tag);
      break;
    case 1:
      accelerate::path::initiate_comm(*(accelerate::nonblocking*)accelerate::list[id],itime,nelem,t,cm,request,is_sender,partner,user_tag);
      break;
    case 2:
      skeletonize::path::initiate_comm(*(skeletonize::nonblocking*)skeletonize::list[id],itime,nelem,t,cm,request,is_sender,partner,user_tag);
      break;
  }
}

void complete_comm(size_t id){
  switch (mechanism){
    case 0:
      profile::path::complete_comm(*(profile::blocking*)profile::list[id]);
      break;
    case 1:
      accelerate::path::complete_comm(*(accelerate::blocking*)accelerate::list[id]);
      break;
    case 2:
      skeletonize::path::complete_comm(*(skeletonize::blocking*)skeletonize::list[id]);
      break;
  }
}

int complete_comm(double curtime, MPI_Request* request, MPI_Status* status, int is_test, int* flag){
  switch (mechanism){
    case 0:
      return profile::path::complete_comm(curtime,request,status,is_test,flag);
    case 1:
      return accelerate::path::complete_comm(curtime,request,status);
    case 2:
      return skeletonize::path::complete_comm(curtime,request,status);
  }
}

int complete_comm(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status, int is_test, int* flag){
  switch (mechanism){
    case 0:
      return profile::path::complete_comm(curtime,count,array_of_requests,indx,status,is_test,flag);
    case 1:
      return accelerate::path::complete_comm(curtime,count,array_of_requests,indx,status);
    case 2:
      return skeletonize::path::complete_comm(curtime,count,array_of_requests,indx,status);
  }
}

int complete_comm(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[], int is_test){
  switch (mechanism){
    case 0:
      return profile::path::complete_comm(curtime,incount,array_of_requests,outcount,array_of_indices,array_of_statuses,is_test);
    case 1:
      return accelerate::path::complete_comm(curtime,incount,array_of_requests,outcount,array_of_indices,array_of_statuses);
    case 2:
      return skeletonize::path::complete_comm(curtime,incount,array_of_requests,outcount,array_of_indices,array_of_statuses);
  }
}

int complete_comm(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[], int is_test, int* flag){
  switch (mechanism){
    case 0:
      return profile::path::complete_comm(curtime,count,array_of_requests,array_of_statuses,is_test,flag);
    case 1:
      return accelerate::path::complete_comm(curtime,count,array_of_requests,array_of_statuses);
    case 2:
      return skeletonize::path::complete_comm(curtime,count,array_of_requests,array_of_statuses);
  }
}

void propagate(MPI_Comm comm){
  switch (mechanism){
    case 0:
      profile::_MPI_Barrier.comm = comm;
      profile::_MPI_Barrier.partner1 = -1;
      profile::_MPI_Barrier.partner2 = -1;
      profile::path::propagate(profile::_MPI_Barrier);
      break;
    case 1:
      // Do nothing: final reduction is performed in 'accelerate::final_accumulate'
      break;
    case 2:
      if (skeletonize::skeleton_type==0){
        skeletonize::_MPI_Barrier.comm = comm;
        skeletonize::_MPI_Barrier.partner1 = -1;
        skeletonize::_MPI_Barrier.partner2 = -1;
        skeletonize::path::propagate_kernels(skeletonize::_MPI_Barrier);
      }
      break;
  }
}

void collect(MPI_Comm comm){
  switch (mechanism){
    case 0:
      profile::volumetric::collect(comm);
      break;
    case 1:
      accelerate::volumetric::collect(comm);
      break;
    case 2:
      skeletonize::volumetric::collect(comm);
      break;
  }
}

void final_accumulate(MPI_Comm comm, double last_time){
  switch (mechanism){
    case 0:
      profile::final_accumulate(comm,last_time);
      break;
    case 1:
      accelerate::final_accumulate(comm,last_time);
      if (accelerate::_MPI_Barrier.aggregate_comp_kernels){
        accelerate::path::comp_state_aggregation(accelerate::_MPI_Barrier);
      }
      if (accelerate::_MPI_Barrier.aggregate_comm_kernels){
        accelerate::path::comm_state_aggregation(accelerate::_MPI_Barrier);
      }
      break;
    case 2:
      skeletonize::final_accumulate(comm,last_time);
      skeletonize::_MPI_Barrier.comm = comm;
      skeletonize::_MPI_Barrier.partner1 = -1;
      skeletonize::_MPI_Barrier.partner2 = -1;
      skeletonize::path::propagate_kernels(skeletonize::_MPI_Barrier);
      break;
  }
}

void set_reference_values(){
  switch (mechanism){
    case 0:
      break;
    case 1:
      accelerate::reference_initiate();
      break;
    case 2:
      break;
  }
}

void save_reference_values(){
  switch (mechanism){
    case 0:
      accelerate::cp_costs_ref[0] = profile::cp_costs[profile::num_cp_measures-5];
      accelerate::cp_costs_ref[1] = profile::cp_costs[profile::num_cp_measures-3];
      accelerate::cp_costs_ref[2] = profile::cp_costs[profile::num_cp_measures-2];
      accelerate::cp_costs_ref[3] = profile::cp_costs[profile::num_cp_measures-1];
      accelerate::max_pp_costs_ref[0] = profile::max_pp_costs[profile::num_pp_measures-5];
      accelerate::max_pp_costs_ref[1] = profile::max_pp_costs[profile::num_pp_measures-3];
      accelerate::max_pp_costs_ref[2] = profile::max_pp_costs[profile::num_pp_measures-2];
      accelerate::max_pp_costs_ref[3] = profile::max_pp_costs[profile::num_pp_measures-1];
      accelerate::vol_costs_ref[0] = profile::vol_costs[profile::num_vol_measures-5];
      accelerate::vol_costs_ref[1] = profile::vol_costs[profile::num_vol_measures-3];
      accelerate::vol_costs_ref[2] = profile::vol_costs[profile::num_vol_measures-2];
      accelerate::vol_costs_ref[3] = profile::vol_costs[profile::num_vol_measures-1];
      break;
    case 1:
      accelerate::reference_transfer();
      break;
    case 2:
      break;
  }
}

void init_symbol(std::vector<std::string>& symbols){
  switch (mechanism){
    case 0:
      profile::init_symbol(symbols);
      break;
    case 1:
      accelerate::init_symbol(symbols);
      break;
    case 2:
      skeletonize::init_symbol(symbols);
      break;
  }
}

void open_symbol(const char* symbol, double curtime){
  switch (mechanism){
    case 0:
      profile::open_symbol(symbol,curtime);
      break;
    case 1:
      accelerate::open_symbol(symbol,curtime);
      break;
    case 2:
      skeletonize::open_symbol(symbol,curtime);
      break;
  }
}

void close_symbol(const char* symbol, double curtime){
  switch (mechanism){
    case 0:
      profile::close_symbol(symbol,curtime);
      break;
    case 1:
      accelerate::close_symbol(symbol,curtime);
      break;
    case 2:
      skeletonize::close_symbol(symbol,curtime);
      break;
  }
}

void clear(int tag_count, int* distribution_tags){
  switch (mechanism){
    case 0:
      profile::clear();
      break;
    case 1:
      accelerate::clear(tag_count,distribution_tags);
      break;
    case 2:
      skeletonize::clear();
      break;
  }
}

void write_file(int variantID, int print_mode, double overhead_time){
  if (std::getenv("CRITTER_VIZ_FILE") == NULL) return;
  switch (mechanism){
    case 0:
      profile::record::write_file(variantID,overhead_time);
      break;
    case 1:
      accelerate::record::write_file(variantID,print_mode,overhead_time);
      break;
    case 2:
      skeletonize::record::write_file(variantID,print_mode,overhead_time);
      break;
  }
}

void print(int variantID, int print_mode, double overhead_time){
  switch (mechanism){
    case 0:
      profile::record::print(variantID,overhead_time);
      break;
    case 1:
      accelerate::record::print(variantID,print_mode,overhead_time);
      break;
    case 2:
      skeletonize::record::print(variantID,print_mode,overhead_time);
      break;
  }
}

int get_critical_path_costs(){
  switch (mechanism){
    case 0:
      return profile::num_cp_measures;
    case 1:
      return accelerate::num_cp_measures;
    case 2:
      return skeletonize::num_cp_measures;
  }
  assert(-1);
  return -1;
}
void get_critical_path_costs(float* costs){
  switch (mechanism){
    case 0:
      std::memcpy(costs,&profile::cp_costs[0],sizeof(float)*profile::num_cp_measures);
      break;
    case 1:
      std::memcpy(costs,&accelerate::cp_costs[0],sizeof(float)*accelerate::num_cp_measures);
      break;
    case 2:
      std::memcpy(costs,&skeletonize::cp_costs[0],sizeof(float)*skeletonize::num_cp_measures);
      break;
  }
  return;
}
int get_max_per_process_costs(){
  switch (mechanism){
    case 0:
      return profile::num_pp_measures;
    case 1:
      return accelerate::num_pp_measures;
    case 2:
      return skeletonize::num_pp_measures;
  }
  assert(-1);
  return -1;
}
void get_max_per_process_costs(float* costs){
  switch (mechanism){
    case 0:
      std::memcpy(costs,&profile::max_pp_costs[0],sizeof(float)*profile::num_pp_measures);
      break;
    case 1:
      std::memcpy(costs,&accelerate::max_pp_costs[0],sizeof(float)*accelerate::num_pp_measures);
      break;
    case 2:
      std::memcpy(costs,&skeletonize::max_pp_costs[0],sizeof(float)*skeletonize::num_pp_measures);
      break;
  }
  return;
}
int get_volumetric_costs(){
  switch (mechanism){
    case 0:
      return profile::num_vol_measures;
    case 1:
      return accelerate::num_vol_measures;
    case 2:
      return skeletonize::num_vol_measures;
  }
  assert(-1);
  return -1;
}
void get_volumetric_costs(float* costs){
  switch (mechanism){
    case 0:
      std::memcpy(costs,&profile::vol_costs[0],sizeof(float)*profile::num_vol_measures);
      break;
    case 1:
      std::memcpy(costs,&accelerate::vol_costs[0],sizeof(float)*accelerate::num_vol_measures);
      break;
    case 2:
      std::memcpy(costs,&skeletonize::vol_costs[0],sizeof(float)*skeletonize::num_vol_measures);
      break;
  }
  return;
}

}
}
