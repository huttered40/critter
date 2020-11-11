#include "comm_tracker.h"
#include "../util/util.h"

namespace critter{
namespace internal{
namespace decomposition{

std::map<MPI_Request,nonblocking*> internal_comm_track;
blocking _MPI_Barrier("MPI_Barrier",0, 
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,0.);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),0.);}
                        );
blocking _MPI_Bcast("MPI_Bcast",1,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n);}
                      );
blocking _MPI_Reduce("MPI_Reduce",2, 
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n);}
                       );
blocking _MPI_Allreduce("MPI_Allreduce",3,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);}, 
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n);}
                          );
blocking _MPI_Gather("MPI_Gather",4,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);}
                       );
blocking _MPI_Allgather("MPI_Allgather",5,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);}
                          );
blocking _MPI_Scatter("MPI_Scatter",6,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);}
                        );
blocking _MPI_Reduce_scatter("MPI_Reduce_scatter",7,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);}
                               );
blocking _MPI_Alltoall("MPI_Alltoall",8,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n);}
                         );
blocking _MPI_Gatherv("MPI_Gatherv",9,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);}
                        );
blocking _MPI_Allgatherv("MPI_Allgatherv",10,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);}
                           );
blocking _MPI_Scatterv("MPI_Scatterv",11,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);}
                         );
blocking _MPI_Alltoallv("MPI_Alltoallv",12,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                           return std::pair<double,double>(log2((double)p),log2((double)p)*n);}
                          );
blocking _MPI_Sendrecv("MPI_Sendrecv",13,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);}
                      );
blocking _MPI_Sendrecv_replace("MPI_Sendrecv_replace",14,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);}
                              );
blocking _MPI_Ssend("MPI_Ssend",15,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);}
                      );
blocking _MPI_Send("MPI_Send",16,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);}
                  );
blocking _MPI_Recv("MPI_Recv",17,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);}
                  );
nonblocking _MPI_Isend("MPI_Isend",18,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);}
                      );
nonblocking _MPI_Irecv("MPI_Irecv",19,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);}
                      );
nonblocking _MPI_Ibcast("MPI_Ibcast",20,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n);}
                       );
nonblocking _MPI_Iallreduce("MPI_Iallreduce",21,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n);}
                           );
nonblocking _MPI_Ireduce("MPI_Ireduce",22,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n);}
                        );
nonblocking _MPI_Igather("MPI_Igather",23,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);}
                        );
nonblocking _MPI_Igatherv("MPI_Igatherv",24,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);}
                         );
nonblocking _MPI_Iallgather("MPI_Iallgather",25,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);}
                           );
nonblocking _MPI_Iallgatherv("MPI_Iallgatherv",26,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);}
                            );
nonblocking _MPI_Iscatter("MPI_Iscatter",27,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);}
                         );
nonblocking _MPI_Iscatterv("MPI_Iscatterv",28,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);}
                          );
nonblocking _MPI_Ireduce_scatter("MPI_Ireduce_scatter",29,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);}
                                );
nonblocking _MPI_Ialltoall("MPI_Ialltoall",30,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n);}
                          );
nonblocking _MPI_Ialltoallv("MPI_Ialltoallv",31,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n);}
                           );
blocking _MPI_Bsend("MPI_Bsend",32,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);}
                      );

comm_tracker* list[list_size] = {
        &_MPI_Barrier,
        &_MPI_Bcast,
        &_MPI_Reduce,
        &_MPI_Allreduce,
        &_MPI_Gather,
        &_MPI_Allgather,
        &_MPI_Scatter,
        &_MPI_Reduce_scatter,
        &_MPI_Alltoall,
        &_MPI_Gatherv,
        &_MPI_Allgatherv,
        &_MPI_Scatterv,
        &_MPI_Alltoallv,
        &_MPI_Sendrecv,
        &_MPI_Sendrecv_replace,
        &_MPI_Ssend,
        &_MPI_Send,
        &_MPI_Recv,
        &_MPI_Isend,
        &_MPI_Irecv,
        &_MPI_Ibcast,
        &_MPI_Iallreduce,
        &_MPI_Ireduce,
        &_MPI_Igather,
        &_MPI_Igatherv,
        &_MPI_Iallgather,
        &_MPI_Iallgatherv,
        &_MPI_Iscatter,
        &_MPI_Iscatterv,
        &_MPI_Ireduce_scatter,
        &_MPI_Ialltoall,
        &_MPI_Ialltoallv,
        &_MPI_Bsend};

void comm_tracker::init(){
  this->set_cost_pointers();
  this->start_time  = -1.;
  this->comp_time   = 0.;
  this->partner1    = -1;
  this->partner2    = -1;
}

void comm_tracker::set_cost_pointers(){
  size_t volume_costs_idx        = num_volume_measures+this->tag*num_tracker_volume_measures;
  this->my_wrd_count             = cost_model_size>0 ? &volume_costs[volume_costs_idx] : &scratch_pad;
  this->my_msg_count             = cost_model_size>0 ? &volume_costs[volume_costs_idx+cost_model_size] : &scratch_pad;
  this->my_comm_time             = &volume_costs[volume_costs_idx+2*cost_model_size];
  this->my_synch_time            = &volume_costs[volume_costs_idx+2*cost_model_size+1];
  if (comm_path_select_size>0){
    size_t critical_path_costs_idx   = num_critical_path_measures+this->tag*comm_path_select_size*num_tracker_critical_path_measures;
    this->critical_path_wrd_count    = cost_model_size>0 ? &critical_path_costs[critical_path_costs_idx] : &scratch_pad;
    this->critical_path_msg_count    = cost_model_size>0 ? &critical_path_costs[critical_path_costs_idx+cost_model_size*comm_path_select_size] : &scratch_pad;
    this->critical_path_comm_time    = &critical_path_costs[critical_path_costs_idx+2*comm_path_select_size*cost_model_size];
    this->critical_path_synch_time   = &critical_path_costs[critical_path_costs_idx+2*comm_path_select_size*cost_model_size+comm_path_select_size];
  } else{
    this->critical_path_wrd_count    = &scratch_pad;
    this->critical_path_msg_count    = &scratch_pad;
    this->critical_path_comm_time    = &scratch_pad;
    this->critical_path_synch_time   = &scratch_pad;
  }
}

void comm_tracker::set_volume_costs(){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if (*this->my_comm_time != 0){
    std::vector<double> vec(num_tracker_volume_measures);
    int save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]=='1'){
        vec[2*save] = *(this->my_wrd_count+save);
        vec[2*save+1] = *(this->my_msg_count+save);
        save++;
      }
    }
    vec[2*cost_model_size] = *this->my_comm_time;
    vec[2*cost_model_size+1] = *this->my_synch_time;
    save_info[this->name] = std::move(vec);
  }
}

void comm_tracker::set_header(){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if (*this->my_comm_time != 0){
    std::vector<double> vec(1);
    save_info[this->name] = std::move(vec);
  }
}

void comm_tracker::set_critical_path_costs(size_t idx){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if ((*this->my_comm_time != 0) && (comm_path_select_size>0)){
    std::vector<double> vec(num_tracker_critical_path_measures);
    int save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]=='1'){
        vec[2*save] = *(this->critical_path_wrd_count+idx+save*comm_path_select_size);
        vec[2*save+1] = *(this->critical_path_msg_count+idx+save*comm_path_select_size);
        save++;
      }
    }
    vec[num_tracker_critical_path_measures-2] = *(this->critical_path_comm_time+idx);
    vec[num_tracker_critical_path_measures-1] = *(this->critical_path_synch_time+idx);
    save_info[this->name] = std::move(vec);
  }
}

void comm_tracker::set_per_process_costs(size_t idx){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if ((*this->my_comm_time != 0) && (comm_path_select_size>0)){
    std::vector<double> vec(num_tracker_per_process_measures);
    int save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]=='1'){
        vec[2*save] = max_per_process_costs[num_per_process_measures+idx*(num_tracker_per_process_measures*list_size+2)+this->tag*num_tracker_per_process_measures+save];
        vec[2*save+1] = max_per_process_costs[num_per_process_measures+idx*(num_tracker_per_process_measures*list_size+2)+this->tag*num_tracker_per_process_measures+cost_model_size+save];
        save++;
      }
    }
    // For now, do not include idle time
    vec[num_tracker_per_process_measures-2] = max_per_process_costs[num_per_process_measures+idx*(num_tracker_per_process_measures*list_size+2)+this->tag*num_tracker_per_process_measures+num_tracker_per_process_measures-2];
    vec[num_tracker_per_process_measures-1] = max_per_process_costs[num_per_process_measures+idx*(num_tracker_per_process_measures*list_size+2)+this->tag*num_tracker_per_process_measures+num_tracker_per_process_measures-1];
    save_info[this->name] = std::move(vec);
  }
}

blocking::blocking(std::string name_, int tag, std::function<std::pair<double,double>(int64_t,int)> cost_func_bsp,
                                             std::function<std::pair<double,double>(int64_t,int)> cost_func_alphabeta){
  this->cost_func_bsp              = cost_func_bsp;
  this->cost_func_alphabeta = cost_func_alphabeta;
  this->name = std::move(name_);
  this->tag = tag;
  this->is_sender = tag < 17 ? true : false;
}

blocking::blocking(blocking const& t){
  this->cost_func_bsp                 = t.cost_func_bsp;
  this->cost_func_alphabeta = t.cost_func_alphabeta;
  this->name = t.name;
  this->tag = t.tag;
  this->is_sender = t.is_sender;
}

nonblocking::nonblocking(std::string name_, int tag, std::function<std::pair<double,double>(int64_t,int)> cost_func_bsp,
                                             std::function<std::pair<double,double>(int64_t,int)> cost_func_alphabeta){
  this->cost_func_bsp                 = cost_func_bsp;
  this->cost_func_alphabeta = cost_func_alphabeta;
  this->name = std::move(name_);
  this->tag = tag;
  this->is_sender = tag==18 ? true : false;
}

nonblocking::nonblocking(nonblocking const& t){
  this->cost_func_bsp                 = t.cost_func_bsp;
  this->cost_func_alphabeta = t.cost_func_alphabeta;
  this->name = t.name;
  this->tag = t.tag;
  this->is_sender = t.is_sender;
}

}
}
}
