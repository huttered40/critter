#include "comm_tracker.h"
#include "../util/util.h"

namespace critter{
namespace internal{
namespace decomposition{

std::map<MPI_Request,nonblocking_info> nonblocking_internal_info;

blocking _MPI_Barrier("MPI_Barrier",0, 
                       [](int64_t n, int p){
                       return std::pair<float,float>(1.,0.);},
                       [](int64_t n, int p){
                       return std::pair<float,float>(log2((float)p),0.);}
                     );
blocking _MPI_Bcast("MPI_Bcast",1,
                     [](int64_t n, int p){
                     return std::pair<float,float>(1.,n);},
                     [](int64_t n, int p){
                     return std::pair<float,float>(2.*log2((float)p),2.*n);}
                   );
blocking _MPI_Reduce("MPI_Reduce",2, 
                      [](int64_t n, int p){
                      return std::pair<float,float>(1.,n);},
                      [](int64_t n, int p){
                      return std::pair<float,float>(2.*log2((float)p),2.*n);}
                    );
blocking _MPI_Allreduce("MPI_Allreduce",3,
                         [](int64_t n, int p){
                         return std::pair<float,float>(1.,n);}, 
                         [](int64_t n, int p){
                         return std::pair<float,float>(2.*log2((float)p),2.*n);}
                       );
blocking _MPI_Gather("MPI_Gather",4,
                      [](int64_t n, int p){
                      return std::pair<float,float>(1.,n);},
                      [](int64_t n, int p){
                      return std::pair<float,float>(log2((float)p),n);}
                    );
blocking _MPI_Allgather("MPI_Allgather",5,
                         [](int64_t n, int p){
                         return std::pair<float,float>(1.,n);},
                         [](int64_t n, int p){
                         return std::pair<float,float>(log2((float)p),n);}
                       );
blocking _MPI_Scatter("MPI_Scatter",6,
                       [](int64_t n, int p){
                       return std::pair<float,float>(1.,n);},
                       [](int64_t n, int p){
                       return std::pair<float,float>(log2((float)p),n);}
                     );
blocking _MPI_Reduce_scatter("MPI_Reduce_scatter",7,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(log2((float)p),n);}
                            );
blocking _MPI_Alltoall("MPI_Alltoall",8,
              [](int64_t n, int p){
              return std::pair<float,float>(1.,n);},
              [](int64_t n, int p){
              return std::pair<float,float>(log2((float)p),log2((float)p)*n);}
                      );
blocking _MPI_Gatherv("MPI_Gatherv",9,
                       [](int64_t n, int p){
                       return std::pair<float,float>(1.,n);},
                       [](int64_t n, int p){
                       return std::pair<float,float>(log2((float)p),n);}
                     );
blocking _MPI_Allgatherv("MPI_Allgatherv",10,
                          [](int64_t n, int p){
                          return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                          return std::pair<float,float>(log2((float)p),n);}
                        );
blocking _MPI_Scatterv("MPI_Scatterv",11,
                        [](int64_t n, int p){
                        return std::pair<float,float>(1.,n);},
                        [](int64_t n, int p){
                        return std::pair<float,float>(log2((float)p),n);}
                      );
blocking _MPI_Alltoallv("MPI_Alltoallv",12,
              [](int64_t n, int p){
              return std::pair<float,float>(1.,n);},
              [](int64_t n, int p){
              return std::pair<float,float>(log2((float)p),log2((float)p)*n);}
                       );
blocking _MPI_Sendrecv("MPI_Sendrecv",13,
                        [](int64_t n, int p){
                        return std::pair<float,float>(1.,n);},
                        [](int64_t n, int p){
                        return std::pair<float,float>(1.,n);}
                      );
blocking _MPI_Sendrecv_replace("MPI_Sendrecv_replace",14,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);}
                              );
blocking _MPI_Ssend("MPI_Ssend",15,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);}
                   );
blocking _MPI_Bsend("MPI_Bsend",16,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);}
                   );
blocking _MPI_Send("MPI_Send",17,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);}
                  );
blocking _MPI_Recv("MPI_Recv",18,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);}
                  );
nonblocking _MPI_Isend("MPI_Isend",19,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);}
                      );
nonblocking _MPI_Irecv("MPI_Irecv",20,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);}
                      );
nonblocking _MPI_Ibcast("MPI_Ibcast",21,
                       [](int64_t n, int p){
                       return std::pair<float,float>(1.,n);},
                       [](int64_t n, int p){
                       return std::pair<float,float>(2.*log2((float)p),2.*n);}
                       );
nonblocking _MPI_Iallreduce("MPI_Iallreduce",22,
                       [](int64_t n, int p){
                       return std::pair<float,float>(1.,n);},
                       [](int64_t n, int p){
                       return std::pair<float,float>(2.*log2((float)p),2.*n);}
                           );
nonblocking _MPI_Ireduce("MPI_Ireduce",23,
                        [](int64_t n, int p){
                        return std::pair<float,float>(1.,n);},
                        [](int64_t n, int p){
                        return std::pair<float,float>(2.*log2((float)p),2.*n);}
                        );
nonblocking _MPI_Igather("MPI_Igather",24,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(log2((float)p),n);}
                        );
nonblocking _MPI_Igatherv("MPI_Igatherv",25,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(log2((float)p),n);}
                         );
nonblocking _MPI_Iallgather("MPI_Iallgather",26,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(log2((float)p),n);}
                           );
nonblocking _MPI_Iallgatherv("MPI_Iallgatherv",27,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(log2((float)p),n);}
                            );
nonblocking _MPI_Iscatter("MPI_Iscatter",28,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(log2((float)p),n);}
                         );
nonblocking _MPI_Iscatterv("MPI_Iscatterv",29,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(log2((float)p),n);}
                          );
nonblocking _MPI_Ireduce_scatter("MPI_Ireduce_scatter",30,
                          [](int64_t n, int p){
                            return std::pair<float,float>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<float,float>(log2((float)p),n);}
                                );
nonblocking _MPI_Ialltoall("MPI_Ialltoall",31,
              [](int64_t n, int p){
              return std::pair<float,float>(1.,n);},
              [](int64_t n, int p){
              return std::pair<float,float>(log2((float)p),log2((float)p)*n);}
                          );
nonblocking _MPI_Ialltoallv("MPI_Ialltoallv",32,
              [](int64_t n, int p){
              return std::pair<float,float>(1.,n);},
              [](int64_t n, int p){
              return std::pair<float,float>(log2((float)p),log2((float)p)*n);}
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
        &_MPI_Bsend,
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
        &_MPI_Ialltoallv};

void comm_tracker::init(){
  this->set_cost_pointers();
  this->start_time = -1.;
  this->comp_time = 0.;
  this->partner1 = -1;
  this->partner2 = -1;
  decomp_text_width = std::max(decomp_text_width,this->name.size());
}

void comm_tracker::set_cost_pointers(){
  this->my_wrd_count = &scratch_pad;
  this->my_msg_count = &scratch_pad;
  this->my_comm_time = &scratch_pad;
  this->my_synch_time = &scratch_pad;
  this->cp_wrd_count = &scratch_pad;
  this->cp_msg_count = &scratch_pad;
  this->cp_comm_time = &scratch_pad;
  this->cp_synch_time = &scratch_pad;
  if (path_decomposition<=1){
    size_t vol_costs_idx = num_vol_measures+this->tag*num_decomp_vol_measures;
    this->my_wrd_count = &vol_costs[vol_costs_idx];
    this->my_msg_count = &vol_costs[vol_costs_idx+1];
    this->my_comm_time = &vol_costs[vol_costs_idx+2];
    this->my_synch_time = &vol_costs[vol_costs_idx+3];
  }
  if (path_decomposition==1 && path_count>0){
    size_t cp_costs_idx = num_cp_measures;
           cp_costs_idx += this->tag*path_count*num_decomp_cp_measures;
    size_t j = 0;
    assert(num_decomp_cp_measures == path_measure_index.size());
    for (int i=0; i<path_measure_index.size(); i++){
      if (path_measure_index[i]==0) { this->cp_wrd_count = &cp_costs[cp_costs_idx+j*path_count]; j++; }
      else if (path_measure_index[i]==1) { this->cp_msg_count = &cp_costs[cp_costs_idx+j*path_count]; j++; }
      else if (path_measure_index[i]==2) { this->cp_comm_time = &cp_costs[cp_costs_idx+j*path_count]; j++; }
      else if (path_measure_index[i]==3) { this->cp_synch_time = &cp_costs[cp_costs_idx+j*path_count]; j++; }
      else assert(0);
    }
  }
}

void comm_tracker::set_vol_costs(std::map<std::string,std::vector<float>>& save_info){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if (*this->my_comm_time != 0 && path_decomposition<=1){
    std::vector<float> vec(num_decomp_vol_measures);
    vec[0] = *this->my_wrd_count;
    vec[1] = *this->my_msg_count;
    vec[2] = *this->my_comm_time;
    vec[3] = *this->my_synch_time;
    save_info[this->name] = std::move(vec);
  }
}

void comm_tracker::set_header(std::map<std::string,std::vector<float>>& save_info){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if (*this->my_comm_time != 0 && path_decomposition<=1){
    std::vector<float> vec(1);
    save_info[this->name] = std::move(vec);
  }
}

void comm_tracker::set_cp_costs(std::map<std::string,std::vector<float>>& save_info, size_t idx){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if ((*this->my_comm_time != 0 && path_decomposition==1 && path_count>0)){
    std::vector<float> vec(4,0);
    for (int i=0; i<num_decomp_cp_measures; i++){
      if (path_measure_index[i]==0) { vec[path_measure_index[i]] = *(this->cp_wrd_count+idx); }
      else if (path_measure_index[i]==1) { vec[path_measure_index[i]] = *(this->cp_msg_count+idx); }
      else if (path_measure_index[i]==2) { vec[path_measure_index[i]] = *(this->cp_comm_time+idx); }
      else if (path_measure_index[i]==3) { vec[path_measure_index[i]] = *(this->cp_synch_time+idx); }
      else assert(0);
    }
    save_info[this->name] = std::move(vec);
  }
}

void comm_tracker::set_pp_costs(std::map<std::string,std::vector<float>>& save_info, size_t idx){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if ((*this->my_comm_time != 0 && path_decomposition<=1 && path_count>0)){
    std::vector<float> vec(num_decomp_pp_measures);
    vec[0] = max_pp_costs[num_pp_measures+idx*(num_decomp_pp_measures*list_size+4)+this->tag*num_decomp_pp_measures];
    vec[1] = max_pp_costs[num_pp_measures+idx*(num_decomp_pp_measures*list_size+4)+this->tag*num_decomp_pp_measures+1];
    vec[2] = max_pp_costs[num_pp_measures+idx*(num_decomp_pp_measures*list_size+4)+this->tag*num_decomp_pp_measures+2];
    vec[3] = max_pp_costs[num_pp_measures+idx*(num_decomp_pp_measures*list_size+4)+this->tag*num_decomp_pp_measures+3];
    save_info[this->name] = std::move(vec);
  }
}

blocking::blocking(std::string name_, int tag, std::function<std::pair<float,float>(int64_t,int)> cost_func_bsp,
                   std::function<std::pair<float,float>(int64_t,int)> cost_func_alphabeta){
  this->cost_func_bsp = cost_func_bsp;
  this->cost_func_alphabeta = cost_func_alphabeta;
  this->name = std::move(name_);
  this->tag = tag;
  this->is_sender = tag < 18 ? true : false;
}

blocking::blocking(blocking const& t){
  this->cost_func_bsp = t.cost_func_bsp;
  this->cost_func_alphabeta = t.cost_func_alphabeta;
  this->name = t.name;
  this->tag = t.tag;
  this->is_sender = t.is_sender;
}

nonblocking::nonblocking(std::string name_, int tag, std::function<std::pair<float,float>(int64_t,int)> cost_func_bsp,
                         std::function<std::pair<float,float>(int64_t,int)> cost_func_alphabeta){
  this->cost_func_bsp = cost_func_bsp;
  this->cost_func_alphabeta = cost_func_alphabeta;
  this->name = std::move(name_);
  this->tag = tag;
  this->is_sender = (tag==19 || tag>20) ? true : false;
}

nonblocking::nonblocking(nonblocking const& t){
  this->cost_func_bsp = t.cost_func_bsp;
  this->cost_func_alphabeta = t.cost_func_alphabeta;
  this->name = t.name;
  this->tag = t.tag;
  this->is_sender = t.is_sender;
}

}
}
}
