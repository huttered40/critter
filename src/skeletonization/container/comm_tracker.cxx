#include "comm_tracker.h"

namespace critter{
namespace internal{
namespace skeletonization{

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
  this->start_time  = -1.;
  this->comp_time   = 0.;
  this->partner1    = -1;
  this->partner2    = -1;
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
