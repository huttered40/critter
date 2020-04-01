#include "critter.h"

namespace critter{
namespace internal{

void add_critical_path_data_op(int_int_double* in, int_int_double* inout, int* len, MPI_Datatype* dtype){
  int_int_double* invec = in;
  int_int_double* inoutvec = inout;
  for (int i=0; i<*len; i++){
    inoutvec[i].first = std::max(inoutvec[i].first,invec[i].first);
    inoutvec[i].second = std::max(inoutvec[i].second,invec[i].second);
    inoutvec[i].third = std::max(inoutvec[i].third,invec[i].third);
  }
}

void update_critical_path(double* in, double* inout, size_t len){
  if (breakdown_size > 0){
    size_t breakdown_idx=0;
    for (int i=0; i<num_critical_path_measures-2*cost_model_size; i++){
      if (breakdown[i]){
        decisions[breakdown_idx++] = inout[2*cost_model_size+i] > in[2*cost_model_size+i];
      }
    }
    for (int i=0; i<num_critical_path_measures; i++){
      inout[i] = std::max(inout[i],in[i]);
    }
    for (int i=num_critical_path_measures; i<len; i++){
      int idx = (i-num_critical_path_measures)%breakdown_size;
      inout[i] = (decisions[idx] ? inout[i] : in[i]);
    }
  } else{
    for (int i=0; i<num_critical_path_measures; i++){
      inout[i] = std::max(inout[i],in[i]);
    }
  }
}

void propagate_critical_path_op(double* in, double* inout, int* len, MPI_Datatype* dtype){
  update_critical_path(in,inout,static_cast<size_t>(*len));
}

synchronous _MPI_Barrier("MPI_Barrier",0, 
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,0.);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),0.);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,0.);}
                        );
synchronous _MPI_Bcast("MPI_Bcast",1,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(is_root ? p : 1.,is_root ? n*p : n);}
                      );
synchronous _MPI_Reduce("MPI_Reduce",2, 
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(is_root ? p : 1.,is_root ? n*p : n);}
                       );
synchronous _MPI_Allreduce("MPI_Allreduce",3,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);}, 
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n);}, 
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n*p);}
                          );
synchronous _MPI_Gather("MPI_Gather",4,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(is_root ? p : 1.,is_root ? n : (n*1.)/(p*1.));}
                       );
synchronous _MPI_Allgather("MPI_Allgather",5,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                          );
synchronous _MPI_Scatter("MPI_Scatter",6,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(is_root ? p : 1.,is_root ? n : (n*1.)/(p*1.));}
                        );
synchronous _MPI_Reduce_scatter("MPI_Reduce_scatter",7,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                               );
synchronous _MPI_Alltoall("MPI_Alltoall",8,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                         );
synchronous _MPI_Gatherv("MPI_Gatherv",9,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(is_root ? p : 1.,is_root ? n : (n*1.)/(p*1.));}
                        );
synchronous _MPI_Allgatherv("MPI_Allgatherv",10,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                           );
synchronous _MPI_Scatterv("MPI_Scatterv",11,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(is_root ? p : 1.,is_root ? n : (n*1.)/(p*1.));}
                         );
synchronous _MPI_Alltoallv("MPI_Alltoallv",12,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                           return std::pair<double,double>(log2((double)p),log2((double)p)*n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                          );
synchronous _MPI_Ssend("MPI_Ssend",13,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(1.,n);} 
                      );
p2p_type _MPI_Sendrecv("MPI_Sendrecv",14,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(1.,n);}
                      );
p2p_type _MPI_Sendrecv_replace("MPI_Sendrecv_replace",15,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(1.,n);}
                              );
p2p_type _MPI_Send("MPI_Send",16,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(1.,n);}
                  );
p2p_type _MPI_Recv("MPI_Recv",17,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(1.,n);}
                  );
nonblocking _MPI_Isend("MPI_Isend",18,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(0.,n);}
                      );
nonblocking _MPI_Irecv("MPI_Irecv",19,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(0.,n);}
                      );
nonblocking _MPI_Ibcast("MPI_Ibcast",20,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                       );
nonblocking _MPI_Iallreduce("MPI_Iallreduce",21,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                           );
nonblocking _MPI_Ireduce("MPI_Ireduce",22,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                        );
nonblocking _MPI_Igather("MPI_Igather",23,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                        );
nonblocking _MPI_Igatherv("MPI_Igatherv",24,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                         );
nonblocking _MPI_Iallgather("MPI_Iallgather",25,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                           );
nonblocking _MPI_Iallgatherv("MPI_Iallgatherv",26,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                            );
nonblocking _MPI_Iscatter("MPI_Iscatter",27,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                         );
nonblocking _MPI_Iscatterv("MPI_Iscatterv",28,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                          );
nonblocking _MPI_Ireduce_scatter("MPI_Ireduce_scatter",29,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                                );
nonblocking _MPI_Ialltoall("MPI_Ialltoall",30,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                          );
nonblocking _MPI_Ialltoallv("MPI_Ialltoallv",31,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n);},
                          [](int64_t n, int p, bool is_root){
                            return std::pair<double,double>(p,n);}
                           );

tracker* list[list_size] = {
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
        &_MPI_Ssend,
        &_MPI_Sendrecv,
        &_MPI_Sendrecv_replace,
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

std::string stream_name,file_name;
std::ofstream stream;
bool flag,is_world_root,is_first_iter,need_new_line;
size_t mode;

double computation_timer;
std::map<MPI_Request,std::pair<MPI_Request,bool>> internal_comm_info;
std::map<MPI_Request,std::pair<MPI_Comm,int>> internal_comm_comm;
std::map<MPI_Request,double*> internal_comm_message;
std::map<MPI_Request,std::pair<double,double>> internal_comm_data;
std::map<MPI_Request,nonblocking*> internal_comm_track;
bool decisions[breakdown_size];
std::array<double,critical_path_costs_size> critical_path_costs;
std::array<double,per_process_costs_size> max_per_process_costs;
std::array<double,volume_costs_size> volume_costs;
std::map<std::string,std::vector<double>> save_info;
double new_cs[critical_path_costs_size];
double scratch_pad;
std::vector<char> synch_pad_send;
std::vector<char> synch_pad_recv;
std::array<char,max_timer_name_length*max_num_symbols> symbol_pad;
std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_local_cp;
std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_global_cp;
std::array<double,(num_ftimer_measures*num_volume_measures+1)*max_num_symbols> symbol_timer_pad_pp;
std::array<double,(num_ftimer_measures*num_volume_measures+1)*max_num_symbols> symbol_timer_pad_vol;
std::unordered_map<std::string,ftimer> symbol_timers;
std::stack<std::string> symbol_stack;
std::array<std::string,max_num_symbols> symbol_order;
std::array<std::string,num_critical_path_measures> critical_path_measure_names;
double_int timer_info_sender[num_volume_measures];
double_int timer_info_receiver[num_volume_measures];
bool wait_id;

void tracker::init(){
  this->last_start_time  = -1.;
  this->save_comp_time   = 0.;
}

void tracker::set_cost_pointers(){
  size_t volume_costs_idx        = num_volume_measures+this->tag*num_tracker_volume_measures;
  this->my_wrd_count             = cost_model_size>0 ? &volume_costs[volume_costs_idx] : &scratch_pad;
  this->my_msg_count             = cost_model_size>0 ? &volume_costs[volume_costs_idx+cost_model_size] : &scratch_pad;
  this->my_bar_time              = &volume_costs[volume_costs_idx+2*cost_model_size];
  this->my_comm_time             = &volume_costs[volume_costs_idx+2*cost_model_size+1];
  this->my_synch_time            = &volume_costs[volume_costs_idx+2*cost_model_size+2];
  this->my_datamvt_time          = &volume_costs[volume_costs_idx+2*cost_model_size+3];
  if (breakdown_size>0){
    size_t critical_path_costs_idx   = num_critical_path_measures+this->tag*breakdown_size*num_tracker_critical_path_measures;
    this->critical_path_wrd_count    = cost_model_size>0 ? &critical_path_costs[critical_path_costs_idx] : &scratch_pad;
    this->critical_path_msg_count    = cost_model_size>0 ? &critical_path_costs[critical_path_costs_idx+cost_model_size*breakdown_size] : &scratch_pad;
    this->critical_path_comm_time    = &critical_path_costs[critical_path_costs_idx+2*breakdown_size*cost_model_size];
    this->critical_path_synch_time   = &critical_path_costs[critical_path_costs_idx+2*breakdown_size*cost_model_size+breakdown_size];
    this->critical_path_datamvt_time = &critical_path_costs[critical_path_costs_idx+2*breakdown_size*cost_model_size+2*breakdown_size];
  } else{
    this->critical_path_wrd_count    = &scratch_pad;
    this->critical_path_msg_count    = &scratch_pad;
    this->critical_path_comm_time    = &scratch_pad;
    this->critical_path_datamvt_time = &scratch_pad;
    this->critical_path_synch_time   = &scratch_pad;
  }
}

synchronous::synchronous(std::string name_, int tag, std::function<std::pair<double,double>(int64_t,int)> cost_func_simple,
                                             std::function<std::pair<double,double>(int64_t,int)> cost_func_alphabeta_butterfly,
                                             std::function<std::pair<double,double>(int64_t,int,bool)> cost_func_bsp){
  this->cost_func_simple              = cost_func_simple;
  this->cost_func_alphabeta_butterfly = cost_func_alphabeta_butterfly;
  this->cost_func_bsp                 = cost_func_bsp;
  this->name = std::move(name_);
  this->tag = tag;
  this->set_cost_pointers();
  this->init();
}

synchronous::synchronous(synchronous const& t){
  this->cost_func_simple              = t.cost_func_simple;
  this->cost_func_alphabeta_butterfly = t.cost_func_alphabeta_butterfly;
  this->cost_func_bsp                 = t.cost_func_bsp;
  this->name = t.name;
  this->tag = t.tag;
  this->set_cost_pointers();
  this->init();
}

blocking::blocking(std::string name_, int tag, std::function<std::pair<double,double>(int64_t,int)> cost_func_simple,
                                             std::function<std::pair<double,double>(int64_t,int)> cost_func_alphabeta_butterfly,
                                             std::function<std::pair<double,double>(int64_t,int,bool)> cost_func_bsp){
  this->cost_func_simple              = cost_func_simple;
  this->cost_func_alphabeta_butterfly = cost_func_alphabeta_butterfly;
  this->cost_func_bsp                 = cost_func_bsp;
  this->name = std::move(name_);
  this->tag = tag;
  this->set_cost_pointers();
  this->init();
}

blocking::blocking(blocking const& t){
  this->cost_func_simple              = t.cost_func_simple;
  this->cost_func_alphabeta_butterfly = t.cost_func_alphabeta_butterfly;
  this->cost_func_bsp                 = t.cost_func_bsp;
  this->name = t.name;
  this->tag = t.tag;
  this->set_cost_pointers();
  this->init();
}

nonblocking::nonblocking(std::string name_, int tag, std::function<std::pair<double,double>(int64_t,int)> cost_func_simple,
                                             std::function<std::pair<double,double>(int64_t,int)> cost_func_alphabeta_butterfly,
                                             std::function<std::pair<double,double>(int64_t,int,bool)> cost_func_bsp){
  this->cost_func_simple              = cost_func_simple;
  this->cost_func_alphabeta_butterfly = cost_func_alphabeta_butterfly;
  this->cost_func_bsp                 = cost_func_bsp;
  this->name = std::move(name_);
  this->tag = tag;
  this->set_cost_pointers();
  this->init();
}

nonblocking::nonblocking(nonblocking const& t){
  this->cost_func_simple              = t.cost_func_simple;
  this->cost_func_alphabeta_butterfly = t.cost_func_alphabeta_butterfly;
  this->cost_func_bsp                 = t.cost_func_bsp;
  this->name = t.name;
  this->tag = t.tag;
  this->set_cost_pointers();
  this->init();
}

void synchronous::start(volatile double curTime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, bool is_root, bool is_sender, int partner1, int partner2){

  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  this->save_comp_time    = curTime - computation_timer;
  critical_path_costs[num_critical_path_measures-2] += this->save_comp_time;	// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += this->save_comp_time;	// update critical path runtime
  volume_costs[num_volume_measures-2]        += this->save_comp_time;		// update local computation time
  volume_costs[num_volume_measures-1]        += this->save_comp_time;		// update local runtime
  for (size_t i=0; i<breakdown_size; i++){
    critical_path_costs[critical_path_costs_size-1-i] += this->save_comp_time;
  }
  if (mode == 2){
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[num_critical_path_measures-1] += (curTime - symbol_timers[symbol_stack.top()].start_timer.top());
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[num_critical_path_measures-2] += (curTime - symbol_timers[symbol_stack.top()].start_timer.top());
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-1] += (curTime - symbol_timers[symbol_stack.top()].start_timer.top());
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-2] += (curTime - symbol_timers[symbol_stack.top()].start_timer.top());
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-1] += (curTime - symbol_timers[symbol_stack.top()].start_timer.top());
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-2] += (curTime - symbol_timers[symbol_stack.top()].start_timer.top());
    //TODO: Need exclusive_measure buckets for tracking the per-process inclusive measures???
  }

  int el_size,p;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(cm, &p);
  this->last_nbytes = nbytes;
  this->last_p = p;
  this->last_is_root = is_root;
  this->last_is_sender = is_sender;
  this->last_partner1 = partner1;
  this->last_partner2 = partner2;

  volatile double init_time = MPI_Wtime();
  if (partner1 == -1){
    PMPI_Barrier(cm);
  }
  else {
    double sbuf=0.; double rbuf=0.;
    PMPI_Sendrecv(&sbuf, 1, MPI_DOUBLE, partner1, internal_tag3, &rbuf, 1, MPI_DOUBLE, partner1, internal_tag3, cm, MPI_STATUS_IGNORE);
    if ((partner2 != -1) && (partner1 != partner2)) PMPI_Sendrecv(&sbuf, 1, MPI_DOUBLE, partner2, internal_tag3, &rbuf, 1, MPI_DOUBLE, partner2, internal_tag3, cm, MPI_STATUS_IGNORE);
  }
  this->last_barrier_time = MPI_Wtime() - init_time;

  // Propogate critical paths for all processes in communicator based on what each process has seen up until now (not including this communication)
  propagate(cm, this->last_is_sender, this->last_partner1, this->last_partner2);
  if ((partner2 != -1) && (partner2 != partner2)) propagate(cm, partner2);
  if (partner1 == -1){
    PMPI_Barrier(cm);
  } else {
    double sbuf=0.; double rbuf=0.;
    PMPI_Sendrecv(&sbuf, 1, MPI_DOUBLE, partner1, internal_tag4, &rbuf, 1, MPI_DOUBLE, partner1, internal_tag4, cm, MPI_STATUS_IGNORE);
    if ((partner2 != -1) && (partner1 != partner2)) PMPI_Sendrecv(&sbuf, 1, MPI_DOUBLE, partner2, internal_tag4, &rbuf, 1, MPI_DOUBLE, partner2, internal_tag4, cm, MPI_STATUS_IGNORE);
  }
  // start synchronization timer for communication routine
  this->last_start_time = MPI_Wtime();
}

void synchronous::intermediate(){
  // Deal with synchronization time
  volatile double synchTime = MPI_Wtime();
  this->last_synch_time = synchTime-this->last_start_time;
  // start communication timer for communication routine
  this->last_start_time = MPI_Wtime();
}

void synchronous::stop(){
  volatile double new_time = MPI_Wtime();
  volatile double dt = new_time - this->last_start_time;	// complete communication time
  double datamvt_time = std::max(0.,(dt-this->last_synch_time));
  std::pair<double,double> dcost_simple    = this->cost_func_simple(this->last_nbytes, this->last_p);
  std::pair<double,double> dcost_alphabeta_butterfly = this->cost_func_alphabeta_butterfly(this->last_nbytes, this->last_p);
  std::pair<double,double> dcost_bsp       = this->cost_func_bsp(this->last_nbytes, this->last_p, this->last_is_root);
  std::vector<std::pair<double,double>> dcosts = {dcost_simple,dcost_alphabeta_butterfly,dcost_bsp};

  if (mode == 1){
    *this->my_synch_time   += this->last_synch_time;
    *this->my_datamvt_time += datamvt_time;
    *this->my_comm_time    += dt;
    int save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]){
        *(this->my_msg_count+save) += dcosts[j].first;
        *(this->my_wrd_count+save) += dcosts[j].second;
        save++;
      }
    }
    *this->my_bar_time     += this->last_barrier_time;
    for (size_t i=0; i<breakdown_size; i++){
      *(this->critical_path_synch_time+i)   += this->last_synch_time;
      *(this->critical_path_datamvt_time+i) += datamvt_time;
      *(this->critical_path_comm_time+i)    += dt;
    }
    save=0;
    for (int j=0; j<cost_models.size(); j++){
      for (size_t i=0; i<breakdown_size; i++){
        if (cost_models[j]){
          *(this->critical_path_msg_count+save*breakdown_size+i) += dcosts[j].first;
          *(this->critical_path_wrd_count+save*breakdown_size+i) += dcosts[j].second;
        }
      }
      save++;
    }
  }
  else if (mode == 2){
    // update all communication-related measures for the top symbol in stack
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[0] += this->last_nbytes;
//    symbol_timers[symbol_stack.top()].exclusive_measure.top()[1] += dcost.second;
//    symbol_timers[symbol_stack.top()].exclusive_measure.top()[2] += dcost.first;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[3] += dt;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[4] += this->last_synch_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[5] += datamvt_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[7] += dt;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[0] += this->last_nbytes;
//    *symbol_timers[symbol_stack.top()].cp_excl_measure[1] += dcost.second;
//    *symbol_timers[symbol_stack.top()].cp_excl_measure[2] += dcost.first;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[3] += dt;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[4] += this->last_synch_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[5] += datamvt_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[7] += dt;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[0] += this->last_nbytes;
//    *symbol_timers[symbol_stack.top()].pp_excl_measure[1] += dcost.second;
 //   *symbol_timers[symbol_stack.top()].pp_excl_measure[2] += dcost.first;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[3] += this->last_barrier_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[4] += dt;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[5] += this->last_synch_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[6] += datamvt_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[8] += (dt+this->last_barrier_time);
  }

  int save=0;
  for (int j=0; j<cost_models.size(); j++){
    if (cost_models[j]){
      critical_path_costs[save]                 += dcosts[j].second;		// update critical path estimated communication cost
      critical_path_costs[cost_model_size+save] += dcosts[j].first;		// update critical path estimated synchronization cost
      save++;
    }
  }
  critical_path_costs[num_critical_path_measures-5] += dt;				// update critical path communication time (for what this process has seen thus far)
  critical_path_costs[num_critical_path_measures-4] += this->last_synch_time;	// update critical path synchronization time
  critical_path_costs[num_critical_path_measures-3] += datamvt_time;		// update critical path data mvt time
  critical_path_costs[num_critical_path_measures-1] += dt;				// update critical path runtime

  save=0;
  for (int j=0; j<cost_models.size(); j++){
    if (cost_models[j]){
      volume_costs[save]                 += dcosts[j].second;		// update local estimated communication cost
      volume_costs[cost_model_size+save] += dcosts[j].first;		// update local estimated synchronization cost
      save++;
    }
  }
  volume_costs[num_volume_measures-6] += this->last_barrier_time;		// update local barrier/idle time
  volume_costs[num_volume_measures-5] += dt;				// update local communication time (not volume until after the completion of the program)
  volume_costs[num_volume_measures-4] += this->last_synch_time;		// update local synchronization time
  volume_costs[num_volume_measures-3] += datamvt_time;			// update local data mvt time
  volume_costs[num_volume_measures-1] += this->last_barrier_time;		// update local runtime with idle time
  volume_costs[num_volume_measures-1] += dt;				// update local runtime

  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
  volume_costs[num_volume_measures-5] = volume_costs[num_volume_measures-5] > critical_path_costs[num_critical_path_measures-5]
                                          ? critical_path_costs[num_critical_path_measures-5] : volume_costs[num_volume_measures-5];
  volume_costs[num_volume_measures-4] = volume_costs[num_volume_measures-4] > critical_path_costs[num_critical_path_measures-4]
                                          ? critical_path_costs[num_critical_path_measures-4] : volume_costs[num_volume_measures-4];
  volume_costs[num_volume_measures-3] = volume_costs[num_volume_measures-3] > critical_path_costs[num_critical_path_measures-3]
                                          ? critical_path_costs[num_critical_path_measures-3] : volume_costs[num_volume_measures-3];
  volume_costs[num_volume_measures-2] = volume_costs[num_volume_measures-2] > critical_path_costs[num_critical_path_measures-2]
                                          ? critical_path_costs[num_critical_path_measures-2] : volume_costs[num_volume_measures-2];
  volume_costs[num_volume_measures-1] = volume_costs[num_volume_measures-1] > critical_path_costs[num_critical_path_measures-1]
                                          ? critical_path_costs[num_critical_path_measures-1] : volume_costs[num_volume_measures-1];

  // Prepare to leave interception and re-enter user code
  this->last_start_time = MPI_Wtime();
  computation_timer = this->last_start_time;
  if (mode == 2){
    symbol_timers[symbol_stack.top()].start_timer.top() = this->last_start_time;
  }
}

void blocking::start(volatile double curTime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, bool is_root, bool is_sender, int partner1, int partner2){
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  if (mode == 2){
    this->save_time = curTime - symbol_timers[symbol_stack.top()].start_timer.top();
  }
  this->save_comp_time = curTime - computation_timer;
  this->last_cm = cm;
  this->last_is_root = is_root;
  this->last_is_sender = is_sender;
  this->last_partner1 = partner1;
  this->last_partner2 = partner2;

  int el_size,p;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(cm, &p);
  this->last_nbytes = nbytes;
  this->last_p = p;

  // start timer for communication routine
  this->last_start_time = MPI_Wtime();
}

void blocking::intermediate(){
  // Deal with synchronization time
  volatile double synchTime = MPI_Wtime();
  this->last_synch_time = 0.;//synchTime-this->last_start_time;
  // start communication timer for communication routine
  this->last_start_time = MPI_Wtime();
}

// Used only for p2p communication. All blocking collectives use sychronous protocol
void blocking::stop(){
  volatile double new_time = MPI_Wtime();
  volatile double dt = new_time - this->last_start_time;	// complete communication time
  double datamvt_time = std::max(0.,(dt-this->last_synch_time));
  std::pair<double,double> dcost_simple    = this->cost_func_simple(this->last_nbytes, this->last_p);
  std::pair<double,double> dcost_alphabeta_butterfly = this->cost_func_alphabeta_butterfly(this->last_nbytes, this->last_p);
  std::pair<double,double> dcost_bsp       = this->cost_func_bsp(this->last_nbytes, this->last_p, false);
  std::vector<std::pair<double,double>> dcosts = {dcost_simple,dcost_alphabeta_butterfly,dcost_bsp};

  if (mode == 1){
    *this->my_synch_time   += this->last_synch_time;
    *this->my_datamvt_time += datamvt_time;
    *this->my_comm_time    += dt;
    int save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]){
        *(this->my_msg_count+save) += dcosts[j].first;
        *(this->my_wrd_count+save) += dcosts[j].second;
        save++;
      }
    }
    *this->my_bar_time     += this->last_barrier_time;
    for (size_t i=0; i<breakdown_size; i++){
      *(this->critical_path_synch_time+i)   += this->last_synch_time;
      *(this->critical_path_datamvt_time+i) += datamvt_time;
      *(this->critical_path_comm_time+i)    += dt;
    }
    save=0;
    for (int j=0; j<cost_models.size(); j++){
      for (size_t i=0; i<breakdown_size; i++){
        if (cost_models[j]){
          *(this->critical_path_msg_count+save*breakdown_size+i) += dcosts[j].first;
          *(this->critical_path_wrd_count+save*breakdown_size+i) += dcosts[j].second;
        }
      }
      save++;
    }
  }
  else if (mode == 2){
    // update all communication-related measures for the top symbol in stack
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[0] += this->last_nbytes;
//    symbol_timers[symbol_stack.top()].exclusive_measure.top()[1] += dcost.second;
//    symbol_timers[symbol_stack.top()].exclusive_measure.top()[2] += dcost.first;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[3] += dt;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[4] += this->last_synch_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[5] += datamvt_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[6] += this->save_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[7] += this->save_time+dt;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[0] += this->last_nbytes;
//    *symbol_timers[symbol_stack.top()].cp_excl_measure[1] += dcost.second;
//    *symbol_timers[symbol_stack.top()].cp_excl_measure[2] += dcost.first;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[3] += dt;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[4] += this->last_synch_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[5] += datamvt_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[6] += this->save_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[7] += this->save_time+dt;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[0] += this->last_nbytes;
//    *symbol_timers[symbol_stack.top()].pp_excl_measure[1] += dcost.second;
//    *symbol_timers[symbol_stack.top()].pp_excl_measure[2] += dcost.first;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[4] += dt;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[5] += this->last_synch_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[6] += datamvt_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[7] += this->save_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[8] += this->save_time+dt;
  }

  critical_path_costs[num_critical_path_measures-2] += this->save_comp_time;	// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += this->save_comp_time+dt;	// update critical path runtime
  volume_costs[num_volume_measures-2] += this->save_comp_time;		// update local computation time
  volume_costs[num_volume_measures-1] += this->save_comp_time+dt;		// update local runtime
  for (size_t i=0; i<breakdown_size; i++){
    critical_path_costs[critical_path_costs_size-1-i] += this->save_comp_time;
  }

  int save=0;
  for (int j=0; j<cost_models.size(); j++){
    if (cost_models[j]){
      critical_path_costs[save]                 += dcosts[j].second;		// update critical path estimated communication cost
      critical_path_costs[cost_model_size+save] += dcosts[j].first;		// update critical path estimated synchronization cost
      save++;
    }
  }
  critical_path_costs[num_critical_path_measures-5] += dt;				// update critical path communication time (for what this process has seen thus far)
  critical_path_costs[num_critical_path_measures-4] += this->last_synch_time;	// update critical path synchronization time
  critical_path_costs[num_critical_path_measures-3] += datamvt_time;		// update critical path data mvt time

  save=0;
  for (int j=0; j<cost_models.size(); j++){
    if (cost_models[j]){
      volume_costs[save]                 += dcosts[j].second;		// update local estimated communication cost
      volume_costs[cost_model_size+save] += dcosts[j].first;		// update local estimated synchronization cost
      save++;
    }
  }
  volume_costs[num_volume_measures-6] += this->last_barrier_time;		// update local barrier/idle time
  volume_costs[num_volume_measures-5] += dt;				// update local communication time (not volume until after the completion of the program)
  volume_costs[num_volume_measures-4] += this->last_synch_time;		// update local synchronization time
  volume_costs[num_volume_measures-3] += datamvt_time;			// update local data mvt time

  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
  volume_costs[num_volume_measures-5] = volume_costs[num_volume_measures-5] > critical_path_costs[num_critical_path_measures-5]
                                          ? critical_path_costs[num_critical_path_measures-5] : volume_costs[num_volume_measures-5];
  volume_costs[num_volume_measures-4] = volume_costs[num_volume_measures-4] > critical_path_costs[num_critical_path_measures-4]
                                          ? critical_path_costs[num_critical_path_measures-4] : volume_costs[num_volume_measures-4];
  volume_costs[num_volume_measures-3] = volume_costs[num_volume_measures-3] > critical_path_costs[num_critical_path_measures-3]
                                          ? critical_path_costs[num_critical_path_measures-3] : volume_costs[num_volume_measures-3];
  volume_costs[num_volume_measures-2] = volume_costs[num_volume_measures-2] > critical_path_costs[num_critical_path_measures-2]
                                          ? critical_path_costs[num_critical_path_measures-2] : volume_costs[num_volume_measures-2];
  volume_costs[num_volume_measures-1] = volume_costs[num_volume_measures-1] > critical_path_costs[num_critical_path_measures-1]
                                          ? critical_path_costs[num_critical_path_measures-1] : volume_costs[num_volume_measures-1];

  // Exchange the tracked routine critical path data
  propagate(this->last_cm,this->last_is_sender,this->last_partner1,this->last_partner2);

  // Prepare to leave interception and re-enter user code
  this->last_start_time = MPI_Wtime();
  computation_timer = this->last_start_time;
  if (mode == 2){
    symbol_timers[symbol_stack.top()].start_timer.top() = this->last_start_time;
  }
}

// Called by both nonblocking p2p and nonblocking collectives
void nonblocking::start(volatile double curTime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, MPI_Request* request, bool is_sender, int partner){
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  this->save_comp_time = curTime - computation_timer;
  critical_path_costs[num_critical_path_measures-2] += this->save_comp_time;		// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += this->save_comp_time;		// update critical path runtime
  volume_costs[num_volume_measures-2]        += this->save_comp_time;		// update local computation time
  volume_costs[num_volume_measures-1]        += this->save_comp_time;		// update local runtime
  for (size_t i=0; i<breakdown_size; i++){
    critical_path_costs[critical_path_costs_size-1-i] += this->save_comp_time;
  }
  if (mode == 2){
    assert(symbol_stack.size()>0);
    assert(symbol_timers[symbol_stack.top()].start_timer.size()>0);
    this->save_time = curTime - symbol_timers[symbol_stack.top()].start_timer.top();
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[num_critical_path_measures-1] += this->save_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[num_critical_path_measures-2] += this->save_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-1] += this->save_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-2] += this->save_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-1] += this->save_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-2] += this->save_time;
  }

  int el_size,p;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(cm, &p);

  // Nonblocking communication to propogate the critical path from sender to receiver. Avoids tricky deadlock in intercepting MPI_Waitall
  // Unlike blocking protocol, Receiver does not need sender's critical path information to include the contribution from this current routine
  MPI_Request internal_request;
  double* data = (double*)malloc(sizeof(double)*critical_path_costs.size());
  // Save local data instead of immediately adding it to critical path, because this communication is not technically completed yet,
  //   and I do not want to corrupt critical path propogation in future communication that may occur before this nonblocking communication completes.
  if (partner == -1){
    for (int i=0; i<critical_path_costs.size(); i++){
      data[i] = critical_path_costs[i];
    }
    MPI_Op op; MPI_Op_create((MPI_User_function*) propagate_critical_path_op,1,&op);
    PMPI_Iallreduce(MPI_IN_PLACE,&data[0],critical_path_costs.size(),MPI_DOUBLE,op,cm,&internal_request);
    //MPI_Op_free(&op);
  } else{
    if (is_sender){
      for (int i=0; i<critical_path_costs.size(); i++){
        data[i] = critical_path_costs[i];
      }
      PMPI_Isend(&data[0],critical_path_costs.size(),MPI_DOUBLE,partner,internal_tag,cm,&internal_request);
    }
    else{
      PMPI_Irecv(&data[0],critical_path_costs.size(),MPI_DOUBLE,partner,internal_tag,cm,&internal_request);
    }
  }
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  internal_comm_info[*request] = std::make_pair(internal_request,is_sender);
  internal_comm_comm[*request] = std::make_pair(cm,partner);
  internal_comm_message[*request] = data;
  internal_comm_data[*request] = std::make_pair((double)nbytes,(double)p);
  internal_comm_track[*request] = this;
}

void nonblocking::stop(MPI_Request* request, double comp_time, double comm_time){
  volatile double new_time = MPI_Wtime();
  auto comm_info_it = internal_comm_info.find(*request);
  auto comm_comm_it = internal_comm_comm.find(*request);
  auto comm_message_it = internal_comm_message.find(*request);
  auto comm_data_it = internal_comm_data.find(*request);
  auto comm_track_it = internal_comm_track.find(*request);
  assert(comm_info_it != internal_comm_info.end());
  assert(comm_comm_it != internal_comm_comm.end());
  assert(comm_message_it != internal_comm_message.end());
  assert(comm_data_it != internal_comm_data.end());
  assert(comm_track_it != internal_comm_track.end());

  // Before accumulating the cost of this communication into our critical path/volume measures, we
  //   must first finish the internal communication, which doesn't take into account the cost of this communication
  // The computation and communication time of the sender between its MPI_Isend and MPI_Wait cannot be tracked. The receiver
  //   will just use its own. This is technically ok I think, because the receiver isn't waiting on the receiver in any capacity.

  MPI_Request internal_request = comm_info_it->second.first;
  bool is_sender = comm_info_it->second.second;
  MPI_Comm cm = comm_comm_it->second.first;
  int partner = comm_comm_it->second.second;
  double nbytes = comm_data_it->second.first;
  double p = comm_data_it->second.second;
  double* data = comm_message_it->second;

  // TODO: For mode==2, I think data will need to be double_int* instead of double*
  propagate(data,internal_request,cm,partner,is_sender);

  // Both sender and receiver will now update its critical path with the data from the communication
  std::pair<double,double> dcost_simple  = this->cost_func_simple(nbytes,p);
  std::pair<double,double> dcost_alphabeta_butterfly = this->cost_func_alphabeta_butterfly(nbytes,p);
  std::pair<double,double> dcost_bsp     = this->cost_func_bsp(nbytes,p,false);
  if ((this->tag<20) && (wait_id)) dcost_bsp.first=1.;	// this is usually zero, but we force it to be 1 in special circumstances (for nonblocking p2p with wait_id one)
  std::vector<std::pair<double,double>> dcosts = {dcost_simple,dcost_alphabeta_butterfly,dcost_bsp};

  if (mode == 1){
    *this->my_synch_time   += 0;			// Nonblocking routines will have no synchronization time component
    *this->my_datamvt_time += comm_time;
    *this->my_comm_time    += comm_time;
    int save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]){
        *(this->my_msg_count+save) += dcosts[j].first;
        *(this->my_wrd_count+save) += dcosts[j].second;
        save++;
      }
    }
    save=0;
    for (size_t i=0; i<breakdown_size; i++){
      *(this->critical_path_synch_time+i)   += this->last_synch_time;
      *(this->critical_path_datamvt_time+i) += comm_time;
      *(this->critical_path_comm_time+i)    += comm_time;
    }
    for (int j=0; j<cost_models.size(); j++){
      for (size_t i=0; i<breakdown_size; i++){
        if (cost_models[j]){
          *(this->critical_path_msg_count+save*breakdown_size+i) += dcosts[j].first;
          *(this->critical_path_wrd_count+save*breakdown_size+i) += dcosts[j].second;
        }
      }
      save++;
    }
  }
  else if (mode == 2){
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[0] += this->last_nbytes;
//    symbol_timers[symbol_stack.top()].exclusive_measure.top()[1] += dcost.second;
//    symbol_timers[symbol_stack.top()].exclusive_measure.top()[2] += dcost.first;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[3] += comm_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[4] += 0;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[5] += comm_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[6] += comp_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[7] += (comp_time+comm_time);
    *symbol_timers[symbol_stack.top()].cp_excl_measure[0] += this->last_nbytes;
//    *symbol_timers[symbol_stack.top()].cp_excl_measure[1] += dcost.second;
//    *symbol_timers[symbol_stack.top()].cp_excl_measure[2] += dcost.first;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[3] += comm_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[4] += 0;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[5] += comm_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[6] += comp_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[7] += (comp_time+comm_time);
    *symbol_timers[symbol_stack.top()].pp_excl_measure[0] += this->last_nbytes;
//    *symbol_timers[symbol_stack.top()].pp_excl_measure[1] += dcost.second;
//    *symbol_timers[symbol_stack.top()].pp_excl_measure[2] += dcost.first;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[4] += comm_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[5] += 0;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[6] += comm_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[7] += comp_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[8] += (comp_time+comm_time);
  }

  int save=0;
  for (int j=0; j<cost_models.size(); j++){
    if (cost_models[j]){
      critical_path_costs[save]                  += dcosts[j].second;		// update critical path estimated communication cost
      critical_path_costs[cost_model_size+save] += dcosts[j].first;		// update critical path estimated synchronization cost
      save++;
    }
  }
  critical_path_costs[num_critical_path_measures-5] += comm_time;				// update critical path communication time (for what this process has seen thus far)
  critical_path_costs[num_critical_path_measures-4] += 0.;				// update critical path synchronization time
  critical_path_costs[num_critical_path_measures-3] += comm_time;				// update critical path runtime
  critical_path_costs[num_critical_path_measures-2] += comp_time;				// update critical path runtime
  critical_path_costs[num_critical_path_measures-1] += comp_time+comm_time;				// update critical path runtime
  for (size_t i=0; i<breakdown_size; i++){
    critical_path_costs[critical_path_costs_size-1-i] += comp_time;
  }

  save=0;
  for (int j=0; j<cost_models.size(); j++){
    if (cost_models[j]){
      volume_costs[save]                  += dcosts[j].second;		// update local estimated communication cost
      volume_costs[cost_model_size+save] += dcosts[j].first;		// update local estimated synchronization cost
      save++;
    }
  }
  volume_costs[num_volume_measures-5] += comm_time;				// update local communication time (not volume until after the completion of the program)
  volume_costs[num_volume_measures-4] += 0.;		// update local synchronization time
  volume_costs[num_volume_measures-3] += comm_time;			// update local data mvt time
  volume_costs[num_volume_measures-2] += comp_time;				// update local runtime
  volume_costs[num_volume_measures-1] += comp_time+comm_time;				// update local runtime

  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
  volume_costs[num_volume_measures-5] = volume_costs[num_volume_measures-5] > critical_path_costs[num_critical_path_measures-5]
                                          ? critical_path_costs[num_critical_path_measures-5] : volume_costs[num_volume_measures-5];
  volume_costs[num_volume_measures-4] = volume_costs[num_volume_measures-4] > critical_path_costs[num_critical_path_measures-4]
                                          ? critical_path_costs[num_critical_path_measures-4] : volume_costs[num_volume_measures-4];
  volume_costs[num_volume_measures-3] = volume_costs[num_volume_measures-3] > critical_path_costs[num_critical_path_measures-3]
                                          ? critical_path_costs[num_critical_path_measures-3] : volume_costs[num_volume_measures-3];
  volume_costs[num_volume_measures-2] = volume_costs[num_volume_measures-2] > critical_path_costs[num_critical_path_measures-2]
                                          ? critical_path_costs[num_critical_path_measures-2] : volume_costs[num_volume_measures-2];
  volume_costs[num_volume_measures-1] = volume_costs[num_volume_measures-1] > critical_path_costs[num_critical_path_measures-1]
                                          ? critical_path_costs[num_critical_path_measures-1] : volume_costs[num_volume_measures-1];

  internal_comm_info.erase(*request);
  internal_comm_comm.erase(*request);
  free(comm_message_it->second);
  internal_comm_message.erase(*request);
  internal_comm_data.erase(*request);
  internal_comm_track.erase(*request);

  this->last_start_time = MPI_Wtime();
  computation_timer = this->last_start_time;
  if (mode == 2){
    symbol_timers[symbol_stack.top()].start_timer.top() = this->last_start_time;
  }
}

void complete_propagation(MPI_Comm cm, bool is_sender, int partner1, int partner2){
  if (mode <= 1){// Note that the only way this routine would be called with mode==0 is after critter::stop
    // First exchange the tracked routine critical path data
    if (partner1 == -1){
      MPI_Op op; MPI_Op_create((MPI_User_function*) propagate_critical_path_op,1,&op);
      PMPI_Allreduce(MPI_IN_PLACE, &critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, op, cm);
      MPI_Op_free(&op);
    }
    else if (partner2 == -1){
      // Sender sends critical path data to receiver. Receiver updates its critical path information.
      if (is_sender){
        PMPI_Send(&critical_path_costs[0],critical_path_costs.size(),MPI_DOUBLE,partner1,internal_tag,cm);
        // Sender needs not wait for handshake with receiver, can continue back into user code
      }
      else{
        PMPI_Recv(&new_cs[0],critical_path_costs.size(),MPI_DOUBLE,partner1,internal_tag,cm,MPI_STATUS_IGNORE);
        update_critical_path(&new_cs[0],&critical_path_costs[0],critical_path_costs_size);
      }
    }
    else{
      PMPI_Sendrecv(&critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, partner1, internal_tag1, &new_cs[0], critical_path_costs.size(),
        MPI_DOUBLE, partner2, internal_tag1, cm, MPI_STATUS_IGNORE);
      update_critical_path(&new_cs[0],&critical_path_costs[0],critical_path_costs_size);
    }
  }
  else if (mode == 2){
/*
    // Note: challenge of updating symbol_stack. For local symbols that are not along the updated critical path, we need to just make their contributions zero, but still keep them in the maps.
    // Next, exchange the critical path metric, together with tracking the rank of the process that determines each critical path
    //TODO: Note that this is missing non-runtime-critical-path breakdown, as is performed above for mode==1 with the MPI_Op or the function 'update_critical_path'

    int rank; MPI_Comm_rank(cm,&rank); int true_rank; MPI_Comm_rank(MPI_COMM_WORLD,&true_rank);
    for (int i=0; i<num_critical_path_measures; i++){
      timer_info_sender[i].first = critical_path_costs[i];
      timer_info_sender[i].second = rank;
    }
    if (partner == -1){
      PMPI_Allreduce(&timer_info_sender[0].first, &timer_info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, MPI_MAXLOC, cm);
    }
    else {
      PMPI_Sendrecv(&timer_info_sender[0].first, num_critical_path_measures, MPI_DOUBLE_INT, partner, internal_tag, &timer_info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, partner, internal_tag, cm, MPI_STATUS_IGNORE);
      for (int i=0; i<num_critical_path_measures; i++){
        if (timer_info_sender[i].first>timer_info_receiver[i].first){timer_info_receiver[i].second = rank;}
        else if (timer_info_sender[i].first==timer_info_receiver[i].first){
          if (timer_info_sender[i].second < timer_info_receiver[i].second){
            timer_info_receiver[i].second = rank;
          } else{
            timer_info_receiver[i].second = partner;
          }
        }
        timer_info_receiver[i].first = std::max(timer_info_sender[i].first, timer_info_receiver[i].first);
      }
    }
    int critical_path_runtime_root_rank = timer_info_receiver[num_critical_path_measures-1].second;
    for (int i=0; i<num_critical_path_measures; i++){
      critical_path_costs[i] = timer_info_receiver[i].first;
    }

    // We consider only critical path runtime
    std::array<int,2> ftimer_size = {0,0};
    if (rank==critical_path_runtime_root_rank){
      ftimer_size[0] = symbol_timers.size();
      ftimer_size[1] = symbol_stack.size() > 0 ? symbol_timers[symbol_stack.top()].exclusive_contributions.top().size() : 0;
    }
    if (partner == -1){
      PMPI_Allreduce(MPI_IN_PLACE,&ftimer_size[0],2,MPI_INT,MPI_SUM,cm);
    }
    else{
      if (rank != partner){
        if (rank==critical_path_runtime_root_rank){
          PMPI_Send(&ftimer_size[0],2,MPI_INT,partner,internal_tag,cm);
        }
        else{
          PMPI_Recv(&ftimer_size[0],2,MPI_INT,partner,internal_tag,cm,MPI_STATUS_IGNORE);
        }
      }
    }

    std::vector<int> symbol_sizes(ftimer_size[0]+ftimer_size[1],0);
    std::vector<double> exclusive_contributions(ftimer_size[1]*num_critical_path_measures,0);
    if (rank==critical_path_runtime_root_rank){
      int symbol_offset = 0;
      for (auto i=0; i<symbol_timers.size(); i++){
        symbol_sizes[i] = symbol_order[i].size();
        for (auto j=0; j<symbol_sizes[i]; j++){
          symbol_pad[symbol_offset+j] = symbol_order[i][j];
        }
        symbol_offset += symbol_sizes[i];
      }
      if (symbol_stack.size()>0){
        int index = 0;
        for (auto& it : symbol_timers[symbol_stack.top()].exclusive_contributions.top()){
          symbol_sizes[index+symbol_timers.size()] = it.first.size();
          for (auto j=0; j<it.first.size(); j++){
            symbol_pad[symbol_offset+j] = it.first[j];
          }
          for (auto j=0; j<num_critical_path_measures; j++){
            exclusive_contributions[index*num_critical_path_measures+j] = it.second[j];
          }
          symbol_offset += it.first.size(); index++;
        }
      }
    }
    if (partner == -1){
      PMPI_Allreduce(MPI_IN_PLACE,&symbol_sizes[0],ftimer_size[0]+ftimer_size[1],MPI_INT,MPI_SUM,cm);
    }
    else{
      if (rank != partner){
        if (rank==critical_path_runtime_root_rank){
          PMPI_Send(&symbol_sizes[0],ftimer_size[0]+ftimer_size[1],MPI_INT,partner,internal_tag,cm);
        }
        else{
          PMPI_Recv(&symbol_sizes[0],ftimer_size[0]+ftimer_size[1],MPI_INT,partner,internal_tag,cm,MPI_STATUS_IGNORE);
        }
      }
    }

    int num_chars = 0;
    for (auto i=0; i<ftimer_size[0]+ftimer_size[1]; i++){
      num_chars += symbol_sizes[i];
    }
    if (rank == critical_path_runtime_root_rank){
      if (partner == -1){
        PMPI_Bcast(&symbol_timer_pad_local_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size[0],MPI_DOUBLE,rank,cm);
        PMPI_Bcast(&exclusive_contributions[0],ftimer_size[1]*num_critical_path_measures,MPI_DOUBLE,rank,cm);
        PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,rank,cm);
      }
      else{
        if (rank != partner){
          PMPI_Send(&symbol_timer_pad_local_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size[0],MPI_DOUBLE,partner,internal_tag,cm);
          PMPI_Send(&exclusive_contributions[0],ftimer_size[1]*num_critical_path_measures,MPI_DOUBLE,partner,internal_tag,cm);
          PMPI_Send(&symbol_pad[0],num_chars,MPI_CHAR,partner,internal_tag,cm);
        }
      }
    }
    else{
      if (partner == -1){
        PMPI_Bcast(&symbol_timer_pad_global_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size[0],MPI_DOUBLE,critical_path_runtime_root_rank,cm);
        PMPI_Bcast(&exclusive_contributions[0], ftimer_size[1]*num_critical_path_measures,MPI_DOUBLE,critical_path_runtime_root_rank,cm);
        PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,critical_path_runtime_root_rank,cm);
      }
      else{
        if (rank != partner){
          PMPI_Recv(&symbol_timer_pad_global_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size[0],MPI_DOUBLE,partner,internal_tag,cm,MPI_STATUS_IGNORE);
          PMPI_Recv(&exclusive_contributions[0],ftimer_size[1]*num_critical_path_measures,MPI_DOUBLE,partner,internal_tag,cm,MPI_STATUS_IGNORE);
          PMPI_Recv(&symbol_pad[0],num_chars,MPI_CHAR,partner,internal_tag,cm,MPI_STATUS_IGNORE);
        }
      }
      if (rank != partner){
        int symbol_offset = 0;
        for (int i=0; i<ftimer_size[0]; i++){
          auto reconstructed_symbol = std::string(symbol_pad.begin()+symbol_offset,symbol_pad.begin()+symbol_offset+symbol_sizes[i]);
          if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
            symbol_timers[reconstructed_symbol] = ftimer(reconstructed_symbol,true);
            symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
          }
          *symbol_timers[reconstructed_symbol].cp_numcalls = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i];
          for (int j=0; j<num_critical_path_measures; j++){
            *symbol_timers[reconstructed_symbol].cp_incl_measure[j] = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i+2*j+1];
            *symbol_timers[reconstructed_symbol].cp_excl_measure[j] = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i+2*(j+1)];
          }
          symbol_timers[reconstructed_symbol].has_been_processed = true;
          symbol_offset += symbol_sizes[i];
        }
        if (symbol_stack.size()>0) { symbol_timers[symbol_stack.top()].exclusive_contributions.top().clear(); }
        for (int i=0; i<ftimer_size[1]; i++){
          auto reconstructed_symbol = std::string(symbol_pad.begin()+symbol_offset,symbol_pad.begin()+symbol_offset+symbol_sizes[ftimer_size[0]+i]);
          symbol_timers[symbol_stack.top()].exclusive_contributions.top()[reconstructed_symbol].fill(0.0);
          for (int j=0; j<num_critical_path_measures; j++){
            symbol_timers[symbol_stack.top()].exclusive_contributions.top()[reconstructed_symbol][j] = exclusive_contributions[i*num_critical_path_measures+j];
          }
          symbol_offset += symbol_sizes[ftimer_size[0]+i];
        }
        // Now cycle through and find the symbols that were not processed and set their accumulated measures to 0
        for (auto& it : symbol_timers){
          if (it.second.has_been_processed){ it.second.has_been_processed = false; }
          else{
            *it.second.cp_numcalls = 0;
            for (int j=0; j<num_critical_path_measures; j++){
              *it.second.cp_incl_measure[j] = 0;
              *it.second.cp_excl_measure[j] = 0;
            }
          }
        }
      }
    }
*/
  }
}

void synchronous::propagate(MPI_Comm cm, bool is_sender, int partner1, int partner2){
  complete_propagation(cm,is_sender,partner1,partner2);
}

void blocking::propagate(MPI_Comm cm, bool is_sender, int partner1, int partner2){
  if (mode == 1){
    if (partner2==-1){
      // Sender sends critical path data to receiver. Receiver updates its critical path information.
      if (is_sender){
        PMPI_Send(&critical_path_costs[0],critical_path_costs.size(),MPI_DOUBLE,partner1,internal_tag,cm);
        // Sender needs not wait for handshake with receiver, can continue back into user code
      }
      else{
        PMPI_Recv(&new_cs[0],critical_path_costs.size(),MPI_DOUBLE,partner1,internal_tag,cm,MPI_STATUS_IGNORE);
        update_critical_path(&new_cs[0],&critical_path_costs[0],critical_path_costs_size);
      }
    }
    else{
      PMPI_Sendrecv(&critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, partner1, internal_tag1, &new_cs[0], critical_path_costs.size(),
        MPI_DOUBLE, partner2, internal_tag1, cm, MPI_STATUS_IGNORE);
      update_critical_path(&new_cs[0],&critical_path_costs[0],critical_path_costs_size);
    }
  }
  else if (mode == 2){
/*
    // Note: challenge of updating symbol_stack. For local symbols that are not along the updated critical path, we need to just make their contributions zero, but still keep them in the maps.
    // Next, exchange the critical path metric, together with tracking the rank of the process that determines each critical path
    //TODO: Note that this is missing non-runtime-critical-path breakdown, as is performed above for mode==1 with the MPI_Op or the function 'update_critical_path'
    if (is_sender){
      PMPI_Send(&critical_path_costs[0],critical_path_costs.size(),MPI_DOUBLE,partner,internal_tag,cm);
      // We consider only critical path runtime
      std::array<int,2> ftimer_size = {0,0};
      ftimer_size[0] = symbol_timers.size();
      ftimer_size[1] = symbol_stack.size() > 0 ? symbol_timers[symbol_stack.top()].exclusive_contributions.top().size() : 0;
      PMPI_Send(&ftimer_size[0],2,MPI_INT,partner,internal_tag,cm);
      std::vector<int> symbol_sizes(ftimer_size[0]+ftimer_size[1],0);
      std::vector<double> exclusive_contributions(ftimer_size[1]*num_critical_path_measures,0);
      int symbol_offset = 0;
      for (auto i=0; i<symbol_timers.size(); i++){
        symbol_sizes[i] = symbol_order[i].size();
        for (auto j=0; j<symbol_sizes[i]; j++){
          symbol_pad[symbol_offset+j] = symbol_order[i][j];
        }
       symbol_offset += symbol_sizes[i];
      }
      if (symbol_stack.size()>0){
        int index = 0;
        for (auto& it : symbol_timers[symbol_stack.top()].exclusive_contributions.top()){
          symbol_sizes[index+symbol_timers.size()] = it.first.size();
          for (auto j=0; j<it.first.size(); j++){
            symbol_pad[symbol_offset+j] = it.first[j];
          }
          for (auto j=0; j<num_critical_path_measures; j++){
            exclusive_contributions[index*num_critical_path_measures+j] = it.second[j];
          }
          symbol_offset += it.first.size(); index++;
        }
      }
      PMPI_Send(&symbol_sizes[0],ftimer_size[0]+ftimer_size[1],MPI_INT,partner,internal_tag,cm);
      int num_chars = 0;
      for (auto i=0; i<ftimer_size[0]+ftimer_size[1]; i++){
        num_chars += symbol_sizes[i];
      }
      PMPI_Send(&symbol_timer_pad_local_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size[0],MPI_DOUBLE,partner,internal_tag,cm);
      PMPI_Send(&exclusive_contributions[0],ftimer_size[1]*num_critical_path_measures,MPI_DOUBLE,partner,internal_tag,cm);
      PMPI_Send(&symbol_pad[0],num_chars,MPI_CHAR,partner,internal_tag,cm);
    }
    else{
      PMPI_Recv(&new_cs[0],critical_path_costs.size(),MPI_DOUBLE,partner,internal_tag,cm,MPI_STATUS_IGNORE);
      bool update_path = critical_path_costs[7] < new_cs[7];	// again, for now, we consider only the runtime critical path, so if receiving process has larger runtime critical path, it does not need to update its critical path
      for (int i=0; i<num_critical_path_measures; i++){
        critical_path_costs[i] = std::max(new_cs[i],critical_path_costs[i]);
      }
      // We consider only critical path runtime
      std::array<int,2> ftimer_size = {0,0};
      PMPI_Recv(&ftimer_size[0],2,MPI_INT,partner,internal_tag,cm,MPI_STATUS_IGNORE);
      std::vector<int> symbol_sizes(ftimer_size[0]+ftimer_size[1],0);
      std::vector<double> exclusive_contributions(ftimer_size[1]*num_critical_path_measures,0);
      PMPI_Recv(&symbol_sizes[0],ftimer_size[0]+ftimer_size[1],MPI_INT,partner,internal_tag,cm,MPI_STATUS_IGNORE);
      int num_chars = 0;
      for (auto i=0; i<ftimer_size[0]+ftimer_size[1]; i++){
        num_chars += symbol_sizes[i];
      }
      PMPI_Recv(&symbol_timer_pad_global_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size[0],MPI_DOUBLE,partner,internal_tag,cm,MPI_STATUS_IGNORE);
      PMPI_Recv(&exclusive_contributions[0],ftimer_size[1]*num_critical_path_measures,MPI_DOUBLE,partner,internal_tag,cm,MPI_STATUS_IGNORE);
      PMPI_Recv(&symbol_pad[0],num_chars,MPI_CHAR,partner,internal_tag,cm,MPI_STATUS_IGNORE);
      if (update_path){
        int symbol_offset = 0;
        for (int i=0; i<ftimer_size[0]; i++){
          auto reconstructed_symbol = std::string(symbol_pad.begin()+symbol_offset,symbol_pad.begin()+symbol_offset+symbol_sizes[i]);

          if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
            symbol_timers[reconstructed_symbol] = ftimer(reconstructed_symbol,true);
            symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
          }
          *symbol_timers[reconstructed_symbol].cp_numcalls = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i];
          for (int j=0; j<num_critical_path_measures; j++){
            *symbol_timers[reconstructed_symbol].cp_incl_measure[j] = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i+2*j+1];
            *symbol_timers[reconstructed_symbol].cp_excl_measure[j] = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i+2*(j+1)];
          }
          symbol_timers[reconstructed_symbol].has_been_processed = true;
          symbol_offset += symbol_sizes[i];
        }
        if (symbol_stack.size()>0) { symbol_timers[symbol_stack.top()].exclusive_contributions.top().clear(); }
        for (int i=0; i<ftimer_size[1]; i++){
          auto reconstructed_symbol = std::string(symbol_pad.begin()+symbol_offset,symbol_pad.begin()+symbol_offset+symbol_sizes[ftimer_size[0]+i]);
          symbol_timers[symbol_stack.top()].exclusive_contributions.top()[reconstructed_symbol].fill(0.0);
          for (int j=0; j<num_critical_path_measures; j++){
            symbol_timers[symbol_stack.top()].exclusive_contributions.top()[reconstructed_symbol][j] = exclusive_contributions[i*num_critical_path_measures+j];
          }
          symbol_offset += symbol_sizes[ftimer_size[0]+i];
        }
        // Now cycle through and find the symbols that were not processed and set their accumulated measures to 0
        for (auto& it : symbol_timers){
          if (it.second.has_been_processed){ it.second.has_been_processed = false; }
          else{
            *it.second.cp_numcalls = 0;
            for (int j=0; j<num_critical_path_measures; j++){
              *it.second.cp_incl_measure[j] = 0;
              *it.second.cp_excl_measure[j] = 0;
            }
          }
        }
      }
    }
*/
  }
}

void nonblocking::propagate(double* data, MPI_Request internal_request, MPI_Comm cm, int partner, bool is_sender){
  // First exchange the tracked routine critical path data
  MPI_Status st;
  PMPI_Wait(&internal_request,&st);
  if (mode == 1){
    if (!is_sender){
      update_critical_path(data,&critical_path_costs[0],critical_path_costs_size);
    }
  }
  else if (mode == 2){
    // Note: challenge of updating symbol_stack. For local symbols that are not along the updated critical path, we need to just make their contributions zero, but still keep them in the maps.
    // Next, exchange the critical path metric, together with tracking the rank of the process that determines each critical path
    //TODO: Note that this is missing non-runtime-critical-path breakdown, as is performed above for mode==1 with the MPI_Op or the function 'update_critical_path'
    int rank; MPI_Comm_rank(cm,&rank); int true_rank; MPI_Comm_rank(MPI_COMM_WORLD,&true_rank);
    for (int i=0; i<num_critical_path_measures; i++){
      timer_info_sender[i].first = critical_path_costs[i];
      timer_info_sender[i].second = rank;
    }
    if (partner == -1){
      PMPI_Allreduce(&timer_info_sender[0].first, &timer_info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, MPI_MAXLOC, cm);
    }
    else {
      PMPI_Sendrecv(&timer_info_sender[0].first, num_critical_path_measures, MPI_DOUBLE_INT, partner, internal_tag2, &timer_info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, partner, internal_tag2, cm, MPI_STATUS_IGNORE);
      for (int i=0; i<num_critical_path_measures; i++){
        if (timer_info_sender[i].first>timer_info_receiver[i].first){timer_info_receiver[i].second = rank;}
        else if (timer_info_sender[i].first==timer_info_receiver[i].first){
          if (timer_info_sender[i].second < timer_info_receiver[i].second){
            timer_info_receiver[i].second = rank;
          } else{
            timer_info_receiver[i].second = partner;
          }
        }
        timer_info_receiver[i].first = std::max(timer_info_sender[i].first, timer_info_receiver[i].first);
      }
    }

    int critical_path_runtime_root_rank = timer_info_receiver[num_critical_path_measures-1].second;
    for (int i=0; i<num_critical_path_measures; i++){
      critical_path_costs[i] = timer_info_receiver[i].first;
    }
    // We consider only critical path runtime

    std::array<int,2> ftimer_size = {0,0};
    if (rank==critical_path_runtime_root_rank){
      ftimer_size[0] = symbol_timers.size();
      ftimer_size[1] = symbol_stack.size() > 0 ? symbol_timers[symbol_stack.top()].exclusive_contributions.top().size() : 0;
    }
    if (partner == -1){
      PMPI_Allreduce(MPI_IN_PLACE,&ftimer_size[0],2,MPI_INT,MPI_SUM,cm);
    }
    else{
      if (rank != partner){
        if (rank==critical_path_runtime_root_rank){
          PMPI_Send(&ftimer_size[0],2,MPI_INT,partner,internal_tag,cm);
        }
        else{
          PMPI_Recv(&ftimer_size[0],2,MPI_INT,partner,internal_tag,cm,MPI_STATUS_IGNORE);
        }
      }
    }

    std::vector<int> symbol_sizes(ftimer_size[0]+ftimer_size[1],0);
    std::vector<double> exclusive_contributions(ftimer_size[1]*num_critical_path_measures,0);
    if (rank==critical_path_runtime_root_rank){
      int symbol_offset = 0;
      for (auto i=0; i<symbol_timers.size(); i++){
        symbol_sizes[i] = symbol_order[i].size();
        for (auto j=0; j<symbol_sizes[i]; j++){
          symbol_pad[symbol_offset+j] = symbol_order[i][j];
        }
        symbol_offset += symbol_sizes[i];
      }
      if (symbol_stack.size()>0){
        int index = 0;
        for (auto& it : symbol_timers[symbol_stack.top()].exclusive_contributions.top()){
          symbol_sizes[index+symbol_timers.size()] = it.first.size();
          for (auto j=0; j<it.first.size(); j++){
            symbol_pad[symbol_offset+j] = it.first[j];
          }
          for (auto j=0; j<num_critical_path_measures; j++){
            exclusive_contributions[index*num_critical_path_measures+j] = it.second[j];
          }
          symbol_offset += it.first.size(); index++;
        }
      }
    }
    if (partner == -1){
      PMPI_Allreduce(MPI_IN_PLACE,&symbol_sizes[0],ftimer_size[0]+ftimer_size[1],MPI_INT,MPI_SUM,cm);
    }
    else{
      if (rank != partner){
        if (rank==critical_path_runtime_root_rank){
          PMPI_Send(&symbol_sizes[0],ftimer_size[0]+ftimer_size[1],MPI_INT,partner,internal_tag,cm);
        }
        else{
          PMPI_Recv(&symbol_sizes[0],ftimer_size[0]+ftimer_size[1],MPI_INT,partner,internal_tag,cm,MPI_STATUS_IGNORE);
        }
      }
    }

    int num_chars = 0;
    for (auto i=0; i<ftimer_size[0]+ftimer_size[1]; i++){
      num_chars += symbol_sizes[i];
    }
    if (rank == critical_path_runtime_root_rank){
      if (partner == -1){
        PMPI_Bcast(&symbol_timer_pad_local_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size[0],MPI_DOUBLE,rank,cm);
        PMPI_Bcast(&exclusive_contributions[0],ftimer_size[1]*num_critical_path_measures,MPI_DOUBLE,rank,cm);
        PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,rank,cm);
      }
      else{
        if (rank != partner){
          PMPI_Send(&symbol_timer_pad_local_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size[0],MPI_DOUBLE,partner,internal_tag,cm);
          PMPI_Send(&exclusive_contributions[0],ftimer_size[1]*num_critical_path_measures,MPI_DOUBLE,partner,internal_tag,cm);
          PMPI_Send(&symbol_pad[0],num_chars,MPI_CHAR,partner,internal_tag,cm);
        }
      }
    }
    else{
      if (partner == -1){
        PMPI_Bcast(&symbol_timer_pad_global_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size[0],MPI_DOUBLE,critical_path_runtime_root_rank,cm);
        PMPI_Bcast(&exclusive_contributions[0], ftimer_size[1]*num_critical_path_measures,MPI_DOUBLE,critical_path_runtime_root_rank,cm);
        PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,critical_path_runtime_root_rank,cm);
      }
      else{
        if (rank != partner){
          PMPI_Recv(&symbol_timer_pad_global_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size[0],MPI_DOUBLE,partner,internal_tag,cm,MPI_STATUS_IGNORE);
          PMPI_Recv(&exclusive_contributions[0],ftimer_size[1]*num_critical_path_measures,MPI_DOUBLE,partner,internal_tag,cm,MPI_STATUS_IGNORE);
          PMPI_Recv(&symbol_pad[0],num_chars,MPI_CHAR,partner,internal_tag,cm,MPI_STATUS_IGNORE);
        }
      }
      if (rank != partner){
        int symbol_offset = 0;
        for (int i=0; i<ftimer_size[0]; i++){
          auto reconstructed_symbol = std::string(symbol_pad.begin()+symbol_offset,symbol_pad.begin()+symbol_offset+symbol_sizes[i]);

          if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
            symbol_timers[reconstructed_symbol] = ftimer(reconstructed_symbol,true);
            symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
          }
          *symbol_timers[reconstructed_symbol].cp_numcalls = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i];
          for (int j=0; j<num_critical_path_measures; j++){
            *symbol_timers[reconstructed_symbol].cp_incl_measure[j] = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i+2*j+1];
            *symbol_timers[reconstructed_symbol].cp_excl_measure[j] = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i+2*(j+1)];
          }
          symbol_timers[reconstructed_symbol].has_been_processed = true;
          symbol_offset += symbol_sizes[i];
        }
        if (symbol_stack.size()>0) { symbol_timers[symbol_stack.top()].exclusive_contributions.top().clear(); }
        for (int i=0; i<ftimer_size[1]; i++){
          auto reconstructed_symbol = std::string(symbol_pad.begin()+symbol_offset,symbol_pad.begin()+symbol_offset+symbol_sizes[ftimer_size[0]+i]);
          symbol_timers[symbol_stack.top()].exclusive_contributions.top()[reconstructed_symbol].fill(0.0);
          for (int j=0; j<num_critical_path_measures; j++){
            symbol_timers[symbol_stack.top()].exclusive_contributions.top()[reconstructed_symbol][j] = exclusive_contributions[i*num_critical_path_measures+j];
          }
          symbol_offset += symbol_sizes[ftimer_size[0]+i];
        }
        // Now cycle through and find the symbols that were not processed and set their accumulated measures to 0
        for (auto& it : symbol_timers){
          if (it.second.has_been_processed){ it.second.has_been_processed = false; }
          else{
            *it.second.cp_numcalls = 0;
            for (int j=0; j<num_critical_path_measures; j++){
              *it.second.cp_incl_measure[j] = 0;
              *it.second.cp_excl_measure[j] = 0;
            }
          }
        }
      }
    }
  }
}


// Note: this function should be called once per start/stop, else it will double count
void find_per_process_max(MPI_Comm cm){
  if (mode<2){
    int cm_rank; MPI_Comm_rank(cm,&cm_rank);
    double_int buffer[num_per_process_measures];
    for (size_t i=0; i<num_per_process_measures; i++){
      max_per_process_costs[i] = volume_costs[i];
      buffer[i].first          = volume_costs[i];
      buffer[i].second         = cm_rank;
    }
    PMPI_Allreduce(MPI_IN_PLACE, &max_per_process_costs[0], num_per_process_measures, MPI_DOUBLE, MPI_MAX, cm);
    PMPI_Allreduce(MPI_IN_PLACE, &buffer[0], num_per_process_measures, MPI_DOUBLE_INT, MPI_MAXLOC, cm);
    size_t save=0;
    for (size_t i=0; i<num_per_process_measures-2*cost_model_size-1; i++){
      if (breakdown[i] == 0) continue;
      if (cm_rank == buffer[2*cost_model_size+i+1].second){
        for (size_t j=0; j<num_tracker_per_process_measures*list_size; j++){
          max_per_process_costs[num_per_process_measures+save*(num_tracker_per_process_measures*list_size+1)+j] = volume_costs[num_volume_measures+j];
        }
        max_per_process_costs[num_per_process_measures+(save+1)*(num_tracker_per_process_measures*list_size+1)-1] = volume_costs[num_volume_measures-2];
      }
      else{
        for (size_t j=0; j<num_tracker_per_process_measures*list_size+1; j++){
          max_per_process_costs[num_per_process_measures+save*(num_tracker_per_process_measures*list_size+1)+j] = 0.;
        }
      }
      PMPI_Allreduce(MPI_IN_PLACE, &max_per_process_costs[num_per_process_measures+save*(num_tracker_per_process_measures*list_size+1)], num_tracker_per_process_measures*list_size+1, MPI_DOUBLE, MPI_MAX, cm);
      save++;
    }
  }
  // else if mode ==2, I need to add (do later)
}

void compute_volume(MPI_Comm cm){
  if (mode<2){
    PMPI_Allreduce(MPI_IN_PLACE, &volume_costs[0], volume_costs.size(), MPI_DOUBLE, MPI_SUM, cm);
    int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    for (int i=0; i<volume_costs.size(); i++){
      volume_costs[i] /= (1.*world_size);
    }
  }
  else if (mode==2){
    int rank; MPI_Comm_rank(cm,&rank);
    for (size_t i=0; i<num_volume_measures; i++){
      max_per_process_costs[i] = volume_costs[i];
      timer_info_sender[i].first = volume_costs[i];
      timer_info_sender[i].second = rank;
    }
    PMPI_Allreduce(&timer_info_sender[0].first, &timer_info_receiver[0].first, max_per_process_costs.size(), MPI_DOUBLE_INT, MPI_MAXLOC, cm);
    PMPI_Allreduce(MPI_IN_PLACE, &max_per_process_costs[0], max_per_process_costs.size(), MPI_DOUBLE, MPI_MAX, cm);
    PMPI_Allreduce(MPI_IN_PLACE, &volume_costs[0], volume_costs.size(), MPI_DOUBLE, MPI_SUM, cm);
    // Using the rank determining the largest per-process max, each process gets its values
    int perproces_runtime_root_rank = timer_info_receiver[num_volume_measures-1].second;
    int ftimer_size = 0;
    if (rank==perproces_runtime_root_rank){
      ftimer_size = symbol_timers.size();
    }
    PMPI_Allreduce(MPI_IN_PLACE,&ftimer_size,1,MPI_INT,MPI_SUM,cm);

    std::vector<int> symbol_sizes(ftimer_size,0);
    if (rank==perproces_runtime_root_rank){
      int symbol_offset = 0;
      for (auto i=0; i<symbol_timers.size(); i++){
        symbol_sizes[i] = symbol_order[i].size();
        for (auto j=0; j<symbol_sizes[i]; j++){
          symbol_pad[symbol_offset+j] = symbol_order[i][j];
        }
        symbol_offset += symbol_sizes[i];
      }
    }
    PMPI_Allreduce(MPI_IN_PLACE,&symbol_sizes[0],ftimer_size,MPI_INT,MPI_SUM,cm);

    int num_chars = 0;
    for (auto i=0; i<ftimer_size; i++){
      num_chars += symbol_sizes[i];
    }
    if (rank == perproces_runtime_root_rank){
      PMPI_Bcast(&symbol_timer_pad_pp[0],(num_ftimer_measures*num_volume_measures+1)*ftimer_size,MPI_DOUBLE,rank,cm);
      PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,rank,cm);
    }
    else{
      PMPI_Bcast(&symbol_timer_pad_pp[0],(num_ftimer_measures*num_volume_measures+1)*ftimer_size,MPI_DOUBLE,perproces_runtime_root_rank,cm);
      PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,perproces_runtime_root_rank,cm);
      int symbol_offset = 0;
      for (int i=0; i<ftimer_size; i++){
        auto reconstructed_symbol = std::string(symbol_pad.begin()+symbol_offset,symbol_pad.begin()+symbol_offset+symbol_sizes[i]);

        if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
          symbol_timers[reconstructed_symbol] = ftimer(reconstructed_symbol);
          symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
        }
        *symbol_timers[reconstructed_symbol].pp_numcalls = symbol_timer_pad_pp[(num_ftimer_measures*num_volume_measures+1)*i];
        for (int j=0; j<num_volume_measures; j++){
          *symbol_timers[reconstructed_symbol].pp_incl_measure[j] = symbol_timer_pad_pp[(num_ftimer_measures*num_volume_measures+1)*i+2*j+1];
          *symbol_timers[reconstructed_symbol].pp_excl_measure[j] = symbol_timer_pad_pp[(num_ftimer_measures*num_volume_measures+1)*i+2*(j+1)];
        }
        symbol_timers[reconstructed_symbol].has_been_processed = true;
        symbol_offset += symbol_sizes[i];
      }
      // Now cycle through and find the symbols that were not processed and set their accumulated measures to 0
      for (auto& it : symbol_timers){
        if (it.second.has_been_processed){ it.second.has_been_processed = false; }
        else{
          *it.second.pp_numcalls = 0;
          for (int j=0; j<num_volume_measures; j++){
            *it.second.pp_incl_measure[j] = 0;
            *it.second.pp_excl_measure[j] = 0;
          }
        }
      }
    }
  }
}

void tracker::set_header(){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if (*this->my_comm_time != 0){
    std::vector<double> vec(1);
    save_info[this->name] = std::move(vec);
  }
}

void tracker::set_critical_path_costs(size_t idx){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if ((*this->my_comm_time != 0) && (breakdown_size>0)){
    std::vector<double> vec(num_tracker_critical_path_measures);
    int save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]){
        vec[2*save] = *(this->critical_path_wrd_count+idx+save*breakdown_size);
        vec[2*save+1] = *(this->critical_path_msg_count+idx+save*breakdown_size);
        save++;
      }
    }
    vec[num_tracker_critical_path_measures-3] = *(this->critical_path_comm_time+idx);
    vec[num_tracker_critical_path_measures-2] = *(this->critical_path_synch_time+idx);
    vec[num_tracker_critical_path_measures-1] = *(this->critical_path_datamvt_time+idx);
    save_info[this->name] = std::move(vec);
  }
}

void tracker::set_per_process_costs(size_t idx){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if ((*this->my_comm_time != 0) && (breakdown_size>0)){
    std::vector<double> vec(num_tracker_per_process_measures-1);	// don't count idle time
    int save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]){
        vec[2*save] = max_per_process_costs[num_per_process_measures+idx*(num_tracker_per_process_measures*list_size+1)+this->tag*num_tracker_per_process_measures+save];
        vec[2*save+1] = max_per_process_costs[num_per_process_measures+idx*(num_tracker_per_process_measures*list_size+1)+this->tag*num_tracker_per_process_measures+cost_model_size+save];
        save++;
      }
    }
    // For now, do not include idle time
    vec[num_tracker_per_process_measures-4] = max_per_process_costs[num_per_process_measures+idx*(num_tracker_per_process_measures*list_size+1)+this->tag*num_tracker_per_process_measures+num_tracker_per_process_measures-3];
    vec[num_tracker_per_process_measures-3] = max_per_process_costs[num_per_process_measures+idx*(num_tracker_per_process_measures*list_size+1)+this->tag*num_tracker_per_process_measures+num_tracker_per_process_measures-2];
    vec[num_tracker_per_process_measures-2] = max_per_process_costs[num_per_process_measures+idx*(num_tracker_per_process_measures*list_size+1)+this->tag*num_tracker_per_process_measures+num_tracker_per_process_measures-1];
    save_info[this->name] = std::move(vec);
  }
}

void tracker::set_volume_costs(){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if (*this->my_comm_time != 0){
    std::vector<double> vec(num_tracker_volume_measures);
    int save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]){
        vec[2*save] = *(this->my_wrd_count+save);
        vec[2*save+1] = *(this->my_msg_count+save);
        save++;
      }
    }
    vec[2*cost_model_size] = *this->my_bar_time;
    vec[2*cost_model_size+1] = *this->my_comm_time;
    vec[2*cost_model_size+2] = *this->my_synch_time;
    vec[2*cost_model_size+3] = *this->my_datamvt_time;
    save_info[this->name] = std::move(vec);
  }
}

ftimer::ftimer(std::string name_, bool internal){
  auto save_time = MPI_Wtime();
  assert(name_.size() <= max_timer_name_length);
  assert(symbol_timers.size() < max_num_symbols);
  this->name = std::move(name_);
  this->cp_numcalls = &symbol_timer_pad_local_cp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)]; *this->cp_numcalls = 0;
  this->pp_numcalls = &symbol_timer_pad_pp[(symbol_timers.size()-1)*(num_ftimer_measures*num_volume_measures+1)]; *this->pp_numcalls = 0;
  this->vol_numcalls = &symbol_timer_pad_vol[(symbol_timers.size()-1)*(num_ftimer_measures*num_volume_measures+1)]; *this->vol_numcalls = 0;
  for (auto i=0; i<num_critical_path_measures; i++){
    this->cp_incl_measure[i] = &symbol_timer_pad_local_cp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)+2*i+1]; *cp_incl_measure[i] = 0.;
    this->cp_excl_measure[i] = &symbol_timer_pad_local_cp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)+2*(i+1)]; *cp_excl_measure[i] = 0.;
  }
  for (auto i=0; i<num_volume_measures; i++){
    this->pp_incl_measure[i] = &symbol_timer_pad_pp[(symbol_timers.size()-1)*(num_ftimer_measures*num_volume_measures+1)+2*i+1]; *pp_incl_measure[i] = 0.;
    this->pp_excl_measure[i] = &symbol_timer_pad_pp[(symbol_timers.size()-1)*(num_ftimer_measures*num_volume_measures+1)+2*(i+1)]; *pp_excl_measure[i] = 0.;
    this->vol_incl_measure[i] = &symbol_timer_pad_vol[(symbol_timers.size()-1)*(num_ftimer_measures*num_volume_measures+1)+2*i+1]; *vol_incl_measure[i] = 0.;
    this->vol_excl_measure[i] = &symbol_timer_pad_vol[(symbol_timers.size()-1)*(num_ftimer_measures*num_volume_measures+1)+2*(i+1)]; *vol_excl_measure[i] = 0.;
  }
  this->has_been_processed = false;
  if (!internal){
    critical_path_costs[num_critical_path_measures-2] += (save_time - computation_timer);		// update critical path computation time
    critical_path_costs[num_critical_path_measures-1] += (save_time - computation_timer);		// update critical path runtime
    volume_costs[num_volume_measures-2]        += (save_time - computation_timer);		// update local computation time
    volume_costs[num_volume_measures-1]        += (save_time - computation_timer);		// update local runtime
    for (size_t i=0; i<breakdown_size; i++){
      critical_path_costs[critical_path_costs_size-1-i] += (save_time - computation_timer);;
    }
    computation_timer = MPI_Wtime();
  }
}

void ftimer::start(){
  auto save_time = MPI_Wtime();
  if (symbol_stack.size()>0){
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[num_critical_path_measures-1] += (save_time-symbol_timers[symbol_stack.top()].start_timer.top());
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[num_critical_path_measures-2] += (save_time-symbol_timers[symbol_stack.top()].start_timer.top());
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-1] += (save_time-symbol_timers[symbol_stack.top()].start_timer.top());
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-2] += (save_time-symbol_timers[symbol_stack.top()].start_timer.top());
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-1] += (save_time-symbol_timers[symbol_stack.top()].start_timer.top());
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-2] += (save_time-symbol_timers[symbol_stack.top()].start_timer.top());
  }
  symbol_stack.push(this->name);
  this->exclusive_contributions.push(typename decltype(this->exclusive_contributions)::value_type());
// FIX LATER! NOT COMPILING  this->exclusive_measure.push({0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0});
  // Explicit measure will never be finalize until the specific symbol's stack is 0 (nontrivial only for recursive nested symbols)
  critical_path_costs[num_critical_path_measures-2] += (save_time - computation_timer);		// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += (save_time - computation_timer);		// update critical path runtime
  volume_costs[num_volume_measures-2]        += (save_time - computation_timer);		// update local computation time
  volume_costs[num_volume_measures-1]        += (save_time - computation_timer);		// update local runtime
  for (size_t i=0; i<breakdown_size; i++){
    critical_path_costs[critical_path_costs_size-1-i] += (save_time - computation_timer);;
  }
  computation_timer = MPI_Wtime();
  this->start_timer.push(computation_timer);
}

void ftimer::stop(){
  auto save_time = MPI_Wtime();
  assert(this->start_timer.size()>0);
  this->exclusive_measure.top()[num_critical_path_measures-1] += (save_time-this->start_timer.top());
  this->exclusive_measure.top()[num_critical_path_measures-2] += (save_time-this->start_timer.top());
  *this->cp_excl_measure[num_critical_path_measures-1] += (save_time-this->start_timer.top());
  *this->cp_excl_measure[num_critical_path_measures-2] += (save_time-this->start_timer.top());
  *this->pp_excl_measure[num_volume_measures-1] += (save_time-this->start_timer.top());
  *this->pp_excl_measure[num_volume_measures-2] += (save_time-this->start_timer.top());

  for (auto i=0; i<num_critical_path_measures; i++){
    //TODO: I guess at this point I will try to simply not add the exclusive contributions until the top recursive symbol is being stopped. I'm doubtful this approach is best
    if (this->start_timer.size()==1){
      *this->cp_incl_measure[i] += this->exclusive_measure.top()[i];
      *this->pp_incl_measure[i] += this->exclusive_measure.top()[i];
      for (auto& it : this->exclusive_contributions.top()){	// Not the best data access pattern. Loops should ideally be switched, but map probably not large enough to matter
        *this->cp_incl_measure[i] += it.second[i];
        *this->pp_incl_measure[i] += it.second[i];
      }
    }
  }
  *this->cp_numcalls = *this->cp_numcalls + 1.;
  *this->pp_numcalls = *this->pp_numcalls + 1.;
  *this->vol_numcalls = *this->vol_numcalls + 1.;

  assert(this->exclusive_contributions.size()>0);
  assert(this->exclusive_measure.size()>0);
  for (auto i=0; i<num_critical_path_measures; i++){
    this->exclusive_contributions.top()[this->name][i] += this->exclusive_measure.top()[i];
  }

  auto save_excl_contributions = this->exclusive_contributions.top();
  this->exclusive_contributions.pop();
  this->start_timer.pop();
  this->exclusive_measure.pop();
  symbol_stack.pop();

  if (symbol_stack.size() > 0){
    for (auto& it : save_excl_contributions){
      for (auto i=0; i<num_critical_path_measures; i++){
        symbol_timers[symbol_stack.top()].exclusive_contributions.top()[it.first][i] += it.second[i];
      }
    }
  } 

  critical_path_costs[num_critical_path_measures-2] += (save_time - computation_timer);		// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += (save_time - computation_timer);		// update critical path runtime
  volume_costs[num_volume_measures-2]        += (save_time - computation_timer);		// update local computation time
  volume_costs[num_volume_measures-1]        += (save_time - computation_timer);		// update local runtime
  for (size_t i=0; i<breakdown_size; i++){
    critical_path_costs[critical_path_costs_size-1-i] += (save_time - computation_timer);;
  }
  computation_timer = MPI_Wtime();
  if (symbol_stack.size()>0){
    symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer;
  }
}

std::vector<std::string> parse_file_string(){
  std::vector<std::string> inputs;
  auto prev=0;
  auto First=false;
  for (auto i=0; i<file_name.size(); i++){
    if (file_name[i]=='+'){
      if (First){
        inputs.emplace_back(file_name.substr(prev,i-prev));
      }
      else{
        First=true;
      }
      prev=i+1;
    }
  }
  return inputs;
}

template<typename StreamType>
void print_inputs(StreamType& Stream, int np, std::vector<std::string>& inputs){
  Stream << np;
  for (auto input_str : inputs){
    Stream << "\t" << input_str;
  }
}

template<typename StreamType>
void print_cost_model_header(StreamType& Stream){
  if (cost_models[0]){
    Stream << std::left << std::setw(25) << "SimpleCommCost";
    Stream << std::left << std::setw(25) << "SimpleSynchCost";
  }
  if (cost_models[1]){
    Stream << std::left << std::setw(25) << "ABbutterflyCommCost";
    Stream << std::left << std::setw(25) << "ABbutterflySynchCost";
  }
  if (cost_models[2]){
    Stream << std::left << std::setw(25) << "BSPCommCost";
    Stream << std::left << std::setw(25) << "BSPSynchCost";
  }
}

template<typename StreamType>
void print_cost_model_header_file(StreamType& Stream){
  if (cost_models[0]){
    Stream << "\tSimpleCommCost";
    Stream << "\tSimpleSynchCost";
  }
  if (cost_models[1]){
    Stream << "\tABbutterflyCommCost";
    Stream << "\tABbutterflySynchCost";
  }
  if (cost_models[2]){
    Stream << "\tBSPCommCost";
    Stream << "\tBSPSynchCost";
  }
}

template<typename StreamType>
void print_header(StreamType& Stream, size_t num_inputs){
  for (size_t idx = 0; idx < (num_inputs+1); idx++){
    if (idx != 0){
      Stream << "\t";
    }
    Stream << "Input";
  }
  print_cost_model_header_file(Stream);
  Stream << "\tCommunicationTime\tSynchronizationTime\tDataMvtTime\tComputationTime\tRunTime";// critical path
  print_cost_model_header_file(Stream);
  Stream << "\tCommunicationTime\tSynchronizationTime\tDataMvtTime\tComputationTime\tRunTime";// per-process
  print_cost_model_header_file(Stream);
  Stream << "\tIdleTime\tCommunicationTime\tSynchronizationTime\tDataMvtTime\tComputationTime\tRunTime";// volume
  for (auto i=0; i<num_tracker_critical_path_measures*breakdown_size+num_tracker_volume_measures;i++){
    for (auto& it : save_info){
     Stream << "\t" << it.first;
    }
  }
}

void record(std::ofstream& Stream, size_t factor){
  assert(internal_comm_info.size() == 0);
  if (mode<2){
    auto np=0; MPI_Comm_size(MPI_COMM_WORLD,&np);
    if (is_world_root){
      auto inputs = parse_file_string();
      for (int i=0; i<list_size; i++){
        list[i]->set_header();
      }
      if (is_first_iter){
        print_header(Stream,inputs.size());
        Stream << "\n";
      }
      print_inputs(Stream,np,inputs);
      for (size_t i=0; i<num_critical_path_measures; i++){
        Stream << "\t" << factor*critical_path_costs[i];
      }
      for (size_t i=0; i<num_volume_measures; i++){
        if (i==2*cost_model_size) continue;// skip idle time (for now?)
        Stream << "\t" << factor*max_per_process_costs[i];
      }
      for (size_t i=0; i<num_volume_measures; i++){
        Stream << "\t" << factor*volume_costs[i];
      }
      for (int i=0; i<list_size; i++){
        list[i]->set_volume_costs();
      }
      for (size_t j=0; j<num_tracker_volume_measures; j++){
        for (auto& it : save_info){
          Stream << "\t" << factor*it.second[j];
        }
      }
      size_t breakdown_idx=0;
      for (auto i=0; i<num_critical_path_measures-2*cost_model_size; i++){
        if (!breakdown[i]) continue;
        Stream << "\t" << factor*critical_path_costs[critical_path_costs_size-breakdown_size+breakdown_idx-1];
        breakdown_idx++;
      }
      breakdown_idx=0;
      for (auto i=0; i<num_critical_path_measures-2*cost_model_size; i++){
        if (!breakdown[i]) continue;
        // Save the critter information before printing
        for (size_t j=0; j<list_size; j++){
          list[j]->set_critical_path_costs(breakdown_idx);
        }
        breakdown_idx++;
        for (size_t j=0; j<num_tracker_critical_path_measures; j++){
          for (auto& it : save_info){
            Stream << "\t" << factor*it.second[j];
          }
        }
      }
    }
  }
}

void record(std::ostream& Stream, size_t factor){
  assert(internal_comm_info.size() == 0);
  if (mode==0){
    if (is_world_root){
      Stream << std::left << std::setw(25) << "Runtime:";
      Stream << "\n";
      Stream << std::left << std::setw(25) << "                  ";
      Stream << std::left << std::setw(25) << factor*critical_path_costs[num_critical_path_measures-1];
      Stream << "\n\n";
    }
  }
  else if (mode==1){
    if (is_world_root){
      Stream << "\n\n";
      Stream << std::left << std::setw(25) << "Critical path:";
      print_cost_model_header(Stream);
      Stream << std::left << std::setw(25) << "";
      Stream << std::left << std::setw(25) << "CommTime";
      Stream << std::left << std::setw(25) << "SynchTime";
      Stream << std::left << std::setw(25) << "DataMvtTime";
      Stream << std::left << std::setw(25) << "CompTime";
      Stream << std::left << std::setw(25) << "RunTime";
      Stream << "\n";
      Stream << std::left << std::setw(25) << "                  ";
      for (size_t i=0; i<num_critical_path_measures+1; i++){
        if (i==(2*cost_model_size)) Stream << std::left << std::setw(25) << "";
        else if (i<(2*cost_model_size)) Stream << std::left << std::setw(25) << factor*critical_path_costs[i];
        else Stream << std::left << std::setw(25) << factor*critical_path_costs[i-1];
      }
      Stream << "\n\n";

      Stream << std::left << std::setw(25) << "Per-process max:";
      print_cost_model_header(Stream);
      Stream << std::left << std::setw(25) << "IdleTime";
      Stream << std::left << std::setw(25) << "CommTime";
      Stream << std::left << std::setw(25) << "SynchTime";
      Stream << std::left << std::setw(25) << "DataMvtTime";
      Stream << std::left << std::setw(25) << "CompTime";
      Stream << std::left << std::setw(25) << "RunTime";
      Stream << "\n";
      Stream << std::left << std::setw(25) << "                  ";
      for (size_t i=0; i<num_volume_measures; i++){
        Stream << std::left << std::setw(25) << factor*max_per_process_costs[i];
      }
      Stream << "\n\n";

      Stream << std::left << std::setw(25) << "Volume:";
      print_cost_model_header(Stream);
      Stream << std::left << std::setw(25) << "IdleTime";
      Stream << std::left << std::setw(25) << "CommTime";
      Stream << std::left << std::setw(25) << "SynchTime";
      Stream << std::left << std::setw(25) << "DataMvtTime";
      Stream << std::left << std::setw(25) << "CompTime";
      Stream << std::left << std::setw(25) << "RunTime";
      Stream << "\n";
      Stream << std::left << std::setw(25) << "                  ";
      for (size_t i=0; i<num_volume_measures; i++){
        Stream << std::left << std::setw(25) << factor*volume_costs[i];
      }
      Stream << "\n\n";

      size_t breakdown_idx=0;
      for (auto i=0; i<num_critical_path_measures-2*cost_model_size; i++){
        if (!breakdown[i]) continue;
        if (i==0){
          Stream << std::left << std::setw(25) << "CommTime max:";
        } else if (i==1){
          Stream << std::left << std::setw(25) << "SynchTime max:";
        } else if (i==2){
          Stream << std::left << std::setw(25) << "DataMvtTime max:";
        } else if (i==3){
          Stream << std::left << std::setw(25) << "CompTime max:";
        } else if (i==4){
          Stream << std::left << std::setw(25) << "RunTime max:";
        }
        Stream << std::left << std::setw(25) << "MeasureType";
        Stream << std::left << std::setw(25) << "CompTime";
        Stream << std::left << std::setw(25) << "IdleTime";
        print_cost_model_header(Stream);
        Stream << std::left << std::setw(25) << "CommTime";
        Stream << std::left << std::setw(25) << "SynchTime";
        Stream << std::left << std::setw(25) << "DataMvtTime";
        Stream << "\n";
        Stream << std::left << std::setw(25) << "Computation";
        Stream << std::left << std::setw(25) << "path";
        Stream << std::left << std::setw(25) << factor*critical_path_costs[critical_path_costs_size-breakdown_size+breakdown_idx];
        Stream << "\n";
        Stream << std::left << std::setw(25) << "Idle";
        Stream << std::left << std::setw(25) << "path";
        Stream << std::left << std::setw(25) << 0.0;
        Stream << std::left << std::setw(25) << 0.0;
        for (int j=0; j<list_size; j++){
          list[j]->set_critical_path_costs(breakdown_idx);
        }
        for (auto& it : save_info){
          Stream << "\n";
          Stream << std::left << std::setw(25) << it.first;
          Stream << std::left << std::setw(25) << "path";
          Stream << std::left << std::setw(25) << 0.0;
          Stream << std::left << std::setw(25) << 0.0;
          for (size_t j=0; j<num_tracker_critical_path_measures; j++){
            Stream << std::left << std::setw(25) << factor*it.second[j];
          }
        }
        Stream << "\n";
        Stream << std::left << std::setw(25) << "Computation";
        Stream << std::left << std::setw(25) << "per-process";
        Stream << std::left << std::setw(25) << factor*max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+1)-1];
        Stream << "\n";
        Stream << std::left << std::setw(25) << "Idle";
        Stream << std::left << std::setw(25) << "per-process";
        Stream << std::left << std::setw(25) << 0.0;	// TODO: Needs replacing
        Stream << std::left << std::setw(25) << 0.0;
        for (int j=0; j<list_size; j++){
          list[j]->set_per_process_costs(breakdown_idx);
        }
        for (auto& it : save_info){
          Stream << "\n";
          Stream << std::left << std::setw(25) << it.first;
          Stream << std::left << std::setw(25) << "per-process";
          Stream << std::left << std::setw(25) << 0.0;
          Stream << std::left << std::setw(25) << 0.0;
          for (size_t j=0; j<num_tracker_per_process_measures-1; j++){
            Stream << std::left << std::setw(25) << factor*it.second[j];
          }
        }
        breakdown_idx++;
        Stream << "\n\n";
      }
      for (int i=0; i<list_size; i++){
        list[i]->set_volume_costs();
      }
      Stream << std::left << std::setw(25) << "Volume:";
      print_cost_model_header(Stream);
      Stream << std::left << std::setw(25) << "IdleTime";
      Stream << std::left << std::setw(25) << "CommTime";
      Stream << std::left << std::setw(25) << "SynchTime";
      Stream << std::left << std::setw(25) << "DataMvtTime";
      for (auto& it : save_info){
        Stream << "\n";
        Stream << std::left << std::setw(25) << it.first;
        for (size_t j=0; j<num_tracker_volume_measures; j++){
          Stream << std::left << std::setw(25) << factor*it.second[j];
        }
      }
      Stream << "\n";
    }
  }
  else if (mode == 2){
    if (is_world_root){
      Stream << "***********************************************************************************************************************";
      std::vector<std::pair<std::string,std::array<double,6>>> sort_info(symbol_timers.size());
      for (auto i=num_critical_path_measures-1; i>=0; i--){
        sort_info.clear(); sort_info.resize(symbol_timers.size());
        // Reset symbol timers and sort
        size_t j=0;
        for (auto& it : symbol_timers){
          assert(it.second.start_timer.size() == 0);
          sort_info[j++] = std::make_pair(it.second.name,std::array<double,6>{*it.second.cp_numcalls,*it.second.cp_excl_measure[i],*it.second.pp_numcalls,*it.second.pp_excl_measure[(i>=3 ? i+1 : i)],*it.second.vol_numcalls,*it.second.vol_excl_measure[(i>=3 ? i+1 : i)]});
        }
        std::sort(sort_info.begin(),sort_info.end(),[](std::pair<std::string,std::array<double,6>>& vec1, std::pair<std::string,std::array<double,6>>& vec2){return vec1.second[1] > vec2.second[1];});
        // Exclusive
        Stream << "\n\n\n\n" << std::left << std::setw(max_timer_name_length) << critical_path_measure_names[i];
        Stream << std::left << std::setw(15) << "cp-#calls";
        Stream << std::left << std::setw(15) << "cp-excl";
        Stream << std::left << std::setw(15) << "cp-excl %";
        Stream << std::left << std::setw(15) << "pp-#calls";
        Stream << std::left << std::setw(15) << "pp-excl";
        Stream << std::left << std::setw(15) << "pp-excl %";
        Stream << std::left << std::setw(15) << "vol-#calls";
        Stream << std::left << std::setw(15) << "vol-excl";
        Stream << std::left << std::setw(15) << "vol-excl %";
        double cp_total_exclusive = 0.;
        double pp_total_exclusive = 0.;
        double vol_total_exclusive = 0.;
        for (auto& it : sort_info){
          Stream << "\n" << std::left << std::setw(max_timer_name_length) << it.first;
          Stream << std::left << std::setw(15) << it.second[0];
          Stream << std::left << std::setw(15) << it.second[1];
          Stream << std::left << std::setw(15) << std::setprecision(4) << 100.*(critical_path_costs[i] == 0. ? 100.0 : it.second[1]/critical_path_costs[i]);
          Stream << std::left << std::setw(15) << it.second[2];
          Stream << std::left << std::setw(15) << it.second[3];
          Stream << std::left << std::setw(15) << std::setprecision(4) << 100.*(critical_path_costs[i] == 0. ? 0.0 : it.second[3]/critical_path_costs[i]);
          Stream << std::left << std::setw(15) << it.second[4];
          Stream << std::left << std::setw(15) << it.second[5];
          Stream << std::left << std::setw(15) << std::setprecision(4) << 100.*(volume_costs[i] == 0. ? 0.0 : it.second[5]/volume_costs[i]);
          cp_total_exclusive += it.second[1];
          pp_total_exclusive += it.second[3];
          vol_total_exclusive += it.second[5];
        }
        Stream << "\n" << std::left << std::setw(max_timer_name_length) << "total";
        Stream << std::left << std::setw(15) << "";
        Stream << std::left << std::setw(15) << cp_total_exclusive;
        Stream << std::left << std::setw(15) << 100.*cp_total_exclusive/critical_path_costs[i];
        Stream << std::left << std::setw(15) << "";
        Stream << std::left << std::setw(15) << pp_total_exclusive;
        Stream << std::left << std::setw(15) << 100.*pp_total_exclusive/critical_path_costs[i];
        Stream << std::left << std::setw(15) << "";
        Stream << std::left << std::setw(15) << vol_total_exclusive;
        Stream << std::left << std::setw(15) << 100.*vol_total_exclusive/volume_costs[i];
        Stream << "\n";

        // Reset symbol timers and sort
        sort_info.clear(); sort_info.resize(symbol_timers.size());
        j=0;
        for (auto& it : symbol_timers){
          assert(it.second.start_timer.size() == 0);
          sort_info[j++] = std::make_pair(it.second.name,std::array<double,6>{*it.second.cp_numcalls,*it.second.cp_incl_measure[i],*it.second.pp_numcalls,*it.second.pp_incl_measure[(i>=3 ? i+1 : i)],*it.second.vol_numcalls,*it.second.vol_incl_measure[(i>=3 ? i+1 : i)]});
        }
        std::sort(sort_info.begin(),sort_info.end(),[](std::pair<std::string,std::array<double,6>>& vec1, std::pair<std::string,std::array<double,6>>& vec2){return vec1.second[1] > vec2.second[1];});
        // Inclusive
        Stream << "\n" << std::left << std::setw(max_timer_name_length) << critical_path_measure_names[i];
        Stream << std::left << std::setw(15) << "cp-#calls";
        Stream << std::left << std::setw(15) << "cp-incl";
        Stream << std::left << std::setw(15) << "cp-incl %";
        Stream << std::left << std::setw(15) << "pp-#calls";
        Stream << std::left << std::setw(15) << "pp-incl";
        Stream << std::left << std::setw(15) << "pp-incl %";
        Stream << std::left << std::setw(15) << "vol-#calls";
        Stream << std::left << std::setw(15) << "vol-incl";
        Stream << std::left << std::setw(15) << "vol-incl %";
        double cp_total_inclusive = 0.;
        double pp_total_inclusive = 0.;
        double vol_total_inclusive = 0.;
        for (auto& it : sort_info){
          Stream << "\n" << std::left << std::setw(max_timer_name_length) << it.first;
          Stream << std::left << std::setw(15) << it.second[0];
          Stream << std::left << std::setw(15) << it.second[1];
          Stream << std::left << std::setw(15) << std::setprecision(4) << 100.*(critical_path_costs[i] == 0. ? 100.0 : it.second[1]/critical_path_costs[i]);
          Stream << std::left << std::setw(15) << it.second[2];
          Stream << std::left << std::setw(15) << it.second[3];
          Stream << std::left << std::setw(15) << std::setprecision(4) << 100.*(critical_path_costs[i] == 0. ? 100.0 : it.second[3]/critical_path_costs[i]);
          Stream << std::left << std::setw(15) << it.second[4];
          Stream << std::left << std::setw(15) << it.second[5];
          Stream << std::left << std::setw(15) << std::setprecision(4) << 100.*(volume_costs[i] == 0. ? 100.0 : it.second[5]/volume_costs[i]);
          cp_total_inclusive = std::max(it.second[1],cp_total_inclusive);
          pp_total_inclusive = std::max(it.second[3],pp_total_inclusive);
          vol_total_inclusive = std::max(it.second[5],vol_total_inclusive);
        }
        Stream << "\n" << std::left << std::setw(max_timer_name_length) << "total";
        Stream << std::left << std::setw(15) << "";
        Stream << std::left << std::setw(15) << cp_total_inclusive;
        Stream << std::left << std::setw(15) << 100.*cp_total_inclusive/critical_path_costs[i];
        Stream << std::left << std::setw(15) << "";
        Stream << std::left << std::setw(15) << pp_total_inclusive;
        Stream << std::left << std::setw(15) << 100.*pp_total_inclusive/critical_path_costs[i];
        Stream << std::left << std::setw(15) << "";
        Stream << std::left << std::setw(15) << vol_total_inclusive;
        Stream << std::left << std::setw(15) << 100.*vol_total_inclusive/volume_costs[i];
        Stream << "\n";
      }
    }
  }
}
};

void start(size_t mode){
  assert(mode>=0 && mode < 3);
  assert(internal::internal_comm_info.size() == 0);
  internal::wait_id=true;
  internal::mode=mode;
  internal::critical_path_measure_names = {"comm time","synch time","datamvt time","comp time","runtime"};
  for (int i=0; i<internal::list_size; i++){
    internal::list[i]->init();
  }
  if (internal::is_world_root){
    if (!internal::is_first_iter){
      if (internal::flag) {internal::stream << "\n";} else {std::cout << "\n";}
    }
  }
  for (auto i=0; i<internal::critical_path_costs.size(); i++){
    internal::critical_path_costs[i]=0.;
  }
  for (auto i=0; i<internal::volume_costs.size(); i++){
    internal::volume_costs[i]=0.;
  }
  /*Initiate new timer*/
  internal::computation_timer=MPI_Wtime();
}

void stop(size_t mode, size_t factor){
  volatile double last_time = MPI_Wtime();
//  assert(mode==internal::mode);
  assert(internal::internal_comm_info.size() == 0);
  internal::critical_path_costs[internal::num_critical_path_measures-2]+=(last_time-internal::computation_timer);	// update critical path computation time
  internal::critical_path_costs[internal::num_critical_path_measures-1]+=(last_time-internal::computation_timer);	// update critical path runtime
  internal::volume_costs[internal::num_volume_measures-2]+=(last_time-internal::computation_timer);	// update computation time volume
  internal::volume_costs[internal::num_volume_measures-1]+=(last_time-internal::computation_timer);	// update runtime volume
  for (size_t i=0; i<breakdown_size; i++){
    internal::critical_path_costs[internal::critical_path_costs_size-1-i] += (last_time-internal::computation_timer);
  }

  internal::complete_propagation(MPI_COMM_WORLD,false,-1,-1);
  internal::find_per_process_max(MPI_COMM_WORLD);
  internal::compute_volume(MPI_COMM_WORLD);

  internal::record(std::cout,factor);
  if (internal::flag) {internal::record(internal::stream,factor);}

  internal::wait_id=false;
  internal::is_first_iter = false;
  internal::mode=0;
  internal::save_info.clear();
  for (auto i=0; i<internal::critical_path_costs.size(); i++){
    internal::critical_path_costs[i]=0.;
  }
  for (auto i=0; i<internal::volume_costs.size(); i++){
    internal::volume_costs[i]=0.;
  }
  internal::need_new_line=false;
  internal::symbol_timers.clear();
}
};
