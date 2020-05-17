#include "critter.h"

namespace critter{
namespace internal{

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
blocking _MPI_Ssend("MPI_Ssend",13,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);}
                      );
blocking _MPI_Sendrecv("MPI_Sendrecv",14,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);},
                          [](int64_t n, int p){
                            return std::pair<double,double>(1.,n);}
                      );
blocking _MPI_Sendrecv_replace("MPI_Sendrecv_replace",15,
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
bool flag,is_world_root,is_first_iter,need_new_line,print_volume_symbol;
size_t mode;

double computation_timer;
std::map<MPI_Request,bool> internal_comm_info;
std::map<MPI_Request,std::pair<MPI_Comm,int>> internal_comm_comm;
std::map<MPI_Request,std::pair<double,double>> internal_comm_data;
std::map<MPI_Request,nonblocking*> internal_comm_track;
std::vector<std::pair<double*,int>> internal_comm_prop;
std::vector<MPI_Request> internal_comm_prop_req;
std::vector<int*> internal_timer_prop_int;
std::vector<double*> internal_timer_prop_double;
std::vector<char*> internal_timer_prop_char;
std::vector<MPI_Request> internal_timer_prop_req;
bool decisions[breakdown_size];
std::array<double,critical_path_costs_size> critical_path_costs;
std::array<double,per_process_costs_size> max_per_process_costs;
std::array<double,volume_costs_size> volume_costs;
std::map<std::string,std::vector<double>> save_info;
double new_cs[critical_path_costs_size];
double scratch_pad;
std::vector<char> synch_pad_send;
std::vector<char> synch_pad_recv;
std::vector<char> barrier_pad_send;
std::vector<char> barrier_pad_recv;
std::array<char,max_timer_name_length*max_num_symbols> symbol_pad;
std::array<int,max_num_symbols> symbol_len_pad;
std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_local_cp;
std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_global_cp;
std::array<double,(num_ftimer_measures*num_volume_measures+1)*max_num_symbols> symbol_timer_pad_local_pp;
std::array<double,(num_ftimer_measures*num_volume_measures+1)*max_num_symbols> symbol_timer_pad_global_pp;
std::array<double,(num_ftimer_measures*num_volume_measures+1)*max_num_symbols> symbol_timer_pad_vol;
std::unordered_map<std::string,ftimer> symbol_timers;
std::stack<std::string> symbol_stack;
std::array<std::string,max_num_symbols> symbol_order;
double_int timer_info_sender[num_critical_path_measures];
double_int timer_info_receiver[num_critical_path_measures];
bool wait_id,waitall_id;
double waitall_comp_time;

void tracker::init(){
  this->last_start_time  = -1.;
  this->save_comp_time   = 0.;
}

void tracker::set_cost_pointers(){
  size_t volume_costs_idx        = num_volume_measures+this->tag*num_tracker_volume_measures;
  this->my_wrd_count             = cost_model_size>0 ? &volume_costs[volume_costs_idx] : &scratch_pad;
  this->my_msg_count             = cost_model_size>0 ? &volume_costs[volume_costs_idx+cost_model_size] : &scratch_pad;
  this->my_comm_time             = &volume_costs[volume_costs_idx+2*cost_model_size];
  this->my_synch_time            = &volume_costs[volume_costs_idx+2*cost_model_size+1];
  this->my_datamvt_time          = &volume_costs[volume_costs_idx+2*cost_model_size+2];
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

blocking::blocking(std::string name_, int tag, std::function<std::pair<double,double>(int64_t,int)> cost_func_bsp,
                                             std::function<std::pair<double,double>(int64_t,int)> cost_func_alphabeta_butterfly){
  this->cost_func_bsp              = cost_func_bsp;
  this->cost_func_alphabeta_butterfly = cost_func_alphabeta_butterfly;
  this->name = std::move(name_);
  this->tag = tag;
  this->set_cost_pointers();
  this->init();
}

blocking::blocking(blocking const& t){
  this->cost_func_bsp                 = t.cost_func_bsp;
  this->cost_func_alphabeta_butterfly = t.cost_func_alphabeta_butterfly;
  this->name = t.name;
  this->tag = t.tag;
  this->set_cost_pointers();
  this->init();
}

nonblocking::nonblocking(std::string name_, int tag, std::function<std::pair<double,double>(int64_t,int)> cost_func_bsp,
                                             std::function<std::pair<double,double>(int64_t,int)> cost_func_alphabeta_butterfly){
  this->cost_func_bsp                 = cost_func_bsp;
  this->cost_func_alphabeta_butterfly = cost_func_alphabeta_butterfly;
  this->name = std::move(name_);
  this->tag = tag;
  this->set_cost_pointers();
  this->init();
}

nonblocking::nonblocking(nonblocking const& t){
  this->cost_func_bsp                 = t.cost_func_bsp;
  this->cost_func_alphabeta_butterfly = t.cost_func_alphabeta_butterfly;
  this->name = t.name;
  this->tag = t.tag;
  this->set_cost_pointers();
  this->init();
}

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
  assert(len == critical_path_costs_size);	// this assert prevents user from obtaining wrong output if MPI implementation cuts up the message.
  if (breakdown_size > 0){
    size_t breakdown_idx=0;
    for (int i=0; i<num_critical_path_measures; i++){
      if (breakdown[i]){ decisions[breakdown_idx++] = inout[i] > in[i]; }
    }
    for (int i=0; i<num_critical_path_measures; i++){
      inout[i] = std::max(inout[i],in[i]);
    }
    for (int i=num_critical_path_measures; i<critical_path_costs_size; i++){
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

void complete_timers(double* remote_path_data, size_t msg_id){
  int* envelope_int[2] = { internal_timer_prop_int[4*msg_id+2], internal_timer_prop_int[4*msg_id+3] };
  double* envelope_double[3] = { remote_path_data, internal_timer_prop_double[4*msg_id+2], internal_timer_prop_double[4*msg_id+3] };
  char* envelope_char = internal_timer_prop_char[2*msg_id+1];
  if (envelope_double[0][num_critical_path_measures-1] > critical_path_costs[num_critical_path_measures-1]){
    int ftimer_size = *envelope_int[0];
    int symbol_offset = 0;
    for (int i=0; i<ftimer_size; i++){
      auto reconstructed_symbol = std::string(envelope_char+symbol_offset,envelope_char+symbol_offset+envelope_int[1][i]);
      if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
        symbol_timers[reconstructed_symbol] = ftimer(reconstructed_symbol);
        symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
      }
      *symbol_timers[reconstructed_symbol].cp_numcalls = envelope_double[1][(num_ftimer_measures*num_critical_path_measures+1)*i];
      for (int j=0; j<num_critical_path_measures; j++){
        *symbol_timers[reconstructed_symbol].cp_incl_measure[j] = envelope_double[1][(num_ftimer_measures*num_critical_path_measures+1)*i+2*j+1];
        *symbol_timers[reconstructed_symbol].cp_excl_measure[j] = envelope_double[1][(num_ftimer_measures*num_critical_path_measures+1)*i+2*(j+1)];
      }
      symbol_timers[reconstructed_symbol].cp_exclusive_measure.fill(0.0);
      for (int j=0; j<num_critical_path_measures; j++){
        symbol_timers[reconstructed_symbol].cp_exclusive_contributions[j] = envelope_double[2][i*num_critical_path_measures+j];
      }
      for (int k=0; k<num_critical_path_measures; k++){
        symbol_timers[reconstructed_symbol].cp_exclusive_measure[k] = envelope_double[2][ftimer_size*num_critical_path_measures+i*num_critical_path_measures+k];
      }
      symbol_timers[reconstructed_symbol].has_been_processed = true;
      symbol_offset += envelope_int[1][i];
    }
  }
}

void complete_path_update(){
  PMPI_Waitall(internal_comm_prop_req.size(), &internal_comm_prop_req[0], MPI_STATUSES_IGNORE);
  if (mode>=2) { PMPI_Waitall(internal_timer_prop_req.size(), &internal_timer_prop_req[0], MPI_STATUSES_IGNORE); }
  size_t msg_id=0;
  for (auto& it : internal_comm_prop){
    if (!it.second){
      if (mode>=2) complete_timers(it.first,msg_id++);
      update_critical_path(it.first,&critical_path_costs[0],critical_path_costs_size);
    }
    free(it.first);
  }
  internal_comm_prop.clear(); internal_comm_prop_req.clear();
  for (auto& it : internal_timer_prop_int){ free(it); }
  for (auto& it : internal_timer_prop_double){ free(it); }
  for (auto& it : internal_timer_prop_char){ free(it); }
  internal_timer_prop_int.clear(); internal_timer_prop_double.clear(); internal_timer_prop_char.clear(); internal_timer_prop_req.clear();
}

void blocking::start(volatile double curTime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, bool is_sender, int partner1, int partner2){
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  this->save_comp_time    = curTime - computation_timer;
  critical_path_costs[num_critical_path_measures-2] += this->save_comp_time;	// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += this->save_comp_time;	// update critical path runtime
  volume_costs[num_volume_measures-2]        += this->save_comp_time;		// update local computation time
  volume_costs[num_volume_measures-1]        += this->save_comp_time;		// update local runtime
  for (size_t i=0; i<breakdown_size; i++){ critical_path_costs[critical_path_costs_size-1-i] += this->save_comp_time; }
  if (mode>=2){
    auto last_symbol_time = curTime - symbol_timers[symbol_stack.top()].start_timer.top();
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-2] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-2] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-1] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-2] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-2] += last_symbol_time;
  }

  int el_size,p,rank;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(cm, &p);
  MPI_Comm_rank(cm, &rank);
  this->last_nbytes = nbytes;
  this->last_cm = cm;
  this->last_p = p;
  this->last_is_sender = is_sender;
  this->last_partner1 = partner1;
  this->last_partner2 = partner2;

  volatile double init_time = MPI_Wtime();
  if (this->last_partner1 == -1){ PMPI_Barrier(cm); }
  else {
    char sbuf='e'; char rbuf='d';
    if ((tag < 14) || (tag > 15)){
      if ((is_sender) && (rank != this->last_partner1)){
        PMPI_Ssend(&sbuf, 1, MPI_CHAR, this->last_partner1, internal_tag3, cm);
      } else if ((!is_sender) && (rank != this->last_partner1)){
        PMPI_Recv(&rbuf, 1, MPI_CHAR, this->last_partner1, internal_tag3, cm, MPI_STATUS_IGNORE);
      }
    } else{
        PMPI_Sendrecv(&sbuf, 1, MPI_CHAR, this->last_partner1, internal_tag3, &rbuf, 1, MPI_CHAR, this->last_partner1, internal_tag3, cm, MPI_STATUS_IGNORE);
    }
/*    if ((is_sender) && (this->last_partner2 != -1)){
      PMPI_Recv(&rbuf, 1, MPI_CHAR, this->last_partner2, internal_tag3, cm, MPI_STATUS_IGNORE);
    }
*/
  }
  this->last_barrier_time = MPI_Wtime() - init_time;
  this->last_start_time = MPI_Wtime();
}

void blocking::intermediate(){
  // start synchronization timer for communication routine
  volatile double synchTime = MPI_Wtime();
  this->last_synch_time = synchTime-this->last_start_time;
  // start communication timer for communication routine
  this->last_start_time = MPI_Wtime();
}

// Used only for p2p communication. All blocking collectives use sychronous protocol
void blocking::stop(){
  volatile double new_time = MPI_Wtime();
  volatile double dt = new_time - this->last_start_time;	// complete communication time
  double datamvt_time = std::max(0.,(dt-this->last_synch_time));
  std::pair<double,double> dcost_bsp    = this->cost_func_bsp(this->last_nbytes, this->last_p);
  std::pair<double,double> dcost_alphabeta_butterfly = this->cost_func_alphabeta_butterfly(this->last_nbytes, this->last_p);
  std::vector<std::pair<double,double>> dcosts = {dcost_bsp,dcost_alphabeta_butterfly};

  if (mode>=1){
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
  if (mode>=2){
    // update all communication-related measures for the top symbol in stack
    size_t save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]){
        symbol_timers[symbol_stack.top()].cp_exclusive_measure[save] += dcosts[j].second;
        symbol_timers[symbol_stack.top()].cp_exclusive_measure[cost_models.size()+save] += dcosts[j].first;
        symbol_timers[symbol_stack.top()].pp_exclusive_measure[save] += dcosts[j].second;
        symbol_timers[symbol_stack.top()].pp_exclusive_measure[cost_models.size()+save] += dcosts[j].first;
        *symbol_timers[symbol_stack.top()].cp_excl_measure[save] += dcosts[j].second;
        *symbol_timers[symbol_stack.top()].cp_excl_measure[cost_models.size()+save] += dcosts[j].first;
        *symbol_timers[symbol_stack.top()].pp_excl_measure[save] += dcosts[j].second;
        *symbol_timers[symbol_stack.top()].pp_excl_measure[cost_models.size()+save] += dcosts[j].first;
      }
      save++;
    }
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-5] += dt;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-4] += this->last_synch_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-3] += datamvt_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-1] += dt;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-6] += this->last_barrier_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-6] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-5] += dt;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-4] += this->last_synch_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-3] += datamvt_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += (dt+this->last_barrier_time);
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-5] += dt;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-4] += this->last_synch_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-3] += datamvt_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-1] += dt;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-6] += this->last_barrier_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-6] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-5] += dt;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-4] += this->last_synch_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-3] += datamvt_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += (dt+this->last_barrier_time);
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
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
  volume_costs[num_volume_measures-6] += this->last_barrier_time;	// update local barrier/idle time
  volume_costs[num_volume_measures-5] += dt;				// update local communication time (not volume until after the completion of the program)
  volume_costs[num_volume_measures-4] += this->last_synch_time;		// update local synchronization time
  volume_costs[num_volume_measures-3] += datamvt_time;			// update local data mvt time
  volume_costs[num_volume_measures-1] += (this->last_barrier_time+dt);	// update local runtime with idle time and comm time
  volume_costs[num_volume_measures-6] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);

  if (mode>=2){
    // Special handling of excessively large idle time caused by suspected tool interference
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-6] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-6] -= std::max(0.,volume_costs[num_volume_measures-1]-critical_path_costs[num_critical_path_measures-1]);
  }

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

  // Propogate critical paths for all processes in communicator based on what each process has seen up until now (not including this communication)
  propagate(this->last_cm, this->tag, this->last_is_sender, this->last_partner1, this->last_partner2);
  // Prepare to leave interception and re-enter user code
  this->last_start_time = MPI_Wtime();
  computation_timer = this->last_start_time;
  if (mode>=2){ symbol_timers[symbol_stack.top()].start_timer.top() = this->last_start_time; }
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
  if (mode>=2){
    assert(symbol_stack.size()>0);
    assert(symbol_timers[symbol_stack.top()].start_timer.size()>0);
    this->save_time = curTime - symbol_timers[symbol_stack.top()].start_timer.top();
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-1] += this->save_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-2] += this->save_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += this->save_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-2] += this->save_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-1] += this->save_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-2] += this->save_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-1] += this->save_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_volume_measures-2] += this->save_time;
  }

  int el_size,p;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(cm, &p);
  int rank; MPI_Comm_rank(cm, &rank);

  internal_comm_info[*request] = is_sender;
  internal_comm_comm[*request] = std::make_pair(cm,partner);
  internal_comm_data[*request] = std::make_pair((double)nbytes,(double)p);
  internal_comm_track[*request] = this;

  this->last_start_time = MPI_Wtime();
  computation_timer = this->last_start_time;
  if (mode>=2){ symbol_timers[symbol_stack.top()].start_timer.top() = this->last_start_time; }
}

void nonblocking::stop(MPI_Request* request, double comp_time, double comm_time){
  auto comm_info_it = internal_comm_info.find(*request);
  auto comm_comm_it = internal_comm_comm.find(*request);
  auto comm_data_it = internal_comm_data.find(*request);
  auto comm_track_it = internal_comm_track.find(*request);
  assert(comm_info_it != internal_comm_info.end());
  assert(comm_comm_it != internal_comm_comm.end());
  assert(comm_data_it != internal_comm_data.end());
  assert(comm_track_it != internal_comm_track.end());

  this->last_is_sender = comm_info_it->second;
  this->last_cm = comm_comm_it->second.first;
  this->last_partner1 = comm_comm_it->second.second;
  this->last_partner2 = -1;
  this->last_nbytes = comm_data_it->second.first;
  this->last_p = comm_data_it->second.second;
  this->last_synch_time=0;
  
  // Both sender and receiver will now update its critical path with the data from the communication
  std::pair<double,double> dcost_bsp  = this->cost_func_bsp(this->last_nbytes,this->last_p);
  std::pair<double,double> dcost_alphabeta_butterfly = this->cost_func_alphabeta_butterfly(this->last_nbytes,this->last_p);
  if ((this->tag<20) && (wait_id)) dcost_bsp.first=1.;	// this is usually zero, but we force it to be 1 in special circumstances (for nonblocking p2p with wait_id one)
  std::vector<std::pair<double,double>> dcosts = {dcost_bsp,dcost_alphabeta_butterfly};

  if (mode >= 1){
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

  if (mode>=2){
    size_t save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]){
        symbol_timers[symbol_stack.top()].cp_exclusive_measure[save] += dcosts[j].second;
        symbol_timers[symbol_stack.top()].cp_exclusive_measure[cost_models.size()+save] += dcosts[j].first;
        symbol_timers[symbol_stack.top()].pp_exclusive_measure[save] += dcosts[j].second;
        symbol_timers[symbol_stack.top()].pp_exclusive_measure[cost_models.size()+save] += dcosts[j].first;
        *symbol_timers[symbol_stack.top()].cp_excl_measure[save] += dcosts[j].second;
        *symbol_timers[symbol_stack.top()].cp_excl_measure[cost_models.size()+save] += dcosts[j].first;
        *symbol_timers[symbol_stack.top()].pp_excl_measure[save] += dcosts[j].second;
        *symbol_timers[symbol_stack.top()].pp_excl_measure[cost_models.size()+save] += dcosts[j].first;
      }
      save++;
    }
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-5] += comm_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-4] += 0;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-3] += comm_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-2] += comp_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-1] += (comp_time+comm_time);
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-5] += comm_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-4] += 0;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-3] += comm_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-2] += comp_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += (comp_time+comm_time);
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-5] += comm_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-4] += 0;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-3] += comm_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-2] += comp_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-1] += (comp_time+comm_time);
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-5] += comm_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-4] += 0;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-3] += comm_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-2] += comp_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += (comp_time+comm_time);
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

  if (mode>=1){
    propagate(this->last_cm,this->tag,this->last_is_sender,this->last_partner1,this->last_partner2);
  }
  internal_comm_info.erase(*request);
  internal_comm_comm.erase(*request);
  internal_comm_data.erase(*request);
  internal_comm_track.erase(*request);

  this->last_start_time = MPI_Wtime();
}

void initiate_timers(int rank, int tag, MPI_Comm cm, int partner1, int partner2){
  assert(mode>=2);
  MPI_Request internal_request[10];
  int* send_envelope1 = nullptr; int* send_envelope2 = nullptr; double* send_envelope3 = nullptr; double* send_envelope4 = nullptr; char* send_envelope5 = nullptr;
  int* recv_envelope1 = nullptr; int* recv_envelope2 = nullptr; double* recv_envelope3 = nullptr; double* recv_envelope4 = nullptr; char* recv_envelope5 = nullptr;
  int ftimer_size = symbol_timers.size();
  int num_chars = 0;
  for (int i=0; i<ftimer_size; i++) { num_chars += symbol_order[i].size(); }
  send_envelope1 = (int*)malloc(sizeof(int)); *send_envelope1 = ftimer_size;
  send_envelope2 = (int*)malloc(sizeof(int)*(ftimer_size));
  send_envelope3 = (double*)malloc(sizeof(double)*(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size);
  send_envelope4 = (double*)malloc(sizeof(double)*2*ftimer_size*num_critical_path_measures);
  send_envelope5 = (char*)malloc(sizeof(char)*num_chars);
  int symbol_offset = 0;
  for (auto i=0; i<ftimer_size; i++){
    send_envelope2[i] = symbol_order[i].size();
    for (auto j=0; j<symbol_order[i].size(); j++){
      send_envelope5[symbol_offset+j] = symbol_order[i][j];
    }
    for (auto j=0; j<num_critical_path_measures; j++){
      send_envelope4[i*num_critical_path_measures+j] = symbol_timers[symbol_order[i]].cp_exclusive_contributions[j];
    }
    for (auto k=0; k<num_critical_path_measures; k++){
      send_envelope4[ftimer_size*num_critical_path_measures + i*num_critical_path_measures+k] = symbol_timers[symbol_order[i]].cp_exclusive_measure[k];
    }
    symbol_offset += symbol_order[i].size();
  }
  for (int i=0; i<(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size; i++){ send_envelope3[i] = symbol_timer_pad_local_cp[i]; }
  PMPI_Isend(&send_envelope1[0],1,MPI_INT,partner1,internal_tag,cm,&internal_request[0]);
  PMPI_Isend(&send_envelope2[0],ftimer_size,MPI_INT,partner1,internal_tag,cm,&internal_request[1]);
  PMPI_Isend(&send_envelope3[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,partner1,internal_tag,cm,&internal_request[2]);
  PMPI_Isend(&send_envelope4[0],2*ftimer_size*num_critical_path_measures,MPI_DOUBLE,partner1,internal_tag,cm,&internal_request[3]);
  PMPI_Isend(&send_envelope5[0],symbol_offset,MPI_CHAR,partner1,internal_tag,cm,&internal_request[4]);

  recv_envelope1 = (int*)malloc(sizeof(int));
  recv_envelope2 = (int*)malloc(sizeof(int)*(max_num_symbols));
  recv_envelope3 = (double*)malloc(sizeof(double)*(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols);
  recv_envelope4 = (double*)malloc(sizeof(double)*2*max_num_symbols*num_critical_path_measures);
  recv_envelope5 = (char*)malloc(sizeof(char)*max_timer_name_length*max_num_symbols);
  PMPI_Irecv(recv_envelope1,1,MPI_INT,partner1,internal_tag,cm,&internal_request[5]);
  PMPI_Irecv(recv_envelope2,max_num_symbols,MPI_INT,partner1,internal_tag,cm,&internal_request[6]);
  PMPI_Irecv(recv_envelope3,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols,MPI_DOUBLE,partner1,internal_tag,cm,&internal_request[7]);
  PMPI_Irecv(recv_envelope4,2*max_num_symbols*num_critical_path_measures,MPI_DOUBLE,partner1,internal_tag,cm,&internal_request[8]);
  PMPI_Irecv(recv_envelope5,max_timer_name_length*max_num_symbols,MPI_CHAR,partner1,internal_tag,cm,&internal_request[9]);

  for (int i=0; i<10; i++) { internal_timer_prop_req.push_back(internal_request[i]); }
  internal_timer_prop_int.push_back(send_envelope1); internal_timer_prop_int.push_back(send_envelope2);
  internal_timer_prop_int.push_back(recv_envelope1); internal_timer_prop_int.push_back(recv_envelope2);
  internal_timer_prop_double.push_back(send_envelope3); internal_timer_prop_double.push_back(send_envelope4);
  internal_timer_prop_double.push_back(recv_envelope3); internal_timer_prop_double.push_back(recv_envelope4);
  internal_timer_prop_char.push_back(send_envelope5); internal_timer_prop_char.push_back(recv_envelope5);
}

void propagate_timers(int rank, int tag, MPI_Comm cm, int partner1, int partner2){
  assert(mode>=2);
  int critical_path_runtime_root_rank = timer_info_receiver[num_critical_path_measures-1].second;
  int ftimer_size = 0;
  if (rank==critical_path_runtime_root_rank){ ftimer_size = symbol_timers.size(); }
  if (partner1 == -1){ PMPI_Allreduce(MPI_IN_PLACE,&ftimer_size,1,MPI_INT,MPI_SUM,cm); }
  else{
    int cm_size; MPI_Comm_size(cm, &cm_size); if (critical_path_runtime_root_rank>= cm_size) { std::cout << cm_size << " " << cm << " " << MPI_COMM_WORLD << " " << rank << " " << partner1 << " " << partner2 << " " << critical_path_runtime_root_rank << " " << tag << std::endl; }
    if (rank==critical_path_runtime_root_rank){ PMPI_Send(&ftimer_size,1,MPI_INT,partner1,internal_tag5,cm);
                                                /*if (partner1 != partner2) { PMPI_Send(&ftimer_size,1,MPI_INT,partner2,internal_tag1,cm); }*/ }
    else { PMPI_Recv(&ftimer_size,1,MPI_INT,critical_path_runtime_root_rank,internal_tag5,cm,MPI_STATUS_IGNORE); }
  }
  symbol_len_pad.fill(0);
  std::vector<double> cp_data(2*ftimer_size*num_critical_path_measures,0);
  if (rank==critical_path_runtime_root_rank){
    int symbol_offset = 0;
    for (auto i=0; i<symbol_timers.size(); i++){
      symbol_len_pad[i] = symbol_order[i].size();
      for (auto j=0; j<symbol_len_pad[i]; j++){
        symbol_pad[symbol_offset+j] = symbol_order[i][j];
      }
      for (auto j=0; j<num_critical_path_measures; j++){
        cp_data[i*num_critical_path_measures+j] = symbol_timers[symbol_order[i]].cp_exclusive_contributions[j];
      }
      for (auto k=0; k<num_critical_path_measures; k++){
        cp_data[ftimer_size*num_critical_path_measures + i*num_critical_path_measures+k] = symbol_timers[symbol_order[i]].cp_exclusive_measure[k];
      }
      symbol_offset += symbol_len_pad[i];
    }
  }
  if (partner1 == -1){ PMPI_Allreduce(MPI_IN_PLACE,&symbol_len_pad[0],ftimer_size,MPI_INT,MPI_SUM,cm); }
  else{
    if (rank==critical_path_runtime_root_rank){ PMPI_Send(&symbol_len_pad[0],ftimer_size,MPI_INT,partner1,internal_tag,cm);
                                                if (partner1 != partner2) PMPI_Send(&symbol_len_pad[0],ftimer_size,MPI_INT,partner2,internal_tag,cm); }
    else{ PMPI_Recv(&symbol_len_pad[0],ftimer_size,MPI_INT,critical_path_runtime_root_rank,internal_tag,cm,MPI_STATUS_IGNORE); }
  }
  int num_chars = 0; int num_records = 0;
  for (auto i=0; i<ftimer_size; i++){ num_chars += symbol_len_pad[i]; }
  if (rank == critical_path_runtime_root_rank){
    if (partner1 == -1){
      PMPI_Bcast(&symbol_timer_pad_local_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,rank,cm);
      PMPI_Bcast(&cp_data[0],2*ftimer_size*num_critical_path_measures,MPI_DOUBLE,rank,cm);
      PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,rank,cm);
    }
    else{
      PMPI_Send(&symbol_timer_pad_local_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,partner1,internal_tag1,cm);
      PMPI_Send(&cp_data[0],2*ftimer_size*num_critical_path_measures,MPI_DOUBLE,partner1,internal_tag1,cm);
      PMPI_Send(&symbol_pad[0],num_chars,MPI_CHAR,partner1,internal_tag1,cm);
      if (partner1 != partner2) PMPI_Send(&symbol_timer_pad_local_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,partner2,internal_tag1,cm);
      if (partner1 != partner2) PMPI_Send(&cp_data[0],2*ftimer_size*num_critical_path_measures,MPI_DOUBLE,partner2,internal_tag1,cm);
      if (partner1 != partner2) PMPI_Send(&symbol_pad[0],num_chars,MPI_CHAR,partner2,internal_tag1,cm);
    }
  }
  else{
    if (partner1 == -1){
      PMPI_Bcast(&symbol_timer_pad_global_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,critical_path_runtime_root_rank,cm);
      PMPI_Bcast(&cp_data[0],2*ftimer_size*num_critical_path_measures,MPI_DOUBLE,critical_path_runtime_root_rank,cm);
      PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,critical_path_runtime_root_rank,cm);
    }
    else{
      PMPI_Recv(&symbol_timer_pad_global_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,critical_path_runtime_root_rank,internal_tag1,cm,MPI_STATUS_IGNORE);
      PMPI_Recv(&cp_data[0],2*ftimer_size*num_critical_path_measures,MPI_DOUBLE,critical_path_runtime_root_rank,internal_tag1,cm,MPI_STATUS_IGNORE);
      PMPI_Recv(&symbol_pad[0],num_chars,MPI_CHAR,critical_path_runtime_root_rank,internal_tag1,cm,MPI_STATUS_IGNORE);
    }
    if (rank != partner1){
      int symbol_offset = 0;
      for (int i=0; i<ftimer_size; i++){
        auto reconstructed_symbol = std::string(symbol_pad.begin()+symbol_offset,symbol_pad.begin()+symbol_offset+symbol_len_pad[i]);
        if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
          symbol_timers[reconstructed_symbol] = ftimer(reconstructed_symbol);
          symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
        }
        *symbol_timers[reconstructed_symbol].cp_numcalls = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i];
        for (int j=0; j<num_critical_path_measures; j++){
          *symbol_timers[reconstructed_symbol].cp_incl_measure[j] = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i+2*j+1];
          *symbol_timers[reconstructed_symbol].cp_excl_measure[j] = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i+2*(j+1)];
        }
        symbol_timers[reconstructed_symbol].cp_exclusive_measure.fill(0.0);
        for (int j=0; j<num_critical_path_measures; j++){
          symbol_timers[reconstructed_symbol].cp_exclusive_contributions[j] = cp_data[i*num_critical_path_measures+j];
        }
        for (int k=0; k<num_critical_path_measures; k++){
          symbol_timers[reconstructed_symbol].cp_exclusive_measure[k] = cp_data[ftimer_size*num_critical_path_measures+i*num_critical_path_measures+k];
        }
        symbol_timers[reconstructed_symbol].has_been_processed = true;
        symbol_offset += symbol_len_pad[i];
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

void propagate(MPI_Comm cm, int tag, bool is_sender, int partner1, int partner2){
  int rank; MPI_Comm_rank(cm,&rank);
  if (mode>=2){
    for (int i=0; i<num_critical_path_measures; i++){
      timer_info_sender[i].first = critical_path_costs[i];
      timer_info_sender[i].second = rank;
    }
    if (partner1 == -1){
      if (tag < 20){
        PMPI_Allreduce(&timer_info_sender[0].first, &timer_info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, MPI_MAXLOC, cm);
      }
    } else{
      if (tag < 18){
        if (rank != partner1){
          PMPI_Sendrecv(&timer_info_sender[0].first, num_critical_path_measures, MPI_DOUBLE_INT, partner1, internal_tag,
                        &timer_info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, partner1, internal_tag, cm, MPI_STATUS_IGNORE);
        }
        for (int i=0; i<num_critical_path_measures; i++){
          if (timer_info_sender[i].first>timer_info_receiver[i].first){timer_info_receiver[i].second = rank;}
          else if (timer_info_sender[i].first==timer_info_receiver[i].first){
            if (timer_info_sender[i].second < timer_info_receiver[i].second){ timer_info_receiver[i].second = rank; }
            else { timer_info_receiver[i].second = partner1; }
          }
          timer_info_receiver[i].first = std::max(timer_info_sender[i].first, timer_info_receiver[i].first);
        }
      }
    }
  }
  // Exchange the tracked routine critical path data
  if (partner1 == -1){
    MPI_Op op; MPI_Op_create((MPI_User_function*) propagate_critical_path_op,0,&op);
    if (tag < 20){
      PMPI_Allreduce(MPI_IN_PLACE, &critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, op, cm);
      MPI_Op_free(&op);
      if (mode >= 2) { propagate_timers(rank,tag,cm,partner1,partner2==-1 ? partner1 : partner2); }
    } else{
      MPI_Request req1;
      double* local_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
      std::memcpy(local_path_data, &critical_path_costs[0], critical_path_costs.size()*sizeof(double));
      PMPI_Iallreduce(MPI_IN_PLACE,local_path_data,critical_path_costs.size(),MPI_DOUBLE,op,cm,&req1);
      //MPI_Op_free(&op);
      internal_comm_prop.push_back(std::make_pair(local_path_data,true));
      internal_comm_prop_req.push_back(req1);
      if (mode>=2) { assert(0); initiate_timers(rank,tag,cm,partner1,partner2==-1 ? partner1 : partner2); }
    }
  }
  else{
   if (tag < 18){
      if (rank != partner1){
        PMPI_Sendrecv(&critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, partner1, internal_tag2, &new_cs[0], critical_path_costs.size(), MPI_DOUBLE, partner1, internal_tag2, cm, MPI_STATUS_IGNORE);
        update_critical_path(&new_cs[0],&critical_path_costs[0],critical_path_costs_size);
        if ((mode>=2) && (rank != partner2)){ propagate_timers(rank,tag,cm,partner1,partner2==-1 ? partner1 : partner2); }
      }
    }
    else{// only possible for Isend and Irecv
      MPI_Request req1,req2;
      if ((is_sender) && (rank != partner1)){
        double* local_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
        std::memcpy(local_path_data, &critical_path_costs[0], critical_path_costs.size()*sizeof(double));
        //TODO: Can I keep sending out `critical_path_costs` or must I make copies and send that out?
        PMPI_Isend(local_path_data, critical_path_costs.size(), MPI_DOUBLE, partner1, internal_tag2, cm, &req1);
        double* remote_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
        PMPI_Irecv(remote_path_data, critical_path_costs.size(), MPI_DOUBLE, partner1, internal_tag2, cm, &req2);
        internal_comm_prop.push_back(std::make_pair(local_path_data,true));
        internal_comm_prop_req.push_back(req1);
        internal_comm_prop.push_back(std::make_pair(remote_path_data,false));
        internal_comm_prop_req.push_back(req2);
      }
      else if ((!is_sender) && (rank != partner1)){
        double* local_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
        std::memcpy(local_path_data, &critical_path_costs[0], critical_path_costs.size()*sizeof(double));
        double* remote_path_data = (double*)malloc(critical_path_costs.size()*sizeof(double));
        PMPI_Irecv(remote_path_data, critical_path_costs.size(), MPI_DOUBLE, partner1, internal_tag2, cm, &req2);
        PMPI_Isend(local_path_data, critical_path_costs.size(), MPI_DOUBLE, partner1, internal_tag2, cm, &req1);
        internal_comm_prop.push_back(std::make_pair(local_path_data,true));
        internal_comm_prop_req.push_back(req1);
        internal_comm_prop.push_back(std::make_pair(remote_path_data,false));
        internal_comm_prop_req.push_back(req2);
      }
      if ((mode>=2) && (rank != partner1) && (rank != partner2)){ initiate_timers(rank,tag,cm,partner1,partner2==-1 ? partner1 : partner2); }
    }
  }
}

// Note: this function should be called once per start/stop, else it will double count
void find_per_process_max(MPI_Comm cm){
  int rank; MPI_Comm_rank(cm,&rank);
  double_int buffer[num_per_process_measures];
  for (size_t i=0; i<num_per_process_measures; i++){
    max_per_process_costs[i] = volume_costs[i];
    buffer[i].first          = volume_costs[i];
    buffer[i].second         = rank;
  }
  PMPI_Allreduce(MPI_IN_PLACE, &max_per_process_costs[0], num_per_process_measures, MPI_DOUBLE, MPI_MAX, cm);
  PMPI_Allreduce(MPI_IN_PLACE, &buffer[0], num_per_process_measures, MPI_DOUBLE_INT, MPI_MAXLOC, cm);
  size_t save=0;
  for (size_t i=0; i<num_per_process_measures-1; i++){// don't consider idle time an option
    if (breakdown[i] == 0) continue;
    size_t z = i<(2*cost_model_size) ? i : i+1;	// careful indexing to avoid idle time
    if (rank == buffer[z].second){
      for (size_t j=0; j<num_tracker_per_process_measures*list_size; j++){
        max_per_process_costs[num_per_process_measures+save*(num_tracker_per_process_measures*list_size+2)+j] = volume_costs[num_volume_measures+j];
      }
      max_per_process_costs[num_per_process_measures+(save+1)*(num_tracker_per_process_measures*list_size+2)-2] = volume_costs[num_volume_measures-2];
      max_per_process_costs[num_per_process_measures+(save+1)*(num_tracker_per_process_measures*list_size+2)-1] = volume_costs[num_volume_measures-6];
    }
    else{
      for (size_t j=0; j<num_tracker_per_process_measures*list_size+2; j++){
        max_per_process_costs[num_per_process_measures+save*(num_tracker_per_process_measures*list_size+2)+j] = 0.;
      }
    }
    PMPI_Allreduce(MPI_IN_PLACE, &max_per_process_costs[num_per_process_measures+save*(num_tracker_per_process_measures*list_size+2)], num_tracker_per_process_measures*list_size+2, MPI_DOUBLE, MPI_MAX, cm);
    save++;
  }
  // For now, buffer[num_per_process_measures-1].second holds the rank of the process with the max per-process runtime
  if (mode>=2){
    int per_process_runtime_root_rank = buffer[num_per_process_measures-1].second;
    // We consider only critical path runtime
    int ftimer_size = 0;
    if (rank==per_process_runtime_root_rank){
      ftimer_size = symbol_timers.size();
    }
    PMPI_Allreduce(MPI_IN_PLACE,&ftimer_size,1,MPI_INT,MPI_SUM,cm);

    symbol_len_pad.fill(0);
    if (rank==per_process_runtime_root_rank){
      int symbol_offset = 0;
      for (auto i=0; i<symbol_timers.size(); i++){
        symbol_len_pad[i] = symbol_order[i].size();
        for (auto j=0; j<symbol_len_pad[i]; j++){
          symbol_pad[symbol_offset+j] = symbol_order[i][j];
        }
        symbol_offset += symbol_len_pad[i];
      }
    }
    PMPI_Allreduce(MPI_IN_PLACE,&symbol_len_pad[0],ftimer_size,MPI_INT,MPI_SUM,cm);

    int num_chars = 0;
    for (auto i=0; i<ftimer_size; i++){
      num_chars += symbol_len_pad[i];
    }
    //TODO: What about the case for PMPI_Sendrecv (so when partner2 != -1)?
    if (rank == per_process_runtime_root_rank){
      PMPI_Bcast(&symbol_timer_pad_local_pp[0],(num_ftimer_measures*num_per_process_measures+1)*ftimer_size,MPI_DOUBLE,rank,cm);
      PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,rank,cm);
    }
    else{
      PMPI_Bcast(&symbol_timer_pad_global_pp[0],(num_ftimer_measures*num_per_process_measures+1)*ftimer_size,MPI_DOUBLE,per_process_runtime_root_rank,cm);
      PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,per_process_runtime_root_rank,cm);
      int symbol_offset = 0;
      for (int i=0; i<ftimer_size; i++){
        auto reconstructed_symbol = std::string(symbol_pad.begin()+symbol_offset,symbol_pad.begin()+symbol_offset+symbol_len_pad[i]);
        if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
          symbol_timers[reconstructed_symbol] = ftimer(reconstructed_symbol);
          symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
        }
        *symbol_timers[reconstructed_symbol].pp_numcalls = symbol_timer_pad_global_pp[(num_ftimer_measures*num_per_process_measures+1)*i];
        for (int j=0; j<num_per_process_measures; j++){
          *symbol_timers[reconstructed_symbol].pp_incl_measure[j] = symbol_timer_pad_global_pp[(num_ftimer_measures*num_per_process_measures+1)*i+2*j+1];
          *symbol_timers[reconstructed_symbol].pp_excl_measure[j] = symbol_timer_pad_global_pp[(num_ftimer_measures*num_per_process_measures+1)*i+2*(j+1)];
        }
        symbol_timers[reconstructed_symbol].has_been_processed = true;
        symbol_offset += symbol_len_pad[i];
      }
      // Now cycle through and find the symbols that were not processed and set their accumulated measures to 0
      for (auto& it : symbol_timers){
        if (it.second.has_been_processed){ it.second.has_been_processed = false; }
        else{
          *it.second.pp_numcalls = 0;
          for (int j=0; j<num_per_process_measures; j++){
            *it.second.pp_incl_measure[j] = 0;
            *it.second.pp_excl_measure[j] = 0;
          }
        }
      }
    }
  }
}

void compute_volume(MPI_Comm cm){
  if (mode>=1){
    PMPI_Allreduce(MPI_IN_PLACE, &volume_costs[0], volume_costs.size(), MPI_DOUBLE, MPI_SUM, cm);
    int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    for (int i=0; i<volume_costs.size(); i++){ volume_costs[i] /= (1.*world_size); }
  }
  if (mode>=2){ print_volume_symbol=false; }
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
    std::vector<double> vec(num_tracker_per_process_measures);
    int save=0;
    for (int j=0; j<cost_models.size(); j++){
      if (cost_models[j]){
        vec[2*save] = max_per_process_costs[num_per_process_measures+idx*(num_tracker_per_process_measures*list_size+2)+this->tag*num_tracker_per_process_measures+save];
        vec[2*save+1] = max_per_process_costs[num_per_process_measures+idx*(num_tracker_per_process_measures*list_size+2)+this->tag*num_tracker_per_process_measures+cost_model_size+save];
        save++;
      }
    }
    // For now, do not include idle time
    vec[num_tracker_per_process_measures-3] = max_per_process_costs[num_per_process_measures+idx*(num_tracker_per_process_measures*list_size+2)+this->tag*num_tracker_per_process_measures+num_tracker_per_process_measures-3];
    vec[num_tracker_per_process_measures-2] = max_per_process_costs[num_per_process_measures+idx*(num_tracker_per_process_measures*list_size+2)+this->tag*num_tracker_per_process_measures+num_tracker_per_process_measures-2];
    vec[num_tracker_per_process_measures-1] = max_per_process_costs[num_per_process_measures+idx*(num_tracker_per_process_measures*list_size+2)+this->tag*num_tracker_per_process_measures+num_tracker_per_process_measures-1];
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
    vec[2*cost_model_size] = *this->my_comm_time;
    vec[2*cost_model_size+1] = *this->my_synch_time;
    vec[2*cost_model_size+2] = *this->my_datamvt_time;
    save_info[this->name] = std::move(vec);
  }
}

ftimer::ftimer(std::string name_){
  assert(name_.size() <= max_timer_name_length);
  assert(symbol_timers.size() < max_num_symbols);
  this->name = std::move(name_);
  this->cp_numcalls = &symbol_timer_pad_local_cp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)]; *this->cp_numcalls = 0;
  this->pp_numcalls = &symbol_timer_pad_local_pp[(symbol_timers.size()-1)*(num_ftimer_measures*num_per_process_measures+1)]; *this->pp_numcalls = 0;
  this->vol_numcalls = &symbol_timer_pad_vol[(symbol_timers.size()-1)*(num_ftimer_measures*num_volume_measures+1)]; *this->vol_numcalls = 0;
  for (auto i=0; i<num_critical_path_measures; i++){
    this->cp_incl_measure[i] = &symbol_timer_pad_local_cp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)+2*i+1]; *cp_incl_measure[i] = 0.;
    this->cp_excl_measure[i] = &symbol_timer_pad_local_cp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)+2*(i+1)]; *cp_excl_measure[i] = 0.;
  }
  for (auto i=0; i<num_per_process_measures; i++){
    this->pp_incl_measure[i] = &symbol_timer_pad_local_pp[(symbol_timers.size()-1)*(num_ftimer_measures*num_per_process_measures+1)+2*i+1]; *pp_incl_measure[i] = 0.;
    this->pp_excl_measure[i] = &symbol_timer_pad_local_pp[(symbol_timers.size()-1)*(num_ftimer_measures*num_per_process_measures+1)+2*(i+1)]; *pp_excl_measure[i] = 0.;
  }
  for (auto i=0; i<num_volume_measures; i++){
    this->vol_incl_measure[i] = &symbol_timer_pad_vol[(symbol_timers.size()-1)*(num_ftimer_measures*num_volume_measures+1)+2*i+1]; *vol_incl_measure[i] = 0.;
    this->vol_excl_measure[i] = &symbol_timer_pad_vol[(symbol_timers.size()-1)*(num_ftimer_measures*num_volume_measures+1)+2*(i+1)]; *vol_excl_measure[i] = 0.;
  }
  this->has_been_processed = false;
  this->cp_exclusive_contributions.fill(0.0); this->pp_exclusive_contributions.fill(0.0);
  this->cp_exclusive_measure.fill(0.0); this->pp_exclusive_measure.fill(0.0);
}

void ftimer::start(double save_time){
  if (symbol_stack.size()>0){
    auto last_symbol_time = save_time-symbol_timers[symbol_stack.top()].start_timer.top();
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-2] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-2] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-2] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-1] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-2] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += last_symbol_time;
  }
  critical_path_costs[num_critical_path_measures-2] += (save_time - computation_timer);		// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += (save_time - computation_timer);		// update critical path runtime
  volume_costs[num_volume_measures-2]        += (save_time - computation_timer);		// update local computation time
  volume_costs[num_volume_measures-1]        += (save_time - computation_timer);		// update local runtime
  for (size_t i=0; i<breakdown_size; i++){
    critical_path_costs[critical_path_costs_size-1-i] += (save_time - computation_timer);;
  }
  symbol_stack.push(this->name);
  std::array<double,num_per_process_measures> temp2; temp2.fill(0.0);
  computation_timer = MPI_Wtime();
  this->start_timer.push(computation_timer);
}

void ftimer::stop(double save_time){
  assert(this->start_timer.size()>0);
  auto last_symbol_time = save_time-this->start_timer.top();
  this->cp_exclusive_measure[num_critical_path_measures-1] += last_symbol_time;
  this->cp_exclusive_measure[num_critical_path_measures-2] += last_symbol_time;
  this->pp_exclusive_measure[num_per_process_measures-1] += last_symbol_time;
  this->pp_exclusive_measure[num_per_process_measures-2] += last_symbol_time;
  *this->cp_excl_measure[num_critical_path_measures-1] += last_symbol_time;
  *this->cp_excl_measure[num_critical_path_measures-2] += last_symbol_time;
  *this->pp_excl_measure[num_per_process_measures-1] += last_symbol_time;
  *this->pp_excl_measure[num_per_process_measures-2] += last_symbol_time;
  *this->cp_numcalls += 1.; *this->pp_numcalls += 1.; *this->vol_numcalls += 1.;

  for (auto i=0; i<num_critical_path_measures; i++){ this->cp_exclusive_contributions[i] += this->cp_exclusive_measure[i]; } this->cp_exclusive_measure.fill(0.0);
  for (auto i=0; i<num_per_process_measures; i++){ this->pp_exclusive_contributions[i] += this->pp_exclusive_measure[i]; } this->pp_exclusive_measure.fill(0.0);

  auto save_symbol = symbol_stack.top();
  this->start_timer.pop(); symbol_stack.pop();
  if (symbol_stack.size() > 0 && (save_symbol != symbol_stack.top())){
    for (auto i=0; i<num_critical_path_measures; i++){ symbol_timers[symbol_stack.top()].cp_exclusive_contributions[i] += this->cp_exclusive_contributions[i];
                                                       *this->cp_incl_measure[i] += this->cp_exclusive_contributions[i]; }
    for (auto i=0; i<num_per_process_measures; i++){ symbol_timers[symbol_stack.top()].pp_exclusive_contributions[i] += this->pp_exclusive_contributions[i];
                                                     *this->pp_incl_measure[i] += this->pp_exclusive_contributions[i]; }
    this->cp_exclusive_contributions.fill(0.0); this->pp_exclusive_contributions.fill(0.0);
  } else if (symbol_stack.size() == 0){
    for (auto i=0; i<num_critical_path_measures; i++){ *this->cp_incl_measure[i] += this->cp_exclusive_contributions[i]; }
    for (auto i=0; i<num_per_process_measures; i++){ *this->pp_incl_measure[i] += this->pp_exclusive_contributions[i]; }
    this->cp_exclusive_contributions.fill(0.0); this->pp_exclusive_contributions.fill(0.0);
  }

  critical_path_costs[num_critical_path_measures-2] += (save_time - computation_timer);		// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += (save_time - computation_timer);		// update critical path runtime
  volume_costs[num_volume_measures-2]        += (save_time - computation_timer);		// update local computation time
  volume_costs[num_volume_measures-1]        += (save_time - computation_timer);		// update local runtime
  for (size_t i=0; i<breakdown_size; i++){ critical_path_costs[critical_path_costs_size-1-i] += (save_time - computation_timer); }
  computation_timer = MPI_Wtime();
  if (symbol_stack.size()>0){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
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

std::string get_measure_title(size_t idx){
  std::array<std::string,per_process_costs_size> measure_titles =\
    {"est_comm_bsp","est_comm_ab","est_synch_bsp","est_synch_ab","idle time","comm time","synch time","datamvt time","comp time","runtime"};
  if ((cost_models[0]==0) && (cost_models[1]==0)){ return measure_titles[idx+2*cost_model_size]; }
  if ((cost_models[0]==0) && (cost_models[1]==1)){ if (idx==0) return measure_titles[0]; if (idx==1) return measure_titles[2]; return measure_titles[2+idx];}
  if ((cost_models[0]==1) && (cost_models[1]==0)){ if (idx==0) return measure_titles[1]; if (idx==1) return measure_titles[3]; return measure_titles[2+idx];}
  if ((cost_models[0]==1) && (cost_models[1]==1)){ return measure_titles[idx]; }
}

template<typename StreamType>
void print_cost_model_header(StreamType& Stream){
  if (cost_models[0]){
    Stream << std::left << std::setw(mode_1_width) << "BSPCommCost";
    Stream << std::left << std::setw(mode_1_width) << "BSPSynchCost";
  }
  if (cost_models[1]){
    Stream << std::left << std::setw(mode_1_width) << "ABbutterflyCommCost";
    Stream << std::left << std::setw(mode_1_width) << "ABbutterflySynchCost";
  }
}

template<typename StreamType>
void print_cost_model_header_file(StreamType& Stream){
  if (cost_models[0]){
    Stream << "\tBSPCommCost";
    Stream << "\tBSPSynchCost";
  }
  if (cost_models[1]){
    Stream << "\tABbutterflyCommCost";
    Stream << "\tABbutterflySynchCost";
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
  Stream << "\tIdleTime\tCommunicationTime\tSynchronizationTime\tDataMvtTime\tComputationTime\tRunTime";// per-process
  print_cost_model_header_file(Stream);
  Stream << "\tIdleTime\tCommunicationTime\tSynchronizationTime\tDataMvtTime\tComputationTime\tRunTime";// volume
  for (auto i=0; i<num_tracker_critical_path_measures*breakdown_size+num_tracker_per_process_measures*breakdown_size+num_tracker_volume_measures;i++){
    for (auto& it : save_info){
     Stream << "\t" << it.first;
    }
  }
}

void record(std::ofstream& Stream, size_t factor){
  assert(internal_comm_info.size() == 0);
  if (mode>=0){
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
      for (size_t i=0; i<num_per_process_measures; i++){
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
      for (auto i=0; i<num_per_process_measures-1; i++){	// no idle time
        if (!breakdown[i]) continue;
        // Save the critter information before printing
        for (size_t j=0; j<list_size; j++){
          list[j]->set_per_process_costs(breakdown_idx);
        }
        for (size_t j=0; j<num_tracker_per_process_measures; j++){
          for (auto& it : save_info){
            Stream << "\t" << factor*it.second[j];
          }
        }
        Stream << "\t" << factor*max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+2)-2];
        Stream << "\t" << factor*max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+2)-1];
        breakdown_idx++;
      }
      breakdown_idx=0;
      for (auto i=0; i<num_critical_path_measures; i++){
        if (!breakdown[i]) continue;
        Stream << "\t" << factor*critical_path_costs[critical_path_costs_size-breakdown_size+breakdown_idx];
        breakdown_idx++;
      }
      breakdown_idx=0;
      for (auto i=0; i<num_critical_path_measures; i++){
        if (!breakdown[i]) continue;
        // Save the critter information before printing
        for (size_t j=0; j<list_size; j++){
          list[j]->set_critical_path_costs(breakdown_idx);
        }
        for (size_t j=0; j<num_tracker_critical_path_measures; j++){
          for (auto& it : save_info){
            Stream << "\t" << factor*it.second[j];
          }
        }
        breakdown_idx++;
      }
    }
  }
}

void record(std::ostream& Stream, size_t factor){
  assert(internal_comm_info.size() == 0);
  if (mode==0){
    if (is_world_root){
      Stream << std::left << std::setw(mode_1_width) << "Runtime:";
      Stream << "\n";
      Stream << std::left << std::setw(mode_1_width) << "                  ";
      Stream << std::left << std::setw(mode_1_width) << factor*critical_path_costs[num_critical_path_measures-1];
      Stream << "\n\n";
    }
  }
  if (mode>=1){
    if (is_world_root){
      Stream << "\n\n";
      Stream << std::left << std::setw(mode_1_width) << "Critical path:";
      print_cost_model_header(Stream);
      Stream << std::left << std::setw(mode_1_width) << "";
      Stream << std::left << std::setw(mode_1_width) << "CommTime";
      Stream << std::left << std::setw(mode_1_width) << "SynchTime";
      Stream << std::left << std::setw(mode_1_width) << "DataMvtTime";
      Stream << std::left << std::setw(mode_1_width) << "CompTime";
      Stream << std::left << std::setw(mode_1_width) << "RunTime";
      Stream << "\n";
      Stream << std::left << std::setw(mode_1_width) << "                  ";
      for (size_t i=0; i<num_critical_path_measures+1; i++){//+1 for blank space
        if (i==(2*cost_model_size)) Stream << std::left << std::setw(mode_1_width) << "";
        else if ((i<(2*cost_model_size)) && (i%2==0)) Stream << std::left << std::setw(mode_1_width) << factor*critical_path_costs[i/2];
        else if ((i<(2*cost_model_size)) && (i%2==1)) Stream << std::left << std::setw(mode_1_width) << factor*critical_path_costs[(i-1)/2+cost_model_size];
        else Stream << std::left << std::setw(mode_1_width) << factor*critical_path_costs[i-1];
      }
      Stream << "\n\n";

      Stream << std::left << std::setw(mode_1_width) << "Per-process max:";
      print_cost_model_header(Stream);
      Stream << std::left << std::setw(mode_1_width) << "IdleTime";
      Stream << std::left << std::setw(mode_1_width) << "CommTime";
      Stream << std::left << std::setw(mode_1_width) << "SynchTime";
      Stream << std::left << std::setw(mode_1_width) << "DataMvtTime";
      Stream << std::left << std::setw(mode_1_width) << "CompTime";
      Stream << std::left << std::setw(mode_1_width) << "RunTime";
      Stream << "\n";
      Stream << std::left << std::setw(mode_1_width) << "                  ";
      for (size_t i=0; i<num_volume_measures; i++){
        if ((i<(2*cost_model_size)) && (i%2==0)) Stream << std::left << std::setw(mode_1_width) << factor*max_per_process_costs[i/2];
        else if ((i<(2*cost_model_size)) && (i%2==1)) Stream << std::left << std::setw(mode_1_width) << factor*max_per_process_costs[(i-1)/2+cost_model_size];
        else Stream << std::left << std::setw(mode_1_width) << factor*max_per_process_costs[i];
      }
      Stream << "\n\n";

      Stream << std::left << std::setw(mode_1_width) << "Volume:";
      print_cost_model_header(Stream);
      Stream << std::left << std::setw(mode_1_width) << "IdleTime";
      Stream << std::left << std::setw(mode_1_width) << "CommTime";
      Stream << std::left << std::setw(mode_1_width) << "SynchTime";
      Stream << std::left << std::setw(mode_1_width) << "DataMvtTime";
      Stream << std::left << std::setw(mode_1_width) << "CompTime";
      Stream << std::left << std::setw(mode_1_width) << "RunTime";
      Stream << "\n";
      Stream << std::left << std::setw(mode_1_width) << "                  ";
      for (size_t i=0; i<num_volume_measures; i++){
        if ((i<(2*cost_model_size)) && (i%2==0)) Stream << std::left << std::setw(mode_1_width) << factor*volume_costs[i/2];
        else if ((i<(2*cost_model_size)) && (i%2==1)) Stream << std::left << std::setw(mode_1_width) << factor*volume_costs[(i-1)/2+cost_model_size];
        else Stream << std::left << std::setw(mode_1_width) << factor*volume_costs[i];
      }
      Stream << "\n\n";

      size_t breakdown_idx=0;
      for (auto i=0; i<num_critical_path_measures; i++){
        if (!breakdown[i]) continue;
        if (i==0){
          Stream << std::left << std::setw(mode_1_width) << "BSPCommCost max:";
        } else if (i==1){
          Stream << std::left << std::setw(mode_1_width) << "ABCommCost max:";
        } else if (i==2){
          Stream << std::left << std::setw(mode_1_width) << "BSPSynchCost max:";
        } else if (i==3){
          Stream << std::left << std::setw(mode_1_width) << "ABSynchCost max:";
        } else if (i==4){
          Stream << std::left << std::setw(mode_1_width) << "CommTime max:";
        } else if (i==5){
          Stream << std::left << std::setw(mode_1_width) << "SynchTime max:";
        } else if (i==6){
          Stream << std::left << std::setw(mode_1_width) << "DataMvtTime max:";
        } else if (i==7){
          Stream << std::left << std::setw(mode_1_width) << "CompTime max:";
        } else if (i==8){
          Stream << std::left << std::setw(mode_1_width) << "RunTime max:";
        }
        Stream << std::left << std::setw(mode_1_width) << "MeasureType";
        Stream << std::left << std::setw(mode_1_width) << "CompTime";
        Stream << std::left << std::setw(mode_1_width) << "IdleTime";
        print_cost_model_header(Stream);
        Stream << std::left << std::setw(mode_1_width) << "CommTime";
        Stream << std::left << std::setw(mode_1_width) << "SynchTime";
        Stream << std::left << std::setw(mode_1_width) << "DataMvtTime";
        Stream << "\n";
        Stream << std::left << std::setw(mode_1_width) << "Computation";
        Stream << std::left << std::setw(mode_1_width) << "path";
        Stream << std::left << std::setw(mode_1_width) << factor*critical_path_costs[critical_path_costs_size-breakdown_size+breakdown_idx];
        Stream << "\n";
        Stream << std::left << std::setw(mode_1_width) << "Idle";
        Stream << std::left << std::setw(mode_1_width) << "path";
        Stream << std::left << std::setw(mode_1_width) << 0.0;
        Stream << std::left << std::setw(mode_1_width) << 0.0;
        for (int j=0; j<list_size; j++){
          list[j]->set_critical_path_costs(breakdown_idx);
        }
        for (auto& it : save_info){
          Stream << "\n";
          Stream << std::left << std::setw(mode_1_width) << it.first;
          Stream << std::left << std::setw(mode_1_width) << "path";
          Stream << std::left << std::setw(mode_1_width) << 0.0;
          Stream << std::left << std::setw(mode_1_width) << 0.0;
          for (size_t j=0; j<num_tracker_critical_path_measures; j++){
            Stream << std::left << std::setw(mode_1_width) << factor*it.second[j];
          }
        }
        Stream << "\n";
        Stream << std::left << std::setw(mode_1_width) << "Computation";
        Stream << std::left << std::setw(mode_1_width) << "per-process";
        Stream << std::left << std::setw(mode_1_width) << factor*max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+2)-2];
        Stream << "\n";
        Stream << std::left << std::setw(mode_1_width) << "Idle";
        Stream << std::left << std::setw(mode_1_width) << "per-process";
        Stream << std::left << std::setw(mode_1_width) << 0.0;
        Stream << std::left << std::setw(mode_1_width) << factor*max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+2)-1];
        for (int j=0; j<list_size; j++){
          list[j]->set_per_process_costs(breakdown_idx);
        }
        for (auto& it : save_info){
          Stream << "\n";
          Stream << std::left << std::setw(mode_1_width) << it.first;
          Stream << std::left << std::setw(mode_1_width) << "per-process";
          Stream << std::left << std::setw(mode_1_width) << 0.0;
          Stream << std::left << std::setw(mode_1_width) << 0.0;
          for (size_t j=0; j<num_tracker_per_process_measures; j++){
            Stream << std::left << std::setw(mode_1_width) << factor*it.second[j];
          }
        }
        breakdown_idx++;
        Stream << "\n\n";
      }
      for (int i=0; i<list_size; i++){
        list[i]->set_volume_costs();
      }
      Stream << std::left << std::setw(mode_1_width) << "Volume:";
      print_cost_model_header(Stream);
      Stream << std::left << std::setw(mode_1_width) << "CommTime";
      Stream << std::left << std::setw(mode_1_width) << "SynchTime";
      Stream << std::left << std::setw(mode_1_width) << "DataMvtTime";
      for (auto& it : save_info){
        Stream << "\n";
        Stream << std::left << std::setw(mode_1_width) << it.first;
        for (size_t j=0; j<num_tracker_volume_measures; j++){
          Stream << std::left << std::setw(mode_1_width) << factor*it.second[j];
        }
      }
      Stream << "\n";
    }
  }
  if (mode>=2){
    if (is_world_root){
      Stream << "***********************************************************************************************************************";
      std::vector<std::pair<std::string,std::array<double,6>>> sort_info(symbol_timers.size());
      for (int i=num_per_process_measures-1; i>=0; i--){
        sort_info.clear(); sort_info.resize(symbol_timers.size());
        // Reset symbol timers and sort
        size_t j=0;
        double cp_ref,pp_ref,vol_ref;
        for (auto& it : symbol_timers){
          assert(it.second.start_timer.size() == 0);
          if (i==2*cost_model_size){
            sort_info[j++] = std::make_pair(it.second.name,std::array<double,6>{*it.second.cp_numcalls,0.,*it.second.pp_numcalls,*it.second.pp_excl_measure[i],*it.second.vol_numcalls,*it.second.vol_excl_measure[i]});
          } else if (i>2*cost_model_size){
            sort_info[j++] = std::make_pair(it.second.name,std::array<double,6>{*it.second.cp_numcalls,*it.second.cp_excl_measure[i-1],*it.second.pp_numcalls,*it.second.pp_excl_measure[i],*it.second.vol_numcalls,*it.second.vol_excl_measure[i]});
          } else{
            sort_info[j++] = std::make_pair(it.second.name,std::array<double,6>{*it.second.cp_numcalls,*it.second.cp_excl_measure[i],*it.second.pp_numcalls,*it.second.pp_excl_measure[i],*it.second.vol_numcalls,*it.second.vol_excl_measure[i]});
          }
        }
        std::sort(sort_info.begin(),sort_info.end(),[](std::pair<std::string,std::array<double,6>>& vec1, std::pair<std::string,std::array<double,6>>& vec2){return vec1.second[1] > vec2.second[1];});
        if (i==2*cost_model_size){
          cp_ref = 1.;
        } else if (i>2*cost_model_size){
          cp_ref = critical_path_costs[i-1];
        } else{
          cp_ref = critical_path_costs[i];
        }
        pp_ref = max_per_process_costs[i];
        vol_ref = volume_costs[i];
        if (cp_ref==0.) cp_ref=1.; if (pp_ref==0.) pp_ref=1.; if (vol_ref==0.) vol_ref=1.;

        // Exclusive
        Stream << "\n\n\n\n" << std::left << std::setw(max_timer_name_length) << get_measure_title(i);
        Stream << std::left << std::setw(mode_2_width) << "cp-#calls";
        Stream << std::left << std::setw(mode_2_width) << "pp-#calls";
        if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << "vol-#calls";
        Stream << std::left << std::setw(mode_2_width) << "cp-excl (s)";
        Stream << std::left << std::setw(mode_2_width) << "pp-excl (s)";
        if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << "vol-excl (s)";
        Stream << std::left << std::setw(mode_2_width) << "cp-excl (%)";
        Stream << std::left << std::setw(mode_2_width) << "pp-excl (%)";
        if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << "vol-excl (%)";
        double cp_total_exclusive = 0.;
        double pp_total_exclusive = 0.;
        double vol_total_exclusive = 0.;
        for (auto& it : sort_info){
          Stream << "\n" << std::left << std::setw(max_timer_name_length) << it.first;
          Stream << std::left << std::setw(mode_2_width) << it.second[0];
          Stream << std::left << std::setw(mode_2_width) << it.second[2];
          if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << it.second[4];
          Stream << std::left << std::setw(mode_2_width) << it.second[1];
          Stream << std::left << std::setw(mode_2_width) << it.second[3];
          if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << it.second[5];
          Stream << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[1]/cp_ref;
          Stream << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[3]/pp_ref;
          if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[5]/vol_ref;
          cp_total_exclusive += it.second[1];
          pp_total_exclusive += it.second[3];
          vol_total_exclusive += it.second[5];
        }
        Stream << "\n" << std::left << std::setw(max_timer_name_length) << "total";
        Stream << std::left << std::setw(mode_2_width) << "";
        Stream << std::left << std::setw(mode_2_width) << "";
        if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << "";
        Stream << std::left << std::setw(mode_2_width) << cp_total_exclusive;
        Stream << std::left << std::setw(mode_2_width) << pp_total_exclusive;
        if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << vol_total_exclusive;
        Stream << std::left << std::setw(mode_2_width) << 100.*cp_total_exclusive/cp_ref;
        Stream << std::left << std::setw(mode_2_width) << 100.*pp_total_exclusive/pp_ref;
        if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << 100.*vol_total_exclusive/vol_ref;
        Stream << "\n";

        // Reset symbol timers and sort
        sort_info.clear(); sort_info.resize(symbol_timers.size());
        j=0;
        for (auto& it : symbol_timers){
          assert(it.second.start_timer.size() == 0);
          if (i==2*cost_model_size){
            sort_info[j++] = std::make_pair(it.second.name,std::array<double,6>{*it.second.cp_numcalls,0.,*it.second.pp_numcalls,*it.second.pp_incl_measure[i],*it.second.vol_numcalls,*it.second.vol_incl_measure[i]});
          } else if (i>2*cost_model_size){
            sort_info[j++] = std::make_pair(it.second.name,std::array<double,6>{*it.second.cp_numcalls,*it.second.cp_incl_measure[i-1],*it.second.pp_numcalls,*it.second.pp_incl_measure[i],*it.second.vol_numcalls,*it.second.vol_incl_measure[i]});
          } else{
            sort_info[j++] = std::make_pair(it.second.name,std::array<double,6>{*it.second.cp_numcalls,*it.second.cp_incl_measure[i],*it.second.pp_numcalls,*it.second.pp_incl_measure[i],*it.second.vol_numcalls,*it.second.vol_incl_measure[i]});
          }
        }
        std::sort(sort_info.begin(),sort_info.end(),[](std::pair<std::string,std::array<double,6>>& vec1, std::pair<std::string,std::array<double,6>>& vec2){return vec1.second[1] > vec2.second[1];});
        if (i==2*cost_model_size){
          cp_ref = 100.;
        } else if (i>2*cost_model_size){
          cp_ref = critical_path_costs[i-1];
        } else{
          cp_ref = critical_path_costs[i];
        }
        pp_ref = max_per_process_costs[i];
        vol_ref = volume_costs[i];
        if (cp_ref==0.) cp_ref=1.; if (pp_ref==0.) pp_ref=1.; if (vol_ref==0.) vol_ref=1.;

        // Inclusive
        Stream << "\n" << std::left << std::setw(max_timer_name_length) << get_measure_title(i);
        Stream << std::left << std::setw(mode_2_width) << "cp-#calls";
        Stream << std::left << std::setw(mode_2_width) << "pp-#calls";
        if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << "vol-#calls";
        Stream << std::left << std::setw(mode_2_width) << "cp-incl (s)";
        Stream << std::left << std::setw(mode_2_width) << "pp-incl (s)";
        if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << "vol-incl (s)";
        Stream << std::left << std::setw(mode_2_width) << "cp-incl (%)";
        Stream << std::left << std::setw(mode_2_width) << "pp-incl (%)";
        if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << "vol-incl (%)";
        double cp_total_inclusive = 0.;
        double pp_total_inclusive = 0.;
        double vol_total_inclusive = 0.;
        for (auto& it : sort_info){
          Stream << "\n" << std::left << std::setw(max_timer_name_length) << it.first;
          Stream << std::left << std::setw(mode_2_width) << it.second[0];
          Stream << std::left << std::setw(mode_2_width) << it.second[2];
          if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << it.second[4];
          Stream << std::left << std::setw(mode_2_width) << it.second[1];
          Stream << std::left << std::setw(mode_2_width) << it.second[3];
          if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << it.second[5];
          Stream << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[1]/cp_ref;
          Stream << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[3]/pp_ref;
          if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[5]/vol_ref;
          cp_total_inclusive = std::max(it.second[1],cp_total_inclusive);
          pp_total_inclusive = std::max(it.second[3],pp_total_inclusive);
          vol_total_inclusive = std::max(it.second[5],vol_total_inclusive);
        }
        Stream << "\n" << std::left << std::setw(max_timer_name_length) << "total";
        Stream << std::left << std::setw(mode_2_width) << "";
        Stream << std::left << std::setw(mode_2_width) << "";
        if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << "";
        Stream << std::left << std::setw(mode_2_width) << cp_total_inclusive;
        Stream << std::left << std::setw(mode_2_width) << pp_total_inclusive;
        if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << vol_total_inclusive;
        Stream << std::left << std::setw(mode_2_width) << 100.*cp_total_inclusive/cp_ref;
        Stream << std::left << std::setw(mode_2_width) << 100.*pp_total_inclusive/pp_ref;
        if (print_volume_symbol) Stream << std::left << std::setw(mode_2_width) << 100.*vol_total_inclusive/vol_ref;
        Stream << "\n";
      }
    }
  }
}
};

void start(size_t mode){
  assert(mode>=0 && mode < 3); assert(internal::internal_comm_info.size() == 0);
  internal::wait_id=true; internal::waitall_id=false; internal::print_volume_symbol=true; internal::mode=mode;
  // TODO: How to allow different number of cost models. Perhaps just put an assert that both cost models must be on? Or don't use these altogether?
  for (int i=0; i<internal::list_size; i++){ internal::list[i]->init(); }
  if (internal::is_world_root){
    if (!internal::is_first_iter){
      if (internal::flag) {internal::stream << "\n";} else {std::cout << "\n";}
    }
  }
  for (auto i=0; i<internal::critical_path_costs.size(); i++){ internal::critical_path_costs[i]=0.; }
  for (auto i=0; i<internal::volume_costs.size(); i++){ internal::volume_costs[i]=0.; }
  /*Initiate new timer*/
  internal::computation_timer=MPI_Wtime();
}

void stop(size_t mode, size_t factor){
  volatile double last_time = MPI_Wtime();
  assert(internal::internal_comm_info.size() == 0);
  internal::critical_path_costs[internal::num_critical_path_measures-2]+=(last_time-internal::computation_timer);	// update critical path computation time
  internal::critical_path_costs[internal::num_critical_path_measures-1]+=(last_time-internal::computation_timer);	// update critical path runtime
  internal::volume_costs[internal::num_volume_measures-2]+=(last_time-internal::computation_timer);			// update computation time volume
  internal::volume_costs[internal::num_volume_measures-1]+=(last_time-internal::computation_timer);			// update runtime volume
  for (size_t i=0; i<breakdown_size; i++){ internal::critical_path_costs[internal::critical_path_costs_size-1-i] += (last_time-internal::computation_timer); }
  internal::propagate(MPI_COMM_WORLD, 0, true, -1, -1);
  internal::find_per_process_max(MPI_COMM_WORLD);
  internal::compute_volume(MPI_COMM_WORLD);

  internal::record(std::cout,factor);
  if (internal::flag) {internal::record(internal::stream,factor);}

  internal::wait_id=false; internal::waitall_id=false; internal::is_first_iter = false;
  internal::mode=0; internal::save_info.clear();
  for (auto i=0; i<internal::critical_path_costs.size(); i++){ internal::critical_path_costs[i]=0.; }
  for (auto i=0; i<internal::max_per_process_costs.size(); i++){ internal::max_per_process_costs[i]=0.; }
  for (auto i=0; i<internal::volume_costs.size(); i++){ internal::volume_costs[i]=0.; }
  internal::need_new_line=false; internal::symbol_timers.clear();
}
};
