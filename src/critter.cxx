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

void propagate_critical_path_op(double* in, double* inout, int* len, MPI_Datatype* dtype){
  if (critical_path_breakdown_size > 0){
    size_t breakdown_idx=0;
    size_t breakdown_size = critical_path_breakdown_size;	// prevents compiler warning
    for (int i=0; i<num_critical_path_measures; i++){
      if (critical_path_breakdown[i]) decisions[breakdown_idx++] = inout[i] > in[i];
      inout[i] = std::max(inout[i],in[i]);
    }
    for (int i=num_critical_path_measures; i<*len; i++){
      int idx = (i-num_critical_path_measures)%breakdown_size;
      inout[i] = (decisions[idx] ? inout[i] : in[i]);
    }
  } else{
    for (int i=0; i<num_critical_path_measures; i++){
      inout[i] = std::max(inout[i],in[i]);
    }
  }
}

tracker _MPI_Barrier("MPI_Barrier",0, 
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),0.); 
                          }),
        _MPI_Bcast("MPI_Bcast",1,
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }),
        _MPI_Reduce("MPI_Reduce",2, 
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }),
        _MPI_Allreduce("MPI_Allreduce",3,
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }),
        _MPI_Gather("MPI_Gather",4,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }),
        _MPI_Gatherv("MPI_Gatherv",5,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }),
        _MPI_Allgather("MPI_Allgather",6,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }),
        _MPI_Allgatherv("MPI_Allgatherv",7,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }),
        _MPI_Scatter("MPI_Scatter",8,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }),
        _MPI_Scatterv("MPI_Scatterv",9,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }),
        _MPI_Reduce_scatter("MPI_Reduce_scatter",10,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }),
        _MPI_Alltoall("MPI_Alltoall",11,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n); 
                          }),
        _MPI_Alltoallv("MPI_Alltoallv",12,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n); 
                          }),
        _MPI_Send("MPI_Send",13,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }),
        _MPI_Recv("MPI_Recv",14,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }),
        _MPI_Isend("MPI_Isend",15,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }),
        _MPI_Irecv("MPI_Irecv",16,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }),
        _MPI_Sendrecv("MPI_Sendrecv",17,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }),
        _MPI_Sendrecv_replace("MPI_Sendrecv_replace",18,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }),
        _MPI_Ibcast("MPI_Ibcast",19,
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }),
        _MPI_Iallreduce("MPI_Iallreduce",20,
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }),
        _MPI_Ireduce("MPI_Ireduce",21,
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }),
        _MPI_Igather("MPI_Igather",22,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }),
        _MPI_Igatherv("MPI_Igatherv",23,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }),
        _MPI_Iallgather("MPI_Iallgather",24,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }),
        _MPI_Iallgatherv("MPI_Iallgatherv",25,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }),
        _MPI_Iscatter("MPI_Iscatter",26,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }),
        _MPI_Iscatterv("MPI_Iscatterv",27,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }),
        _MPI_Ireduce_scatter("MPI_Ireduce_scatter",28,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }),
        _MPI_Ialltoall("MPI_Ialltoall",29,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n); 
                          }),
        _MPI_Ialltoallv("MPI_Ialltoallv",30,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n); 
                          }),
        _MPI_Ssend("MPI_Ssend",31,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          });


tracker* list[list_size] = {
        &_MPI_Barrier,
        &_MPI_Bcast,
        &_MPI_Reduce,
        &_MPI_Allreduce,
        &_MPI_Gather,
        &_MPI_Gatherv,
        &_MPI_Allgather,
        &_MPI_Allgatherv,
        &_MPI_Scatter,
        &_MPI_Scatterv,
        &_MPI_Reduce_scatter,
        &_MPI_Alltoall,
        &_MPI_Alltoallv,
        &_MPI_Send,
        &_MPI_Recv,
        &_MPI_Irecv,
        &_MPI_Isend,
        &_MPI_Sendrecv,
        &_MPI_Sendrecv_replace,
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
        &_MPI_Ssend};

std::string stream_name,file_name;
std::ofstream stream;
bool flag,is_world_root,is_first_iter,need_new_line;
size_t mode;

double computation_timer;
std::map<MPI_Request,std::pair<MPI_Request,bool>> internal_comm_info;
std::map<MPI_Request,double*> internal_comm_message;
std::map<MPI_Request,std::pair<double,double>> internal_comm_data;
std::map<MPI_Request,tracker*> internal_comm_track;
bool decisions[critical_path_breakdown_size];
std::array<double,critical_path_costs_size> critical_path_costs;
std::array<double,volume_costs_size> volume_costs;
std::array<double,num_critical_path_measures> max_per_process_costs;
std::map<std::string,std::vector<double>> save_info;
double new_cs[critical_path_costs_size];
double scratch_pad;
std::vector<char> synch_pad_send;
std::vector<char> synch_pad_recv;
std::array<char,max_timer_name_length*max_num_symbols> symbol_pad;
std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_local;
std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_global;
std::unordered_map<std::string,ftimer> symbol_timers;
std::stack<std::string> symbol_stack;
std::array<std::string,max_num_symbols> symbol_order;
std::array<std::string,num_critical_path_measures> critical_path_measure_names;
double_int timer_cp_info_sender[num_critical_path_measures];
double_int timer_cp_info_receiver[num_critical_path_measures];

void tracker::init(){
  this->last_start_time  = -1.;
  this->save_comp_time   = 0.;
}

void tracker::set_cost_pointers(){
  size_t volume_costs_idx        = num_volume_measures+this->tag*num_tracker_volume_measures;
  this->my_bytes                 = &volume_costs[volume_costs_idx];
  this->my_msg                   = &volume_costs[volume_costs_idx+1];
  this->my_wrd                   = &volume_costs[volume_costs_idx+2];
  this->my_bar_time              = &volume_costs[volume_costs_idx+3];
  this->my_comm_time             = &volume_costs[volume_costs_idx+4];
  this->my_synch_time            = &volume_costs[volume_costs_idx+5];
  this->my_datamvt_time          = &volume_costs[volume_costs_idx+6];
  if (this->tag*critical_path_breakdown_size>0){
    size_t critical_path_costs_idx   = num_critical_path_measures+this->tag*critical_path_breakdown_size*num_tracker_critical_path_measures;
    this->critical_path_bytes        = &critical_path_costs[critical_path_costs_idx];
    this->critical_path_wrd          = &critical_path_costs[critical_path_costs_idx+critical_path_breakdown_size];
    this->critical_path_msg          = &critical_path_costs[critical_path_costs_idx+critical_path_breakdown_size*2];
    this->critical_path_comm_time    = &critical_path_costs[critical_path_costs_idx+critical_path_breakdown_size*3];
    this->critical_path_synch_time   = &critical_path_costs[critical_path_costs_idx+critical_path_breakdown_size*4];
    this->critical_path_datamvt_time = &critical_path_costs[critical_path_costs_idx+critical_path_breakdown_size*5];
  } else{
    this->critical_path_bytes        = &scratch_pad;
    this->critical_path_synch_time   = &scratch_pad;
    this->critical_path_datamvt_time = &scratch_pad;
    this->critical_path_comm_time    = &scratch_pad;
    this->critical_path_msg          = &scratch_pad;
    this->critical_path_wrd          = &scratch_pad;
  }
}

tracker::tracker(std::string name_, int tag, std::function<std::pair<double,double>(int64_t,int)> cost_func_){
  this->cost_func = cost_func_;
  this->name = std::move(name_);
  this->tag = tag;
  this->set_cost_pointers();
  this->init();
}

tracker::tracker(tracker const& t){
  this->cost_func = t.cost_func;
  this->name = t.name;
  this->tag = t.tag;
  this->set_cost_pointers();
  this->init();
}

tracker::~tracker(){}

void tracker::start_synch(volatile double curTime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, int nbr_pe, int nbr_pe2){
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  this->save_time         = curTime;
  this->save_comp_time    = curTime - computation_timer;
  critical_path_costs[6] += this->save_comp_time;		// update critical path computation time
  critical_path_costs[7] += this->save_comp_time;		// update critical path runtime
  volume_costs[7]        += this->save_comp_time;		// update local computation time
  volume_costs[8]        += this->save_comp_time;		// update local runtime

  int el_size,p;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(cm, &p);
  this->last_nbytes = nbytes;
  this->last_p = p;

  volatile double init_time = MPI_Wtime();
  if (nbr_pe == -1){
    PMPI_Barrier(cm);
  }
  else {
    double sbuf, rbuf;
    sbuf = 0.;
    PMPI_Sendrecv(&sbuf, 1, MPI_DOUBLE, nbr_pe, internal_tag, &rbuf, 1, MPI_DOUBLE, nbr_pe, internal_tag, cm, MPI_STATUS_IGNORE);
  }
  this->last_barrier_time = MPI_Wtime() - init_time;

  // Propogate critical paths for all processes in communicator based on what each process has seen up until now (not including this communication)
  propagate_critical_path(cm, nbr_pe, nbr_pe2);
  if (last_nbr_pe == -1){
    PMPI_Barrier(cm);
  } else {
    double sbuf, rbuf;
    sbuf = 0.;
    PMPI_Sendrecv(&sbuf, 1, MPI_DOUBLE, nbr_pe, internal_tag, &rbuf, 1, MPI_DOUBLE, nbr_pe, internal_tag, cm, MPI_STATUS_IGNORE);
  }
  // start synchronization timer for communication routine
  this->last_start_time = MPI_Wtime();
}

void tracker::start_synch(){
  // Deal with synchronization time
  volatile double synchTime = MPI_Wtime();
  this->last_synch_time = synchTime-this->last_start_time;
  // start communication timer for communication routine
  this->last_start_time = MPI_Wtime();
  if (mode == 2){
    auto overhead_time = this->last_start_time - this->save_time;
    symbol_timers[symbol_stack.top()].exclusive_overhead_time.top() += overhead_time;
  }
}

void tracker::stop_synch(){
  volatile double new_time = MPI_Wtime();
  volatile double dt = new_time - this->last_start_time;	// complete communication time
  double datamvt_time = std::max(0.,(dt-this->last_synch_time));
  std::pair<double,double> dcost = cost_func(this->last_nbytes, this->last_p);

  if (mode == 1){
    *this->my_bytes        += this->last_nbytes;
    *this->my_synch_time   += this->last_synch_time;
    *this->my_datamvt_time += datamvt_time;
    *this->my_comm_time    += dt;
    *this->my_msg          += dcost.first;
    *this->my_wrd          += dcost.second;
    *this->my_bar_time     += this->last_barrier_time;
    for (size_t i=0; i<critical_path_breakdown_size; i++){
      *(this->critical_path_bytes+i)        += this->last_nbytes;
      *(this->critical_path_synch_time+i)   += this->last_synch_time;
      *(this->critical_path_datamvt_time+i) += datamvt_time;
      *(this->critical_path_comm_time+i)    += dt;
      *(this->critical_path_msg+i)          += dcost.first;
      *(this->critical_path_wrd+i)          += dcost.second;
    }
  }
  else if (mode == 2){
    // update all communication-related measures for the top symbol in stack
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[0] += this->last_nbytes;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[1] += dcost.second;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[2] += dcost.first;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[3] += dt;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[4] += this->last_synch_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[5] += datamvt_time;
  }

  critical_path_costs[0] += this->last_nbytes;		// update critical path bytes communicated
  critical_path_costs[1] += dcost.second;		// update critical path estimated communication cost
  critical_path_costs[2] += dcost.first;		// update critical path estimated synchronization cost
  critical_path_costs[3] += dt;				// update critical path communication time (for what this process has seen thus far)
  critical_path_costs[4] += this->last_synch_time;	// update critical path synchronization time
  critical_path_costs[5] += datamvt_time;		// update critical path data mvt time
  critical_path_costs[7] += dt;				// update critical path runtime

  volume_costs[0] += this->last_nbytes;			// update local bytes communication
  volume_costs[1] += dcost.second;			// update local estimated communication cost
  volume_costs[2] += dcost.first;			// update local estimated synchronization cost
  volume_costs[3] += this->last_barrier_time;		// update local barrier/idle time
  volume_costs[4] += dt;				// update local communication time (not volume until after the completion of the program)
  volume_costs[5] += this->last_synch_time;		// update local synchronization time
  volume_costs[6] += datamvt_time;			// update local data mvt time
  volume_costs[8] += this->last_barrier_time;		// update local runtime with idle time
  volume_costs[8] += dt;				// update local runtime

  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
  volume_costs[4] = volume_costs[4] > critical_path_costs[3] ? critical_path_costs[3] : volume_costs[4];
  volume_costs[5] = volume_costs[5] > critical_path_costs[4] ? critical_path_costs[4] : volume_costs[5];
  volume_costs[6] = volume_costs[6] > critical_path_costs[5] ? critical_path_costs[5] : volume_costs[6];
  volume_costs[7] = volume_costs[7] > critical_path_costs[6] ? critical_path_costs[6] : volume_costs[7];
  volume_costs[8] = volume_costs[8] > critical_path_costs[7] ? critical_path_costs[7] : volume_costs[8];

  // Prepare to leave interception and re-enter user code
  this->last_start_time = MPI_Wtime();
  computation_timer = this->last_start_time;
  if (mode == 2){
    auto overhead_time = this->last_start_time - new_time;
    symbol_timers[symbol_stack.top()].exclusive_overhead_time.top() += overhead_time;
  }
}

void tracker::start_block(volatile double curTime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, int nbr_pe, int nbr_pe2){
  
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  this->save_time = curTime;
  this->save_comp_time = curTime - computation_timer;
  this->last_cm = cm;
  this->last_nbr_pe = nbr_pe;
  this->last_nbr_pe2 = nbr_pe2;

  int el_size,p;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(cm, &p);
  this->last_nbytes = nbytes;
  this->last_p = p;

  // start timer for communication routine
  this->last_start_time = MPI_Wtime();
}

void tracker::start_block(){
  // Deal with synchronization time
  volatile double synchTime = MPI_Wtime();
  this->last_synch_time = synchTime-this->last_start_time;
  // start communication timer for communication routine
  this->last_start_time = MPI_Wtime();
  if (mode == 2){
    auto overhead_time = this->last_start_time - this->save_time;
    symbol_timers[symbol_stack.top()].exclusive_overhead_time.top() += overhead_time;
  }
}

// Used only for p2p communication. All blocking collectives use sychronous protocol
void tracker::stop_block(bool is_sender){
  volatile double new_time = MPI_Wtime();
  volatile double dt = new_time - this->last_start_time;	// complete communication time
  double datamvt_time = std::max(0.,(dt-this->last_synch_time));
  std::pair<double,double> dcost = cost_func(this->last_nbytes, this->last_p);

  if (mode == 1){
    *this->my_bytes        += this->last_nbytes;
    *this->my_synch_time   += this->last_synch_time;
    *this->my_datamvt_time += datamvt_time;
    *this->my_comm_time    += dt;
    *this->my_msg          += dcost.first;
    *this->my_wrd          += dcost.second;
    for (size_t i=0; i<critical_path_breakdown_size; i++){
      *(this->critical_path_bytes+i)        += this->last_nbytes;
      *(this->critical_path_synch_time+i)   += this->last_synch_time;
      *(this->critical_path_datamvt_time+i) += datamvt_time;
      *(this->critical_path_comm_time+i)    += dt;
      *(this->critical_path_msg+i)          += dcost.first;
      *(this->critical_path_wrd+i)          += dcost.second;
    }
  }
  else if (mode == 2){
    // update all communication-related measures for the top symbol in stack
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[0] += this->last_nbytes;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[1] += dcost.second;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[2] += dcost.first;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[3] += dt;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[4] += this->last_synch_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[5] += datamvt_time;
  }

  critical_path_costs[6] += this->save_comp_time;	// update critical path computation time
  critical_path_costs[7] += this->save_comp_time+dt;	// update critical path runtime
  volume_costs[7] += this->save_comp_time;		// update local computation time
  volume_costs[8] += this->save_comp_time+dt;		// update local runtime

  critical_path_costs[0] += this->last_nbytes;		// update critical path bytes communicated
  critical_path_costs[1] += dcost.second;		// update critical path estimated communication cost
  critical_path_costs[2] += dcost.first;		// update critical path estimated synchronization cost
  critical_path_costs[3] += dt;				// update critical path communication time (for what this process has seen thus far)
  critical_path_costs[4] += this->last_synch_time;	// update critical path synchronization time (for what this process has seen thus far)
  critical_path_costs[5] += datamvt_time;		// update critical path data movement time (for what this process has seen thus far)

  volume_costs[0] += this->last_nbytes;			// update local bytes communication
  volume_costs[1] += dcost.second;			// update local estimated communication cost
  volume_costs[2] += dcost.first;			// update local estimated synchronization cost
  volume_costs[4] += dt;				// update local communication time (not volume until after the completion of the program)
  volume_costs[5] += this->last_synch_time;		// update local synchronzation time (not volume until after the completion of the program)
  volume_costs[6] += datamvt_time;			// update local data movement time (not volume until after the completion of the program)

  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
  volume_costs[4] = volume_costs[4] > critical_path_costs[3] ? critical_path_costs[3] : volume_costs[4];
  volume_costs[5] = volume_costs[5] > critical_path_costs[4] ? critical_path_costs[4] : volume_costs[5];
  volume_costs[6] = volume_costs[6] > critical_path_costs[5] ? critical_path_costs[5] : volume_costs[6];
  volume_costs[7] = volume_costs[7] > critical_path_costs[6] ? critical_path_costs[6] : volume_costs[7];
  volume_costs[8] = volume_costs[8] > critical_path_costs[7] ? critical_path_costs[7] : volume_costs[8];

  // Sender sends critical path data to receiver. Receiver updates its critical path information.
  if (is_sender){
    PMPI_Send(&critical_path_costs[0],critical_path_costs.size(),MPI_DOUBLE,this->last_nbr_pe,internal_tag,this->last_cm);
    // Sender needs not wait for handshake with receiver, can continue back into user code
  }
  else{
    MPI_Status st;
    PMPI_Recv(&new_cs[0],critical_path_costs.size(),MPI_DOUBLE,this->last_nbr_pe,internal_tag,this->last_cm,&st);
    update_critical_path(&new_cs[0]);
  }

  // Prepare to leave interception and re-enter user code
  this->last_start_time = MPI_Wtime();
  computation_timer = this->last_start_time;
  if (mode == 2){
    auto overhead_time = this->last_start_time - new_time;
    symbol_timers[symbol_stack.top()].exclusive_overhead_time.top() += overhead_time;
  }
}

// Called by both nonblocking p2p and nonblocking collectives
void tracker::start_nonblock(volatile double curTime, MPI_Request* request, int64_t nelem, MPI_Datatype t, MPI_Comm cm, bool is_sender, int nbr_pe, int nbr_pe2){
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  this->save_time = curTime;
  this->save_comp_time = curTime - computation_timer;
  critical_path_costs[6] += this->save_comp_time;		// update critical path computation time
  critical_path_costs[7] += this->save_comp_time;		// update critical path runtime
  volume_costs[7]        += this->save_comp_time;		// update local computation time
  volume_costs[8]        += this->save_comp_time;		// update local runtime

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
  if (nbr_pe == -1){
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
      PMPI_Isend(&data[0],critical_path_costs.size(),MPI_DOUBLE,nbr_pe,internal_tag,cm,&internal_request);
    }
    else{
      PMPI_Irecv(&data[0],critical_path_costs.size(),MPI_DOUBLE,nbr_pe,internal_tag,cm,&internal_request);
    }
  }
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  internal_comm_info[*request] = std::make_pair(internal_request,is_sender);
  internal_comm_message[*request] = data;
  internal_comm_data[*request] = std::make_pair((double)nbytes,(double)p);
  internal_comm_track[*request] = this;
  if (mode == 2){
    auto overhead_time = this->last_start_time - this->save_time;
    symbol_timers[symbol_stack.top()].exclusive_overhead_time.top() += overhead_time;
  }
}

void tracker::stop_nonblock(MPI_Request* request, double comp_time, double comm_time){
  volatile double new_time = MPI_Wtime();
  auto comm_info_it = internal_comm_info.find(*request);
  auto comm_message_it = internal_comm_message.find(*request);
  auto comm_data_it = internal_comm_data.find(*request);
  auto comm_track_it = internal_comm_track.find(*request);
  assert(comm_info_it != internal_comm_info.end());
  assert(comm_message_it != internal_comm_message.end());
  assert(comm_data_it != internal_comm_data.end());
  assert(comm_track_it != internal_comm_track.end());

  // Before accumulating the cost of this communication into our critical path/volume measures, we
  //   must first finish the internal communication, which doesn't take into account the cost of this communication
  // The computation and communication time of the sender between its MPI_Isend and MPI_Wait cannot be tracked. The receiver
  //   will just use its own. This is technically ok I think, because the receiver isn't waiting on the receiver in any capacity.

  MPI_Request internal_request = comm_info_it->second.first;
  bool is_sender = comm_info_it->second.second;
  double nbytes = comm_data_it->second.first;
  double p = comm_data_it->second.second;
  double* data = comm_message_it->second;
  MPI_Status st;
  PMPI_Wait(&internal_request,&st);
  // note: nonblocking collectives must update
  if (!is_sender){
    update_critical_path(data);
  }
  // Both sender and receiver will now update its critical path with the data from the communication
  std::pair<double,double> dcost = cost_func(nbytes, p);

  if (mode == 1){
    *this->my_bytes        += nbytes;
    *this->my_synch_time   += 0;			// Nonblocking routines will have no synchronization time component
    *this->my_datamvt_time += comm_time;
    *this->my_comm_time    += comm_time;
    *this->my_msg          += dcost.first;
    *this->my_wrd          += dcost.second;
    for (size_t i=0; i<critical_path_breakdown_size; i++){
      *(this->critical_path_bytes+i)     +=nbytes;
      *(this->critical_path_comm_time+i) += comm_time;
      *(this->critical_path_msg+i)       += dcost.first;
      *(this->critical_path_wrd+i)       += dcost.second;
    }
  }
  else if (mode == 2){
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[0] += this->last_nbytes;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[1] += dcost.second;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[2] += dcost.first;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[3] += comm_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[4] += 0;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[5] += comm_time;
  }

  critical_path_costs[0] += nbytes;
  critical_path_costs[1] += dcost.second;
  critical_path_costs[2] += dcost.first;
  critical_path_costs[3] += comm_time;
  critical_path_costs[4] += 0;
  critical_path_costs[5] += comm_time;
  critical_path_costs[6] += comp_time;
  critical_path_costs[7] += comp_time+comm_time;

  volume_costs[0] += nbytes;
  volume_costs[1] += dcost.second;
  volume_costs[2] += dcost.first;
  volume_costs[4] += comm_time;
  volume_costs[5] += 0;
  volume_costs[6] += comm_time;
  volume_costs[7] += comp_time;
  volume_costs[8] += comp_time+comm_time;

  // Due to granularity of timing, if a per-process measure ever gets more expensive than a critical path measure, we set the per-process measure to the cp measure
  volume_costs[4] = volume_costs[4] > critical_path_costs[3] ? critical_path_costs[3] : volume_costs[4];
  volume_costs[5] = volume_costs[5] > critical_path_costs[4] ? critical_path_costs[4] : volume_costs[5];
  volume_costs[6] = volume_costs[6] > critical_path_costs[5] ? critical_path_costs[5] : volume_costs[6];
  volume_costs[7] = volume_costs[7] > critical_path_costs[6] ? critical_path_costs[6] : volume_costs[7];
  volume_costs[8] = volume_costs[8] > critical_path_costs[7] ? critical_path_costs[7] : volume_costs[8];

  internal_comm_info.erase(*request);
  free(comm_message_it->second);
  internal_comm_message.erase(*request);
  internal_comm_data.erase(*request);
  internal_comm_track.erase(*request);

  this->last_start_time = MPI_Wtime();
  computation_timer = this->last_start_time;
  if (mode == 2){
    auto overhead_time = this->last_start_time - new_time;
    symbol_timers[symbol_stack.top()].exclusive_overhead_time.top() += overhead_time;
  }
}

void update_critical_path(double* data){
  if (critical_path_breakdown_size>0){
    bool decisions[critical_path_breakdown_size];
    size_t breakdown_idx=0;
    size_t breakdown_size = critical_path_breakdown_size;	// prevents compiler warning
    for (int i=0; i<num_critical_path_measures; i++){
      if (critical_path_breakdown[i]) decisions[breakdown_idx++] = data[i] > critical_path_costs[i];
      critical_path_costs[i] = std::max(data[i],critical_path_costs[i]);
    }
    for (int i=num_critical_path_measures; i<critical_path_costs.size(); i++){
      int idx = (i-num_critical_path_measures)%breakdown_size;
      critical_path_costs[i] = (decisions[idx] ? data[i] : critical_path_costs[i]);
    }
  } else{
    for (int i=0; i<num_critical_path_measures; i++){
      critical_path_costs[i] = std::max(data[i],critical_path_costs[i]);
    }
  }
}

// Only used with synchronous communication protocol
void propagate_critical_path(MPI_Comm cm, int nbr_pe, int nbr_pe2){
  if (mode == 1){
    // First exchange the tracked routine critical path data
    if (nbr_pe == -1){
      MPI_Op op; MPI_Op_create((MPI_User_function*) propagate_critical_path_op,1,&op);
      PMPI_Allreduce(MPI_IN_PLACE, &critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, op, cm);
      MPI_Op_free(&op);
    }
    else {
      PMPI_Sendrecv(&critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, nbr_pe, internal_tag, &new_cs[0], critical_path_costs.size(),
        MPI_DOUBLE, nbr_pe, internal_tag, cm, MPI_STATUS_IGNORE);
      update_critical_path(&new_cs[0]);
      if (nbr_pe2 != -1 && nbr_pe2 != nbr_pe){
        PMPI_Sendrecv(&critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, nbr_pe2, internal_tag, &new_cs[0], critical_path_costs.size(),
          MPI_DOUBLE, nbr_pe2, internal_tag, cm, MPI_STATUS_IGNORE);
        update_critical_path(&new_cs[0]);
     }
    }
  }
  else if (mode == 2){
    // Note: challenge of updating symbol_stack. For local symbols that are not along the updated critical path, we need to just make their contributions zero, but still keep them in the maps.
    // Next, exchange the critical path metric, together with tracking the rank of the process that determines each critical path
    int rank; MPI_Comm_rank(cm,&rank); int true_rank; MPI_Comm_rank(MPI_COMM_WORLD,&true_rank);
    for (int i=0; i<num_critical_path_measures; i++){
      timer_cp_info_sender[i].first = critical_path_costs[i];
      timer_cp_info_sender[i].second = rank;
    }
    if (nbr_pe == -1){
      PMPI_Allreduce(&timer_cp_info_sender[0].first, &timer_cp_info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, MPI_MAXLOC, cm);
    }
    else {
      PMPI_Sendrecv(&timer_cp_info_sender[0].first, num_critical_path_measures, MPI_DOUBLE_INT, nbr_pe, internal_tag, &timer_cp_info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, nbr_pe, internal_tag, cm, MPI_STATUS_IGNORE);
      for (int i=0; i<num_critical_path_measures; i++){
        if (timer_cp_info_sender[i].first>timer_cp_info_receiver[i].first){timer_cp_info_receiver[i].second = rank;}
        else if (timer_cp_info_sender[i].first==timer_cp_info_receiver[i].first){
          if (timer_cp_info_sender[i].second < timer_cp_info_receiver[i].second){
            timer_cp_info_receiver[i].second = rank;
          } else{
            timer_cp_info_receiver[i].second = nbr_pe;
          }
        }
        timer_cp_info_receiver[i].first = std::max(timer_cp_info_sender[i].first, timer_cp_info_receiver[i].first);
      }
/*
      if (nbr_pe2 != -1 && nbr_pe2 != nbr_pe){
        PMPI_Sendrecv(&timer_cp_info_sender[0].first, num_critical_path_measures, MPI_DOUBLE_INT, nbr_pe2, internal_tag, &timer_cp_info_receiver[0].first, num_critical_path_measures, MPI_DOUBLE_INT, nbr_pe2, internal_tag, cm, MPI_STATUS_IGNORE);
        for (int i=0; i<7; i++){
          if (timer_cp_info_sender[i].first<timer_cp_info_receiver[i].first){timer_cp_info_receiver[i].second = nbr_pe2;}
          timer_cp_info_receiver[i].first = std::max(timer_cp_info_sender[i].first, timer_cp_info_receiver[i].first);
        }
      }
*/
    }

    for (int i=0; i<num_critical_path_measures; i++){
      critical_path_costs[i] = timer_cp_info_receiver[i].first;
    }
    // We consider only critical path runtime
    int ftimer_size = 0;
    if (rank==timer_cp_info_receiver[num_critical_path_measures-1].second){
      ftimer_size = symbol_timers.size();
    }
    if (nbr_pe == -1){
      PMPI_Allreduce(MPI_IN_PLACE,&ftimer_size,1,MPI_INT,MPI_SUM,cm);
    }
    else{
      if (rank != nbr_pe){
        if (rank==timer_cp_info_receiver[num_critical_path_measures-1].second){
          PMPI_Send(&ftimer_size,1,MPI_INT,nbr_pe,internal_tag,cm);
        }
        else{
          PMPI_Recv(&ftimer_size,1,MPI_INT,nbr_pe,internal_tag,cm,MPI_STATUS_IGNORE);
        }
      }
/*
      if (nbr_pe2 != -1 && nbr_pe2 != nbr_pe){
        if (rank==timer_cp_info_receiver[num_critical_path_measures-1].second){
          PMPI_Send(&ftimer_size,1,MPI_INT,nbr_pe2,internal_tag,cm);
        }
        else{
          PMPI_Recv(&ftimer_size,1,MPI_INT,nbr_pe2,internal_tag,cm,MPI_STATUS_IGNORE);
        }
      }
*/
    }

    std::vector<int> symbol_sizes(ftimer_size);
    if (rank==timer_cp_info_receiver[num_critical_path_measures-1].second){
      int symbol_offset = 0; int i=0;
      for (auto i=0; i<symbol_timers.size(); i++){
        symbol_sizes[i] = symbol_order[i].size();
        for (auto j=0; j<symbol_sizes[i]; j++){
          symbol_pad[symbol_offset+j] = symbol_order[i][j];
        }
        symbol_offset += symbol_sizes[i];
      }
    }
    else{
      for (auto i=0; i<ftimer_size; i++){
        symbol_sizes[i] = 0;
      }
    }
    if (nbr_pe == -1){
      PMPI_Allreduce(MPI_IN_PLACE,&symbol_sizes[0],ftimer_size,MPI_INT,MPI_SUM,cm);
    }
    else{
      if (rank != nbr_pe){
        if (rank==timer_cp_info_receiver[num_critical_path_measures-1].second){
          PMPI_Send(&symbol_sizes[0],ftimer_size,MPI_INT,nbr_pe,internal_tag,cm);
        }
        else{
          PMPI_Recv(&symbol_sizes[0],ftimer_size,MPI_INT,nbr_pe,internal_tag,cm,MPI_STATUS_IGNORE);
        }
      }
/*
      if (nbr_pe2 != -1 && nbr_pe2 != nbr_pe){
        if (rank==timer_cp_info_receiver[num_critical_path_measures-1].second){
          PMPI_Send(&symbol_sizes[0],ftimer_size,MPI_INT,nbr_pe2,internal_tag,cm);
        }
        else{
          PMPI_Recv(&symbol_sizes[0],ftimer_size,MPI_INT,nbr_pe2,internal_tag,cm,MPI_STATUS_IGNORE);
        }
      }
*/
    }

    int num_chars = 0;
    for (auto i=0; i<ftimer_size; i++){
      num_chars += symbol_sizes[i];
    }
    if (rank == timer_cp_info_receiver[num_critical_path_measures-1].second){
      if (nbr_pe == -1){
        PMPI_Bcast(&symbol_timer_pad_local[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,rank,cm);
        PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,rank,cm);
      }
      else{
        if (rank != nbr_pe){
          PMPI_Send(&symbol_timer_pad_local[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,nbr_pe,internal_tag,cm);
          PMPI_Send(&symbol_pad[0],num_chars,MPI_CHAR,nbr_pe,internal_tag,cm);
/*
        if (nbr_pe2 != -1 && nbr_pe2 != nbr_pe){
          PMPI_Send(&symbol_timer_pad_local[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,nbr_pe2,internal_tag,cm);
          PMPI_Send(&symbol_pad[0],num_chars,MPI_CHAR,nbr_pe2,internal_tag,cm);
        }
*/
        }
      }
    }
    else{
      if (nbr_pe == -1){
        PMPI_Bcast(&symbol_timer_pad_global[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,timer_cp_info_receiver[num_critical_path_measures-1].second,cm);
        PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,timer_cp_info_receiver[num_critical_path_measures-1].second,cm);
      }
      else{
        if (rank != nbr_pe){
          PMPI_Recv(&symbol_timer_pad_global[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,nbr_pe,internal_tag,cm,MPI_STATUS_IGNORE);
          PMPI_Recv(&symbol_pad[0],num_chars,MPI_CHAR,nbr_pe,internal_tag,cm,MPI_STATUS_IGNORE);
/*
        if (nbr_pe2 != -1 && nbr_pe2 != nbr_pe){
          PMPI_Recv(&symbol_timer_pad_local[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,nbr_pe2,internal_tag,cm,MPI_STATUS_IGNORE);
          PMPI_Recv(&symbol_pad[0],num_chars,MPI_CHAR,nbr_pe2,internal_tag,cm,MPI_STATUS_IGNORE);
        }
*/
        }
      }
      if (rank != nbr_pe){
        int symbol_offset = 0;
        for (int i=0; i<ftimer_size; i++){
          auto reconstructed_symbol = std::string(symbol_pad.begin()+symbol_offset,symbol_pad.begin()+symbol_offset+symbol_sizes[i]);
          if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
            symbol_timers[reconstructed_symbol] = ftimer(reconstructed_symbol);
            symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
          }
          *symbol_timers[reconstructed_symbol].acc_numcalls = symbol_timer_pad_global[(num_ftimer_measures*num_critical_path_measures+1)*i];
          for (int j=0; j<num_critical_path_measures; j++){
            *symbol_timers[reconstructed_symbol].acc_measure[j] = symbol_timer_pad_global[(num_ftimer_measures*num_critical_path_measures+1)*i+2*j+1];
            *symbol_timers[reconstructed_symbol].acc_excl_measure[j] = symbol_timer_pad_global[(num_ftimer_measures*num_critical_path_measures+1)*i+2*(j+1)];
          }
          symbol_timers[reconstructed_symbol].has_been_processed = true;
          symbol_offset += symbol_sizes[i];
        }
        // Now cycle through and find the symbols that were not processed and set their accumulated measures to 0
        for (auto& it : symbol_timers){
          if (it.second.has_been_processed){ it.second.has_been_processed = false; }
          else{
            *it.second.acc_numcalls = 0;
            for (int j=0; j<num_critical_path_measures; j++){
              *it.second.acc_measure[j] = 0;
              *it.second.acc_excl_measure[j] = 0;
            }
          }
        }
      }
    }
  }
}

// Note: this function should be called once per start/stop, else it will double count
void compute_volume(MPI_Comm cm){
  size_t j=0;
  for (size_t i=0; i<num_volume_measures; i++){
    if (i!=3){// skip idle time
      max_per_process_costs[j++] = volume_costs[i];
    }
  }
  PMPI_Allreduce(MPI_IN_PLACE, &max_per_process_costs[0], max_per_process_costs.size(), MPI_DOUBLE, MPI_MAX, cm);
  PMPI_Allreduce(MPI_IN_PLACE, &volume_costs[0], volume_costs.size(), MPI_DOUBLE, MPI_SUM, cm);

  if (mode==2){
    //.. communicate timer stuff. Find which dude has max
  }
}

void tracker::set_header(){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if (*this->my_bytes != 0.){
    std::vector<double> vec(1);
    save_info[this->name] = std::move(vec);
  }
}

void tracker::set_critical_path_costs(size_t idx){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if ((*this->my_bytes != 0.) && (critical_path_breakdown_size>0)){
    std::vector<double> vec(num_tracker_critical_path_measures);
    vec[0] = *(this->critical_path_bytes+idx);
    vec[1] = *(this->critical_path_wrd+idx);
    vec[2] = *(this->critical_path_msg+idx);
    vec[3] = *(this->critical_path_comm_time+idx);
    vec[4] = *(this->critical_path_synch_time+idx);
    vec[5] = *(this->critical_path_datamvt_time+idx);
    save_info[this->name] = std::move(vec);
  }
}

void tracker::set_volume_costs(){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if (*this->my_bytes != 0.){
    std::vector<double> vec(num_tracker_volume_measures);
    vec[0] = *this->my_bytes;
    vec[1] = *this->my_wrd;
    vec[2] = *this->my_msg;
    vec[3] = *this->my_bar_time;
    vec[4] = *this->my_comm_time;
    vec[5] = *this->my_synch_time;
    vec[6] = *this->my_datamvt_time;
    save_info[this->name] = std::move(vec);
  }
}

ftimer::ftimer(std::string name_){
  assert(critter::internal::mode==2);
  assert(name_.size() <= max_timer_name_length);
  this->name = std::move(name_);
  this->acc_numcalls = &symbol_timer_pad_local[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)]; *this->acc_numcalls = 0;
  for (auto i=0; i<num_critical_path_measures; i++){
    this->acc_measure[i]      = &symbol_timer_pad_local[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)+2*i+1]; *acc_measure[i] = 0.;
    this->acc_excl_measure[i] = &symbol_timer_pad_local[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)+2*(i+1)]; *acc_excl_measure[i] = 0.;
  }
  this->has_been_processed = false;
}

void ftimer::start(){
  assert(critter::internal::mode==2);
  symbol_stack.push(this->name);
  this->start_timer.push(MPI_Wtime());
  this->exclusive_overhead_time.push(0.0);
  this->inclusive_overhead_time.push(0.0);
  this->inclusive_measure.push({0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0});
  this->exclusive_measure.push({0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0});
}

void ftimer::stop(){
  assert(critter::internal::mode==2);
  assert(this->start_timer.size()>0);
  // delta_time represents the symbol's exclusive time (even if the symbol is stacked recursively)
  this->inclusive_overhead_time.top() += this->exclusive_overhead_time.top();
  volatile double delta_time = MPI_Wtime() - this->start_timer.top() - this->inclusive_measure.top()[num_critical_path_measures-1] - this->inclusive_overhead_time.top();
  this->exclusive_measure.top()[num_critical_path_measures-1] += delta_time;
  this->exclusive_measure.top()[num_critical_path_measures-2] += (this->exclusive_measure.top()[num_critical_path_measures-1]-this->exclusive_measure.top()[num_critical_path_measures-5]);
  for (auto i=0; i<num_critical_path_measures; i++){
    this->inclusive_measure.top()[i] += this->exclusive_measure.top()[i];
    *this->acc_excl_measure[i] += this->exclusive_measure.top()[i];
    // branch below prevents overcounting of recursive symbols
    if (this->start_timer.size() == 1){
      *this->acc_measure[i]      += this->inclusive_measure.top()[i];
    }
  }
  *this->acc_numcalls = *this->acc_numcalls + 1.;
  auto save_inclusive_info = this->inclusive_measure.top();
  auto save_overhead_time  = this->inclusive_overhead_time.top();
  auto old_name = symbol_stack.top();
  this->start_timer.pop();
  this->exclusive_overhead_time.pop();
  this->inclusive_overhead_time.pop();
  this->inclusive_measure.pop();
  this->exclusive_measure.pop();
  symbol_stack.pop();
  if (symbol_stack.size() > 0){
    for (auto i=0; i<num_critical_path_measures; i++){
      symbol_timers[symbol_stack.top()].inclusive_measure.top()[i] += save_inclusive_info[i];
    }
    symbol_timers[symbol_stack.top()].inclusive_overhead_time.top() += save_overhead_time;
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
void print_header(StreamType& Stream, size_t num_inputs){
  for (size_t idx = 0; idx < (num_inputs+1); idx++){
    if (idx != 0){
      Stream << "\t";
    }
    Stream << "Input";
  }
  Stream << "\tNumBytes\tEstimatedCommCost\tEstimatedSynchCost\tCommunicationTime\tSynchronizationTime\tDataMvtTime\tComputationTime\tRunTime";// critical path
  Stream << "\tNumBytes\tEstimatedCommCost\tEstimatedSynchCost\tCommunicationTime\tSynchronizationTime\tDataMvtTime\tComputationTime\tRunTime";// per-process
  Stream << "\tNumBytes\tEstimatedCommCost\tEstimatedSynchCost\tIdleTime\tCommunicationTime\tSynchronizationTime\tDataMvtTime\tComputationTime\tRunTime";// volume
  for (auto i=0; i<num_tracker_critical_path_measures*critical_path_breakdown_size+num_tracker_volume_measures;i++){
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
      for (size_t i=0; i<max_per_process_costs.size(); i++){
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
      for (auto i=0; i<num_critical_path_measures; i++){
        if (!critical_path_breakdown[i]) continue;
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
  if (mode<2){
    if (is_world_root){
      Stream << "\n\n";
      Stream << std::left << std::setw(25) << "Critical path:";
      Stream << std::left << std::setw(25) << "NumBytes";
      Stream << std::left << std::setw(25) << "EstCommCost";
      Stream << std::left << std::setw(25) << "EstSynchCost";
      Stream << std::left << std::setw(25) << "CommTime";
      Stream << std::left << std::setw(25) << "SynchTime";
      Stream << std::left << std::setw(25) << "DataMvtTime";
      Stream << std::left << std::setw(25) << "CompTime";
      Stream << std::left << std::setw(25) << "RunTime";
      Stream << "\n";
      Stream << std::left << std::setw(25) << "                  ";
      for (size_t i=0; i<num_critical_path_measures; i++){
        Stream << std::left << std::setw(25) << factor*critical_path_costs[i];
      }
      Stream << "\n\n";

      Stream << std::left << std::setw(25) << "Per-process max:";
      Stream << std::left << std::setw(25) << "NumBytes";
      Stream << std::left << std::setw(25) << "EstCommCost";
      Stream << std::left << std::setw(25) << "EstSynchCost";
      Stream << std::left << std::setw(25) << "CommTime";
      Stream << std::left << std::setw(25) << "SynchTime";
      Stream << std::left << std::setw(25) << "DataMvtTime";
      Stream << std::left << std::setw(25) << "CompTime";
      Stream << std::left << std::setw(25) << "RunTime";
      Stream << "\n";
      Stream << std::left << std::setw(25) << "                  ";
      for (size_t i=0; i<max_per_process_costs.size(); i++){
        Stream << std::left << std::setw(25) << factor*max_per_process_costs[i];
      }
      Stream << "\n\n";

      Stream << std::left << std::setw(25) << "Volume:";
      Stream << std::left << std::setw(25) << "NumBytes";
      Stream << std::left << std::setw(25) << "EstCommCost";
      Stream << std::left << std::setw(25) << "EstSynchCost";
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
      for (auto i=0; i<num_critical_path_measures; i++){
        if (!critical_path_breakdown[i]) continue;
        for (int j=0; j<list_size; j++){
          list[j]->set_critical_path_costs(breakdown_idx);
        }
        breakdown_idx++;
        if (i==0){
          Stream << std::left << std::setw(25) << "NumBytes:";
        } else if (i==1){
          Stream << std::left << std::setw(25) << "EstCommCost:";
        } else if (i==2){
          Stream << std::left << std::setw(25) << "EstSynchCost:";
        } else if (i==3){
          Stream << std::left << std::setw(25) << "CommTime:";
        } else if (i==4){
          Stream << std::left << std::setw(25) << "SynchTime:";
        } else if (i==5){
          Stream << std::left << std::setw(25) << "DataMvtTime:";
        } else if (i==6){
          Stream << std::left << std::setw(25) << "CompTime:";
        } else if (i==7){
          Stream << std::left << std::setw(25) << "RunTime:";
        }
        Stream << std::left << std::setw(25) << "NumBytes";
        Stream << std::left << std::setw(25) << "EstCommCost";
        Stream << std::left << std::setw(25) << "EstSynchCost";
        Stream << std::left << std::setw(25) << "CommTime";
        Stream << std::left << std::setw(25) << "SynchTime";
        Stream << std::left << std::setw(25) << "DataMvtTime";
        for (auto& it : save_info){
          Stream << "\n";
          Stream << std::left << std::setw(25) << it.first;
          for (size_t j=0; j<num_tracker_critical_path_measures; j++){
            Stream << std::left << std::setw(25) << factor*it.second[j];
          }
        }
        Stream << "\n\n";
      }
      for (int i=0; i<list_size; i++){
        list[i]->set_volume_costs();
      }
      Stream << std::left << std::setw(25) << "Volume:";
      Stream << std::left << std::setw(25) << "NumBytes";
      Stream << std::left << std::setw(25) << "EstCommCost";
      Stream << std::left << std::setw(25) << "EstSynchCost";
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
      for (auto i=num_critical_path_measures-1; i>=0; i--){
        Stream << "\n" << std::left << std::setw(max_timer_name_length) << critical_path_measure_names[i];
        Stream << std::left << std::setw(25) << "number of calls";
        Stream << std::left << std::setw(25) << "exclusive";
        Stream << std::left << std::setw(25) << "exclusive %";
        Stream << std::left << std::setw(25) << "inclusive";
        Stream << std::left << std::setw(25) << "inclusive %";
        for (auto& it : symbol_timers){
          assert(it.second.start_timer.size() == 0);
          Stream << "\n" << std::left << std::setw(max_timer_name_length) << it.second.name;
          Stream << std::left << std::setw(25) << *it.second.acc_numcalls;
          Stream << std::left << std::setw(25) << *it.second.acc_excl_measure[i];
          Stream << std::left << std::setw(25) << std::setprecision(4) << 100.*(critical_path_costs[i] == 0. ? 0.00 : *it.second.acc_excl_measure[i]/critical_path_costs[i]);
          Stream << std::left << std::setw(25) << *it.second.acc_measure[i];
          Stream << std::left << std::setw(25) << std::setprecision(4) << 100.*(critical_path_costs[i] == 0. ? 0.00 : *it.second.acc_measure[i]/critical_path_costs[i]);
        }
        Stream << "\n" << std::left << std::setw(max_timer_name_length) << "total";
        Stream << std::left << std::setw(25) << "";
        Stream << std::left << std::setw(25) << critical_path_costs[i];
        Stream << std::left << std::setw(25) << "";
        Stream << std::left << std::setw(25) << critical_path_costs[i];
        Stream << std::left << std::setw(25) << "";
        Stream << "\n";
      }
    }
  }
}
};

void start(size_t mode){
  assert(mode>=0 && mode < 3);
  assert(internal::internal_comm_info.size() == 0);
  internal::mode=mode;
  internal::critical_path_measure_names = {"num bytes","est comm cost","est synch cost","comm time","synch time","datamvt time","comp time","runtime"};
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
  assert(mode==internal::mode);
  assert(internal::internal_comm_info.size() == 0);
  internal::critical_path_costs[6]+=(last_time-internal::computation_timer);	// update critical path computation time
  internal::critical_path_costs[7]+=(last_time-internal::computation_timer);	// update critical path runtime
  internal::volume_costs[7]+=(last_time-internal::computation_timer);	// update computation time volume
  internal::volume_costs[8]+=(last_time-internal::computation_timer);	// update runtime volume

  internal::propagate_critical_path(MPI_COMM_WORLD,-1,-1);
  internal::compute_volume(MPI_COMM_WORLD);

  internal::record(std::cout,factor);
  if (internal::flag) {internal::record(internal::stream,factor);}

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
