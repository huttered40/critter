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
std::map<MPI_Request,std::pair<MPI_Comm,int>> internal_comm_comm;
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
std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_local_cp;
std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_global_cp;
std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_pp;
std::array<double,(num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols> symbol_timer_pad_vol;
std::unordered_map<std::string,ftimer> symbol_timers;
std::stack<std::string> symbol_stack;
std::array<std::string,max_num_symbols> symbol_order;
std::array<std::string,num_critical_path_measures> critical_path_measure_names;
double_int timer_info_sender[num_critical_path_measures];
double_int timer_info_receiver[num_critical_path_measures];

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

void tracker::start_synch(volatile double curTime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, int partner1, int partner2){

  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  this->save_comp_time    = curTime - computation_timer;
  critical_path_costs[6] += this->save_comp_time;		// update critical path computation time
  critical_path_costs[7] += this->save_comp_time;		// update critical path runtime
  volume_costs[7]        += this->save_comp_time;		// update local computation time
  volume_costs[8]        += this->save_comp_time;		// update local runtime
  if (mode == 2){
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[num_critical_path_measures-1] += (curTime - symbol_timers[symbol_stack.top()].start_timer.top());
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[num_critical_path_measures-2] += (curTime - symbol_timers[symbol_stack.top()].start_timer.top());
  }

  int el_size,p;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(cm, &p);
  this->last_nbytes = nbytes;
  this->last_p = p;

  volatile double init_time = MPI_Wtime();
  if (partner1 == -1){
    PMPI_Barrier(cm);
  }
  else {
    double sbuf=0.; double rbuf=0.;
    PMPI_Sendrecv(&sbuf, 1, MPI_DOUBLE, partner1, internal_tag, &rbuf, 1, MPI_DOUBLE, partner1, internal_tag, cm, MPI_STATUS_IGNORE);
    if ((partner2 != -1) && (partner2 != partner2)){
      PMPI_Sendrecv(&sbuf, 1, MPI_DOUBLE, partner2, internal_tag, &rbuf, 1, MPI_DOUBLE, partner2, internal_tag, cm, MPI_STATUS_IGNORE);
    }
  }
  this->last_barrier_time = MPI_Wtime() - init_time;

  // Propogate critical paths for all processes in communicator based on what each process has seen up until now (not including this communication)
  propagate_critical_path_synch(cm, partner1);
  if ((partner2 != -1) && (partner2 != partner2)) propagate_critical_path_synch(cm, partner2);
  if (partner1 == -1){
    PMPI_Barrier(cm);
  } else {
    double sbuf=0.; double rbuf=0.;
    PMPI_Sendrecv(&sbuf, 1, MPI_DOUBLE, partner1, internal_tag, &rbuf, 1, MPI_DOUBLE, partner1, internal_tag, cm, MPI_STATUS_IGNORE);
    if (partner2 != -1 && partner2 != partner2) PMPI_Sendrecv(&sbuf, 1, MPI_DOUBLE, partner2, internal_tag, &rbuf, 1, MPI_DOUBLE, partner2, internal_tag, cm, MPI_STATUS_IGNORE);
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
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[7] += dt;
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
    symbol_timers[symbol_stack.top()].start_timer.top() = this->last_start_time;
  }
}

void tracker::start_block(volatile double curTime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, int partner){
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  if (mode == 2){
    this->save_time = curTime - symbol_timers[symbol_stack.top()].start_timer.top();
  }
  this->save_comp_time = curTime - computation_timer;
  this->last_cm = cm;
  this->last_partner = partner;

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
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[6] += this->save_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[7] += this->save_time+dt;
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

  // Exchange the tracked routine critical path data
  propagate_critical_path_blocking(this->last_cm,this->last_partner,is_sender);

  // Prepare to leave interception and re-enter user code
  this->last_start_time = MPI_Wtime();
  computation_timer = this->last_start_time;
  if (mode == 2){
    symbol_timers[symbol_stack.top()].start_timer.top() = this->last_start_time;
  }
}

// Called by both nonblocking p2p and nonblocking collectives
void tracker::start_nonblock(volatile double curTime, MPI_Request* request, int64_t nelem, MPI_Datatype t, MPI_Comm cm, bool is_sender, int partner){
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical path-path yet or that will screw up calculation of overlap!
  if (mode == 2){
    this->save_time = curTime - symbol_timers[symbol_stack.top()].start_timer.top();
  }
  this->save_comp_time = curTime - computation_timer;
  critical_path_costs[6] += this->save_comp_time;		// update critical path computation time
  critical_path_costs[7] += this->save_comp_time;		// update critical path runtime
  volume_costs[7]        += this->save_comp_time;		// update local computation time
  volume_costs[8]        += this->save_comp_time;		// update local runtime
  if (mode == 2){
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[num_critical_path_measures-1] += this->save_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[num_critical_path_measures-2] += this->save_time;
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

void tracker::stop_nonblock(MPI_Request* request, double comp_time, double comm_time){
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
  propagate_critical_path_nonblocking(data,internal_request,cm,partner,is_sender);

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
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[6] += comp_time;
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[7] += (comp_time+comm_time);
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

void propagate_critical_path_synch(MPI_Comm cm, int partner){
  if (mode == 1){
    // First exchange the tracked routine critical path data
    if (partner == -1){
      MPI_Op op; MPI_Op_create((MPI_User_function*) propagate_critical_path_op,1,&op);
      PMPI_Allreduce(MPI_IN_PLACE, &critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, op, cm);
      MPI_Op_free(&op);
    }
    else {
      PMPI_Sendrecv(&critical_path_costs[0], critical_path_costs.size(), MPI_DOUBLE, partner, internal_tag, &new_cs[0], critical_path_costs.size(),
        MPI_DOUBLE, partner, internal_tag, cm, MPI_STATUS_IGNORE);
      update_critical_path(&new_cs[0]);
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

    for (int i=0; i<num_critical_path_measures; i++){
      critical_path_costs[i] = timer_info_receiver[i].first;
    }
    // We consider only critical path runtime
    std::array<int,2> ftimer_size = {0,0};
    if (rank==timer_info_receiver[num_critical_path_measures-1].second){
      ftimer_size[0] = symbol_timers.size();
      ftimer_size[1] = symbol_stack.size() > 0 ? symbol_timers[symbol_stack.top()].exclusive_contributions.top().size() : 0;
    }
    if (partner == -1){
      PMPI_Allreduce(MPI_IN_PLACE,&ftimer_size[0],2,MPI_INT,MPI_SUM,cm);
    }
    else{
      if (rank != partner){
        if (rank==timer_info_receiver[num_critical_path_measures-1].second){
          PMPI_Send(&ftimer_size[0],2,MPI_INT,partner,internal_tag,cm);
        }
        else{
          PMPI_Recv(&ftimer_size[0],2,MPI_INT,partner,internal_tag,cm,MPI_STATUS_IGNORE);
        }
      }
    }

    std::vector<int> symbol_sizes(ftimer_size[0]+ftimer_size[1],0);
    std::vector<double> exclusive_contributions(ftimer_size[1]*num_critical_path_measures,0);
    if (rank==timer_info_receiver[num_critical_path_measures-1].second){
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
        if (rank==timer_info_receiver[num_critical_path_measures-1].second){
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
    if (rank == timer_info_receiver[num_critical_path_measures-1].second){
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
        PMPI_Bcast(&symbol_timer_pad_global_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size[0],MPI_DOUBLE,timer_info_receiver[num_critical_path_measures-1].second,cm);
        PMPI_Bcast(&exclusive_contributions[0], ftimer_size[1]*num_critical_path_measures,MPI_DOUBLE,timer_info_receiver[num_critical_path_measures-1].second,cm);
        PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,timer_info_receiver[num_critical_path_measures-1].second,cm);
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
            symbol_timers[reconstructed_symbol] = ftimer(reconstructed_symbol);
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
          symbol_timers[symbol_stack.top()].exclusive_contributions.top()[reconstructed_symbol] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
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

void propagate_critical_path_blocking(MPI_Comm cm, int partner, bool is_sender){
  if (mode == 1){
    // Sender sends critical path data to receiver. Receiver updates its critical path information.
    if (is_sender){
      PMPI_Send(&critical_path_costs[0],critical_path_costs.size(),MPI_DOUBLE,partner,internal_tag,cm);
      // Sender needs not wait for handshake with receiver, can continue back into user code
    }
    else{
      PMPI_Recv(&new_cs[0],critical_path_costs.size(),MPI_DOUBLE,partner,internal_tag,cm,MPI_STATUS_IGNORE);
      update_critical_path(&new_cs[0]);
    }
  }
  else if (mode == 2){
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
            symbol_timers[reconstructed_symbol] = ftimer(reconstructed_symbol);
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
          symbol_timers[symbol_stack.top()].exclusive_contributions.top()[reconstructed_symbol] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
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

void propagate_critical_path_nonblocking(double* data, MPI_Request internal_request, MPI_Comm cm, int partner, bool is_sender){
  // First exchange the tracked routine critical path data
  MPI_Status st;
  PMPI_Wait(&internal_request,&st);
  if (mode == 1){
    if (!is_sender){
      update_critical_path(data);
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

    for (int i=0; i<num_critical_path_measures; i++){
      critical_path_costs[i] = timer_info_receiver[i].first;
    }
    // We consider only critical path runtime

    std::array<int,2> ftimer_size = {0,0};
    if (rank==timer_info_receiver[num_critical_path_measures-1].second){
      ftimer_size[0] = symbol_timers.size();
      ftimer_size[1] = symbol_stack.size() > 0 ? symbol_timers[symbol_stack.top()].exclusive_contributions.top().size() : 0;
    }
    if (partner == -1){
      PMPI_Allreduce(MPI_IN_PLACE,&ftimer_size[0],2,MPI_INT,MPI_SUM,cm);
    }
    else{
      if (rank != partner){
        if (rank==timer_info_receiver[num_critical_path_measures-1].second){
          PMPI_Send(&ftimer_size[0],2,MPI_INT,partner,internal_tag,cm);
        }
        else{
          PMPI_Recv(&ftimer_size[0],2,MPI_INT,partner,internal_tag,cm,MPI_STATUS_IGNORE);
        }
      }
    }

    std::vector<int> symbol_sizes(ftimer_size[0]+ftimer_size[1],0);
    std::vector<double> exclusive_contributions(ftimer_size[1]*num_critical_path_measures,0);
    if (rank==timer_info_receiver[num_critical_path_measures-1].second){
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
        if (rank==timer_info_receiver[num_critical_path_measures-1].second){
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
    if (rank == timer_info_receiver[num_critical_path_measures-1].second){
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
        PMPI_Bcast(&symbol_timer_pad_global_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size[0],MPI_DOUBLE,timer_info_receiver[num_critical_path_measures-1].second,cm);
        PMPI_Bcast(&exclusive_contributions[0], ftimer_size[1]*num_critical_path_measures,MPI_DOUBLE,timer_info_receiver[num_critical_path_measures-1].second,cm);
        PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,timer_info_receiver[num_critical_path_measures-1].second,cm);
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
            symbol_timers[reconstructed_symbol] = ftimer(reconstructed_symbol);
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
          symbol_timers[symbol_stack.top()].exclusive_contributions.top()[reconstructed_symbol] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
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
void compute_volume(MPI_Comm cm){
  if (mode<2){
    size_t j=0;
    for (size_t i=0; i<num_volume_measures; i++){
      if (i!=3){// skip idle time
        max_per_process_costs[j++] = volume_costs[i];
      }
    }
    PMPI_Allreduce(MPI_IN_PLACE, &max_per_process_costs[0], max_per_process_costs.size(), MPI_DOUBLE, MPI_MAX, cm);
    PMPI_Allreduce(MPI_IN_PLACE, &volume_costs[0], volume_costs.size(), MPI_DOUBLE, MPI_SUM, cm);
  }
  else if (mode==2){
    int rank; MPI_Comm_rank(cm,&rank);
    size_t j=0;
    for (size_t i=0; i<num_volume_measures; i++){
      if (i!=3){// skip idle time
        max_per_process_costs[j] = volume_costs[i];
        timer_info_sender[j].first = volume_costs[i];
        timer_info_sender[j].second = rank;
        j++;
      }
    }
    PMPI_Allreduce(&timer_info_sender[0].first, &timer_info_receiver[0].first, max_per_process_costs.size(), MPI_DOUBLE_INT, MPI_MAXLOC, cm);
    PMPI_Allreduce(MPI_IN_PLACE, &max_per_process_costs[0], max_per_process_costs.size(), MPI_DOUBLE, MPI_MAX, cm);
    PMPI_Allreduce(MPI_IN_PLACE, &volume_costs[0], volume_costs.size(), MPI_DOUBLE, MPI_SUM, cm);
    // Using the rank determining the largest per-process max, each process gets its values
    int ftimer_size = 0;
    if (rank==timer_info_receiver[num_critical_path_measures-1].second){
      ftimer_size = symbol_timers.size();
    }
    PMPI_Allreduce(MPI_IN_PLACE,&ftimer_size,1,MPI_INT,MPI_SUM,cm);

    std::vector<int> symbol_sizes(ftimer_size,0);
    if (rank==timer_info_receiver[num_critical_path_measures-1].second){
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
    if (rank == timer_info_receiver[num_critical_path_measures-1].second){
      PMPI_Bcast(&symbol_timer_pad_local_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,rank,cm);
      PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,rank,cm);
    }
    else{
      PMPI_Bcast(&symbol_timer_pad_global_cp[0],(num_ftimer_measures*num_critical_path_measures+1)*ftimer_size,MPI_DOUBLE,timer_info_receiver[num_critical_path_measures-1].second,cm);
      PMPI_Bcast(&symbol_pad[0],num_chars,MPI_CHAR,timer_info_receiver[num_critical_path_measures-1].second,cm);
      int symbol_offset = 0;
      for (int i=0; i<ftimer_size; i++){
        auto reconstructed_symbol = std::string(symbol_pad.begin()+symbol_offset,symbol_pad.begin()+symbol_offset+symbol_sizes[i]);

        if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
          symbol_timers[reconstructed_symbol] = ftimer(reconstructed_symbol);
          symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
        }
        *symbol_timers[reconstructed_symbol].pp_numcalls = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i];
        for (int j=0; j<num_critical_path_measures; j++){
          *symbol_timers[reconstructed_symbol].pp_incl_measure[j] = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i+2*j+1];
          *symbol_timers[reconstructed_symbol].pp_excl_measure[j] = symbol_timer_pad_global_cp[(num_ftimer_measures*num_critical_path_measures+1)*i+2*(j+1)];
        }
        symbol_timers[reconstructed_symbol].has_been_processed = true;
        symbol_offset += symbol_sizes[i];
      }
      // Now cycle through and find the symbols that were not processed and set their accumulated measures to 0
      for (auto& it : symbol_timers){
        if (it.second.has_been_processed){ it.second.has_been_processed = false; }
        else{
          *it.second.pp_numcalls = 0;
          for (int j=0; j<num_critical_path_measures; j++){
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
  auto save_time = MPI_Wtime();
  assert(name_.size() <= max_timer_name_length);
  assert(symbol_timers.size() < max_num_symbols);
  this->name = std::move(name_);
  this->cp_numcalls = &symbol_timer_pad_local_cp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)]; *this->cp_numcalls = 0;
  this->pp_numcalls = &symbol_timer_pad_pp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)]; *this->pp_numcalls = 0;
  this->vol_numcalls = &symbol_timer_pad_vol[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)]; *this->vol_numcalls = 0;
  for (auto i=0; i<num_critical_path_measures; i++){
    this->cp_incl_measure[i] = &symbol_timer_pad_local_cp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)+2*i+1]; *cp_incl_measure[i] = 0.;
    this->cp_excl_measure[i] = &symbol_timer_pad_local_cp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)+2*(i+1)]; *cp_excl_measure[i] = 0.;
    this->pp_incl_measure[i] = &symbol_timer_pad_pp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)+2*i+1]; *pp_incl_measure[i] = 0.;
    this->pp_excl_measure[i] = &symbol_timer_pad_pp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)+2*(i+1)]; *pp_excl_measure[i] = 0.;
    this->vol_incl_measure[i] = &symbol_timer_pad_vol[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)+2*i+1]; *vol_incl_measure[i] = 0.;
    this->vol_excl_measure[i] = &symbol_timer_pad_vol[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)+2*(i+1)]; *vol_excl_measure[i] = 0.;
  }
  this->has_been_processed = false;
  critical_path_costs[6] += (save_time - computation_timer);		// update critical path computation time
  critical_path_costs[7] += (save_time - computation_timer);		// update critical path runtime
  volume_costs[7]        += (save_time - computation_timer);		// update local computation time
  volume_costs[8]        += (save_time - computation_timer);		// update local runtime
  computation_timer = MPI_Wtime();
}

void ftimer::start(){
  auto save_time = MPI_Wtime();
  if (symbol_stack.size()>0){
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[num_critical_path_measures-1] += (save_time-symbol_timers[symbol_stack.top()].start_timer.top());
    symbol_timers[symbol_stack.top()].exclusive_measure.top()[num_critical_path_measures-2] += (save_time-symbol_timers[symbol_stack.top()].start_timer.top());
  }
  symbol_stack.push(this->name);
  this->exclusive_contributions.push(typename decltype(this->exclusive_contributions)::value_type());
  this->exclusive_measure.push({0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0});
  // Explicit measure will never be finalize until the specific symbol's stack is 0 (nontrivial only for recursive nested symbols)
  critical_path_costs[6] += (save_time - computation_timer);		// update critical path computation time
  critical_path_costs[7] += (save_time - computation_timer);		// update critical path runtime
  volume_costs[7]        += (save_time - computation_timer);		// update local computation time
  volume_costs[8]        += (save_time - computation_timer);		// update local runtime
  computation_timer = MPI_Wtime();
  this->start_timer.push(computation_timer);
}

void ftimer::stop(){
  auto save_time = MPI_Wtime();
  assert(this->start_timer.size()>0);
  this->exclusive_measure.top()[num_critical_path_measures-1] += (save_time-this->start_timer.top());
  this->exclusive_measure.top()[num_critical_path_measures-2] += (save_time-this->start_timer.top());

  for (auto i=0; i<num_critical_path_measures; i++){
    *this->cp_excl_measure[i] += this->exclusive_measure.top()[i];
    *this->pp_excl_measure[i] += this->exclusive_measure.top()[i];
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

  critical_path_costs[6] += (save_time - computation_timer);		// update critical path computation time
  critical_path_costs[7] += (save_time - computation_timer);		// update critical path runtime
  volume_costs[7]        += (save_time - computation_timer);		// update local computation time
  volume_costs[8]        += (save_time - computation_timer);		// update local runtime
  computation_timer = MPI_Wtime();
  if (symbol_stack.size()>0) symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer;
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
      Stream << "***********************************************************************************************************************";
      for (auto i=num_critical_path_measures-1; i>=0; i--){
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
        for (auto& it : symbol_timers){
          assert(it.second.start_timer.size() == 0);
          Stream << "\n" << std::left << std::setw(max_timer_name_length) << it.second.name;
          Stream << std::left << std::setw(15) << *it.second.cp_numcalls;
          Stream << std::left << std::setw(15) << *it.second.cp_excl_measure[i];
          Stream << std::left << std::setw(15) << std::setprecision(4) << 100.*(critical_path_costs[i] == 0. ? 100.0 : *it.second.cp_excl_measure[i]/critical_path_costs[i]);
          Stream << std::left << std::setw(15) << *it.second.pp_numcalls;
          Stream << std::left << std::setw(15) << *it.second.pp_excl_measure[i];
          Stream << std::left << std::setw(15) << std::setprecision(4) << 100.*(max_per_process_costs[i] == 0. ? 0.0 : *it.second.pp_excl_measure[i]/max_per_process_costs[i]);
          Stream << std::left << std::setw(15) << *it.second.vol_numcalls;
          Stream << std::left << std::setw(15) << *it.second.vol_excl_measure[i];
          Stream << std::left << std::setw(15) << std::setprecision(4) << 100.*(volume_costs[i] == 0. ? 0.0 : *it.second.vol_excl_measure[i]/volume_costs[i]);
          cp_total_exclusive += *it.second.cp_excl_measure[i];
          pp_total_exclusive += *it.second.pp_excl_measure[i];
          vol_total_exclusive += *it.second.vol_excl_measure[i];
        }
        Stream << "\n" << std::left << std::setw(max_timer_name_length) << "total";
        Stream << std::left << std::setw(15) << "";
        Stream << std::left << std::setw(15) << cp_total_exclusive;
        Stream << std::left << std::setw(15) << 100.*cp_total_exclusive/critical_path_costs[i];
        Stream << std::left << std::setw(15) << "";
        Stream << std::left << std::setw(15) << pp_total_exclusive;
        Stream << std::left << std::setw(15) << 100.*pp_total_exclusive/max_per_process_costs[i];
        Stream << std::left << std::setw(15) << "";
        Stream << std::left << std::setw(15) << vol_total_exclusive;
        Stream << std::left << std::setw(15) << 100.*vol_total_exclusive/volume_costs[i];
        Stream << "\n";

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
        for (auto& it : symbol_timers){
          assert(it.second.start_timer.size() == 0);
          Stream << "\n" << std::left << std::setw(max_timer_name_length) << it.second.name;
          Stream << std::left << std::setw(15) << *it.second.cp_numcalls;
          Stream << std::left << std::setw(15) << *it.second.cp_incl_measure[i];
          Stream << std::left << std::setw(15) << std::setprecision(4) << 100.*(critical_path_costs[i] == 0. ? 100.0 : *it.second.cp_incl_measure[i]/critical_path_costs[i]);
          Stream << std::left << std::setw(15) << *it.second.pp_numcalls;
          Stream << std::left << std::setw(15) << *it.second.pp_incl_measure[i];
          Stream << std::left << std::setw(15) << std::setprecision(4) << 100.*(max_per_process_costs[i] == 0. ? 100.0 : *it.second.pp_incl_measure[i]/max_per_process_costs[i]);
          Stream << std::left << std::setw(15) << *it.second.vol_numcalls;
          Stream << std::left << std::setw(15) << *it.second.vol_incl_measure[i];
          Stream << std::left << std::setw(15) << std::setprecision(4) << 100.*(volume_costs[i] == 0. ? 100.0 : *it.second.vol_incl_measure[i]/volume_costs[i]);
          cp_total_inclusive = std::max(*it.second.cp_incl_measure[i],cp_total_inclusive);
          pp_total_inclusive = std::max(*it.second.pp_incl_measure[i],pp_total_inclusive);
          vol_total_inclusive = std::max(*it.second.vol_incl_measure[i],vol_total_inclusive);
        }
        Stream << "\n" << std::left << std::setw(max_timer_name_length) << "total";
        Stream << std::left << std::setw(15) << "";
        Stream << std::left << std::setw(15) << cp_total_inclusive;
        Stream << std::left << std::setw(15) << 100.*cp_total_inclusive/critical_path_costs[i];
        Stream << std::left << std::setw(15) << "";
        Stream << std::left << std::setw(15) << pp_total_inclusive;
        Stream << std::left << std::setw(15) << 100.*pp_total_inclusive/max_per_process_costs[i];
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

  internal::propagate_critical_path_synch(MPI_COMM_WORLD,-1);
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
