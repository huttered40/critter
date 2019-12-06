#include "critter.h"

namespace critter{
namespace internal{

void add_critical_path_data(int_int_double* in, int_int_double* inout, int* len, MPI_Datatype* dtype){
  int_int_double* invec = in;
  int_int_double* inoutvec = inout;
  for (int i=0; i<*len; i++){
    inoutvec[i].first = std::max(inoutvec[i].first,invec[i].first);
    inoutvec[i].second = std::max(inoutvec[i].second,invec[i].second);
    inoutvec[i].third = std::max(inoutvec[i].third,invec[i].third);
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
        &_MPI_Ialltoallv};

std::string stream_name,/*stream_track_name,*/file_name;
std::ofstream stream/*, stream_track*/;
bool track,flag,is_world_root,is_first_iter,need_new_line;

double computation_timer;
std::map<MPI_Request,std::pair<MPI_Request,bool>> internal_comm_info;
std::map<MPI_Request,double*> internal_comm_message;
std::map<MPI_Request,std::pair<double,double>> internal_comm_data;
std::map<MPI_Request,tracker*> internal_comm_track;
//std::vector<std::vector<int_int_double>> critical_paths(7);
std::array<double,14> costs;	// NumBytes,CommTime,IdleTime,EstCommCost,EstSynchCost,CompTime,RunTime
// Instead of printing out each Critter for each iteration individually, I will save them for each iteration, print out the iteration, and then clear before next iteration
std::map<std::string,std::tuple<double,double,double,double,double,double,double,double,double,double>> save_info;

double old_cs[5*list_size];
double new_cs[5*list_size];
double_int old_cp[7];
double_int new_cp[7];
double old_vp[7];
double new_vp[7];
int root_array[7];
int crit_path_size_array[7];


void tracker::init(){
  this->last_start_time  = -1.;
  this->my_bytes         = 0.;
  this->my_comm_time     = 0.;
  this->my_bar_time      = 0.;
  this->my_msg           = 0.;
  this->my_wrd           = 0.;
  this->critical_path_bytes       = 0.;
  this->critical_path_comm_time   = 0.;
  this->critical_path_bar_time    = 0.;
  this->critical_path_msg         = 0.;
  this->critical_path_wrd         = 0.;
  this->save_comp_time   = 0.;
  this->volume_bytes     = 0.;
  this->volume_comm_time = 0.;
  this->volume_bar_time  = 0.;
  this->volume_msg       = 0.;
  this->volume_wrd       = 0.;
}

tracker::tracker(std::string name_, int tag, std::function<std::pair<double,double>(int64_t,int)> cost_func_){
  this->cost_func = cost_func_;
  this->name = std::move(name_);
  this->tag = tag;
  this->init();
}

tracker::tracker(tracker const& t){
  this->cost_func = t.cost_func;
  this->name = t.name;
  this->tag = t.tag;
  this->init();
}

tracker::~tracker(){}

void tracker::start_synch(int64_t nelem, MPI_Datatype t, MPI_Comm cm, int nbr_pe, int nbr_pe2){
  //assert(this->last_start_time == -1.); //assert timer was not started twice without first being stopped
  
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical_pathical-path yet or that will screw up calculation of overlap!
  volatile double curTime = MPI_Wtime();
  this->save_comp_time = curTime - computation_timer;
  costs[5] += this->save_comp_time;		// update critical_pathical path computation time
  costs[6] += this->save_comp_time;		// update critical_pathical path runtime
  costs[12] += this->save_comp_time;		// update local computation time
  costs[13] += this->save_comp_time;		// update local runtime

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

  // Propogate critical_pathical paths for all processes in communicator based on what each process has seen up until now (not including this communication)
  propagate_critical_path(cm, nbr_pe, nbr_pe2);
  if (last_nbr_pe == -1){
    PMPI_Barrier(cm);
  } else {
    double sbuf, rbuf;
    sbuf = 0.;
    PMPI_Sendrecv(&sbuf, 1, MPI_DOUBLE, nbr_pe, internal_tag, &rbuf, 1, MPI_DOUBLE, nbr_pe, internal_tag, cm, MPI_STATUS_IGNORE);
  }
  // start timer for communication routine
  this->last_start_time = MPI_Wtime();
}

void tracker::stop_synch(){
  double dt = MPI_Wtime() - this->last_start_time;	// complete communication time
  std::pair<double,double> dcost = cost_func(this->last_nbytes, this->last_p);

  this->my_comm_time += dt;
  this->critical_path_comm_time += dt;
  this->my_bytes += this->last_nbytes;
  this->critical_path_bytes += this->last_nbytes;
  this->critical_path_msg += dcost.first;
  this->my_msg += dcost.first;
  this->critical_path_wrd += dcost.second;
  this->my_wrd += dcost.second;
  this->critical_path_bar_time += this->last_barrier_time;
  this->my_bar_time += this->last_barrier_time;

  costs[0] += this->last_nbytes;		// update critical_pathical path bytes communicated
  costs[1] += dt;				// update critical_pathical path communication time (for what this process has seen thus far)
  costs[2] += this->last_barrier_time;		// update critical_pathical path barrier/idle time
  costs[3] += dcost.second;			// update critical_pathical path estimated communication cost
  costs[4] += dcost.first;			// update critical_pathical path estimated synchronization cost
  costs[6] += dt;		// update critical_pathical path runtime
  costs[7] += this->last_nbytes;		// update local bytes communication
  costs[8] += dt;				// update local communication time (not volume until after the completion of the program)
  costs[9] += this->last_barrier_time;		// update local barrier/idle time
  costs[10] += dcost.second;			// update local estimated communication cost
  costs[11] += dcost.first;			// update local estimated synchronization cost
  costs[13] += dt;		// update local runtime

/*
  // Mark the local synchronization point before exchanging with its neighbors in the communicator
  critical_pathical_paths[0].emplace_back(int_int_double(this->tag,p,nbytes));
  critical_pathical_paths[1].emplace_back(int_int_double(this->tag,p,nbytes));
  critical_pathical_paths[2].emplace_back(int_int_double(this->tag,p,nbytes));
  critical_pathical_paths[3].emplace_back(int_int_double(this->tag,p,nbytes));
  critical_pathical_paths[4].emplace_back(int_int_double(this->tag,p,nbytes));
  critical_pathical_paths[5].emplace_back(int_int_double(this->tag,p,nbytes));
  critical_pathical_paths[6].emplace_back(int_int_double(this->tag,p,nbytes));
*/

  // Prepare to leave interception and re-enter user code
  this->last_start_time = MPI_Wtime();
  computation_timer = this->last_start_time;
}

void tracker::start_block(int64_t nelem, MPI_Datatype t, MPI_Comm cm, int nbr_pe, int nbr_pe2){
  //assert(this->last_start_time == -1.); //assert timer was not started twice without first being stopped
  
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical_pathical-path yet or that will screw up calculation of overlap!
  volatile double curTime = MPI_Wtime();
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

void tracker::stop_block(bool is_sender){
  double dt = MPI_Wtime() - this->last_start_time;	// complete communication time
  std::pair<double,double> dcost = cost_func(this->last_nbytes, this->last_p);

  this->my_comm_time += dt;
  this->critical_path_comm_time += dt;
  this->my_bytes += this->last_nbytes;
  this->critical_path_bytes += this->last_nbytes;
  this->critical_path_msg += dcost.first;
  this->my_msg += dcost.first;
  this->critical_path_wrd += dcost.second;
  this->my_wrd += dcost.second;

  costs[0] += this->last_nbytes;		// update critical_pathical path bytes communicated
  costs[1] += dt;				// update critical_pathical path communication time (for what this process has seen thus far)
  costs[3] += dcost.second;			// update critical_pathical path estimated communication cost
  costs[4] += dcost.first;			// update critical_pathical path estimated synchronization cost
  costs[5] += this->save_comp_time;		// update critical_pathical path computation time
  costs[6] += this->save_comp_time+dt;		// update critical_pathical path runtime
  costs[7] += this->last_nbytes;		// update local bytes communication
  costs[8] += dt;				// update local communication time (not volume until after the completion of the program)
  costs[10] += dcost.second;			// update local estimated communication cost
  costs[11] += dcost.first;			// update local estimated synchronization cost
  costs[12] += this->save_comp_time;		// update local computation time
  costs[13] += this->save_comp_time+dt;		// update local runtime

  // Sender sends critical_pathical path data to receiver. Receiver updates its critical_pathical path information.
  std::vector<double> data(10);
  if (is_sender){
    data[0]=this->critical_path_bytes; data[1]=this->critical_path_comm_time; data[2]=this->critical_path_msg; data[3]=this->critical_path_wrd;
    data[4]=costs[0]; data[5]=costs[1]; data[6]=costs[3]; data[7]=costs[4]; data[8]=costs[5]; data[9]=costs[6];
    PMPI_Send(&data[0],10,MPI_DOUBLE,this->last_nbr_pe,internal_tag,this->last_cm);
    // Sender needs not wait for handshake with receiver, can continue back into user code
  }
  else{
    MPI_Status st;
    PMPI_Recv(&data[0],10,MPI_DOUBLE,this->last_nbr_pe,internal_tag,this->last_cm,&st);
    this->critical_path_bytes=std::max(data[0],this->critical_path_bytes); this->critical_path_comm_time=std::max(data[1],this->critical_path_comm_time);
    this->critical_path_msg=std::max(data[2],this->critical_path_msg); this->critical_path_wrd=std::max(data[3],this->critical_path_wrd);
    costs[0]=std::max(data[4],costs[0]); costs[1]=std::max(data[5],costs[1]); costs[3]=std::max(data[6],costs[3]);
    costs[4]=std::max(data[7],costs[4]); costs[5]=std::max(data[8],costs[5]); costs[6]=std::max(data[9],costs[6]);
  }

  // Prepare to leave interception and re-enter user code
  this->last_start_time = MPI_Wtime();
  computation_timer = this->last_start_time;
}

void tracker::start_nonblock(MPI_Request* request, int64_t nelem, MPI_Datatype t, MPI_Comm cm, bool is_sender, int nbr_pe, int nbr_pe2){
  
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical_pathical-path yet or that will screw up calculation of overlap!
  volatile double curTime = MPI_Wtime();
  this->save_comp_time = curTime - computation_timer;
  costs[5] += this->save_comp_time;		// update critical_pathical path computation time
  costs[6] += this->save_comp_time;		// update critical_pathical path runtime
  costs[12] += this->save_comp_time;		// update local computation time
  costs[13] += this->save_comp_time;		// update local runtime

  int el_size,p;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  MPI_Comm_size(cm, &p);

  // Nonblocking communication to propogate the critical_pathical path from sender to receiver. Avoids tricky deadlock in intercepting MPI_Waitall
  // Unlike blocking protocol, Receiver does not need sender's critical_pathical path information to include the contribution from this current routine
  MPI_Request internal_request;
  double* data = (double*)malloc(sizeof(double)*10);
  // Save local data instead of immediately adding it to critical_pathical path, because this communication is not technically completed yet,
  //   and I do not want to corrupt critical_pathical path propogation in future communication that may occur before this nonblocking communication completes.
  if (nbr_pe == -1){
    data[0]=this->critical_path_bytes; data[1]=this->critical_path_comm_time; data[2]=this->critical_path_msg; data[3]=this->critical_path_wrd;
    data[4]=costs[0]; data[5]=costs[1]; data[6]=costs[3]; data[7]=costs[4]; data[8]=costs[5]; data[9]=costs[6];
    PMPI_Iallreduce(MPI_IN_PLACE,&data[0],10,MPI_DOUBLE,MPI_MAX,cm,&internal_request);
  } else{
    if (is_sender){
      data[0]=this->critical_path_bytes; data[1]=this->critical_path_comm_time; data[2]=this->critical_path_msg; data[3]=this->critical_path_wrd;
      data[4]=costs[0]; data[5]=costs[1]; data[6]=costs[3]; data[7]=costs[4]; data[8]=costs[5]; data[9]=costs[6];
      PMPI_Isend(&data[0],10,MPI_DOUBLE,nbr_pe,internal_tag,cm,&internal_request);
    }
    else{
      PMPI_Irecv(&data[0],10,MPI_DOUBLE,nbr_pe,internal_tag,cm,&internal_request);
    }
  }
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  internal_comm_info[*request] = std::make_pair(internal_request,is_sender);
  internal_comm_message[*request] = data;
  internal_comm_data[*request] = std::make_pair((double)nbytes,(double)p);
  internal_comm_track[*request] = this;
}

void tracker::stop_nonblock(MPI_Request* request, double comp_time, double comm_time){
  auto comm_info_it = internal_comm_info.find(*request);\
  auto comm_message_it = internal_comm_message.find(*request);\
  auto comm_data_it = internal_comm_data.find(*request);\
  auto comm_track_it = internal_comm_track.find(*request);\
  assert(comm_info_it != internal_comm_info.end());\
  assert(comm_message_it != internal_comm_message.end());\
  assert(comm_data_it != internal_comm_data.end());\
  assert(comm_track_it != internal_comm_track.end());\

  // Before accumulating the cost of this communication into our critical_pathical path/volume measures, we
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
  if (!is_sender){
    this->critical_path_bytes=std::max(data[0],this->critical_path_bytes); this->critical_path_comm_time=std::max(data[1],this->critical_path_comm_time);
    this->critical_path_msg=std::max(data[2],this->critical_path_msg); this->critical_path_wrd=std::max(data[3],this->critical_path_wrd);
    costs[0]=std::max(data[4],costs[0]); costs[1]=std::max(data[5],costs[1]); costs[3]=std::max(data[6],costs[3]);
    costs[4]=std::max(data[7],costs[4]); costs[5]=std::max(data[8],costs[5]); costs[6]=std::max(data[9],costs[6]);
  }
  // Both sender and receiver will now update its critical_pathical path with the data from the communication
  std::pair<double,double> dcost = cost_func(nbytes, p);
  this->my_comm_time += comm_time;
  this->critical_path_comm_time += comm_time;
  this->my_bytes+=nbytes;
  this->critical_path_bytes+=nbytes;
  this->critical_path_msg += dcost.first;
  this->my_msg += dcost.first;
  this->critical_path_wrd += dcost.second;
  this->my_wrd += dcost.second;
  costs[0] += nbytes;
  costs[1] += comm_time;
  costs[3] += dcost.second;
  costs[4] += dcost.first;
  costs[5] += comp_time;
  costs[6] += comp_time+comm_time;
  costs[7] += nbytes;
  costs[8] += comm_time;
  costs[10] += dcost.second;
  costs[11] += dcost.first;
  costs[12] += comp_time;
  costs[13] += comp_time+comm_time;

  internal_comm_info.erase(*request);
  free(comm_message_it->second);
  internal_comm_message.erase(*request);
  internal_comm_data.erase(*request);
  internal_comm_track.erase(*request);

  this->last_start_time = MPI_Wtime();
  computation_timer = this->last_start_time;
}

void tracker::get_my_data(double* container){
  container[0] = this->my_bytes;
  container[1] = this->my_comm_time;
  container[2] = this->my_bar_time;
  container[3] = this->my_msg;
  container[4] = this->my_wrd;
}

void tracker::get_critical_path_data(double* container){
  container[0] = this->critical_path_bytes;
  container[1] = this->critical_path_comm_time;
  container[2] = this->critical_path_bar_time;
  container[3] = this->critical_path_msg;
  container[4] = this->critical_path_wrd;
}

void tracker::set_critical_path_data(double* container){
  this->critical_path_bytes     = container[0];
  this->critical_path_comm_time = container[1];
  this->critical_path_bar_time  = container[2];
  this->critical_path_msg       = container[3];
  this->critical_path_wrd       = container[4];
}

void tracker::get_volume_data(double* container){
  container[0] = this->volume_bytes;
  container[1] = this->volume_comm_time;
  container[2] = this->volume_bar_time;
  container[3] = this->volume_msg;
  container[4] = this->volume_wrd;
}

void tracker::set_volume_data(double* container){
  this->volume_bytes     = container[0];
  this->volume_comm_time = container[1];
  this->volume_bar_time  = container[2];
  this->volume_msg       = container[3];
  this->volume_wrd       = container[4];
}

void tracker::set_costs(){
  // This branch ensures that we produce data only for the MPI routines actually called over the course of the program
  if (this->critical_path_bytes != 0.){
    save_info[this->name] = std::make_tuple(this->critical_path_bytes, this->critical_path_comm_time, this->critical_path_bar_time, this->critical_path_msg, this->critical_path_wrd,
                                            this->volume_bytes,        this->volume_comm_time,        this->volume_bar_time,        this->volume_msg,        this->volume_wrd);
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

// Only used with synchronous communication protocol
void propagate_critical_path(MPI_Comm cm, int nbr_pe, int nbr_pe2){
  // First exchange the tracked routine critical path data
  constexpr auto num_costs = 5*list_size;
  for (int i=0; i<list_size; i++){
    list[i]->get_critical_path_data(&old_cs[5*i]);
  }
  if (nbr_pe == -1)
    PMPI_Allreduce(&old_cs[0], &new_cs[0], num_costs, MPI_DOUBLE, MPI_MAX, cm);
  else {
    PMPI_Sendrecv(&old_cs[0], num_costs, MPI_DOUBLE, nbr_pe, internal_tag, &new_cs[0], num_costs, MPI_DOUBLE, nbr_pe, internal_tag, cm, MPI_STATUS_IGNORE);
    for (int i=0; i<num_costs; i++){
      new_cs[i] = std::max(old_cs[i], new_cs[i]);
    }
    if (nbr_pe2 != -1 && nbr_pe2 != nbr_pe){
      PMPI_Sendrecv(&new_cs[0], num_costs, MPI_DOUBLE, nbr_pe2, internal_tag, &old_cs[0], num_costs, MPI_DOUBLE, nbr_pe2, internal_tag, cm, MPI_STATUS_IGNORE);
      for (int i=0; i<num_costs; i++){
        new_cs[i] = std::max(old_cs[i], new_cs[i]);
      }
    }
  }
  for (int i=0; i<list_size; i++){
    list[i]->set_critical_path_data(&new_cs[5*i]);
  }

  // Next, exchange the critical path metric, together with tracking the rank of the process that determines each critical path
  int rank; MPI_Comm_rank(cm,&rank);
  for (int i=0; i<7; i++){
    old_cp[i].first = costs[i];
    old_cp[i].second = rank;
  }
  if (nbr_pe == -1)
    PMPI_Allreduce(old_cp, new_cp, 7, MPI_DOUBLE_INT, MPI_MAXLOC, cm);
  else {
    PMPI_Sendrecv(old_cp, 7, MPI_DOUBLE_INT, nbr_pe, internal_tag, new_cp, 7, MPI_DOUBLE_INT, nbr_pe, internal_tag, cm, MPI_STATUS_IGNORE);
    for (int i=0; i<7; i++){
      new_cp[i].first = std::max(old_cp[i].first, new_cp[i].first);
      if (old_cp[i].first<new_cp[i].first){new_cp[i].second = nbr_pe;}
    }
    if (nbr_pe2 != -1 && nbr_pe2 != nbr_pe){
      PMPI_Sendrecv(new_cp, 7, MPI_DOUBLE_INT, nbr_pe2, internal_tag, old_cp, 7, MPI_DOUBLE_INT, nbr_pe2, internal_tag, cm, MPI_STATUS_IGNORE);
      for (int i=0; i<7; i++){
        new_cp[i].first = std::max(old_cp[i].first, new_cp[i].first);
        if (old_cp[i].first<new_cp[i].first){new_cp[i].second = nbr_pe2;}
      }
    }
  }
  for (int i=0; i<7; i++){
    costs[i] = new_cp[i].first;
    root_array[i] = new_cp[i].second;
  }
  if (internal::flag){
    /* TODO: Fix later: notice that the MPI routines below use the communicator, but with send/recv, this is the wrong thing to do
    for (int i=0; i<7; i++){
      if (rank==root_array[i]){
        crit_path_size_array[i] = CritterPaths[i].size();
      } else{ crit_path_size_array[i]=0; }
    }
    // Note that instead of using a MPI_Bcast, we can use an AllReduce and that will require fewer messages (only 1)
    // set up new vectors to handle whats about to come
    PMPI_Allreduce(MPI_IN_PLACE,&crit_path_size_array[0],7,MPI_INT,MPI_MAX,cm);
    int crit_length=0;
    for (int i=0; i<7; i++){
      crit_length+=crit_path_size_array[i];
    }
    std::vector<int_int_double> crit_buffer(crit_length);
    int offset=0;
    for (int i=0; i<7; i++){
      if (rank==root_array[i]){
        assert(CritterPaths[i].size() == crit_path_size_array[i]);
        for (auto j=0; j<crit_path_size_array[i]; j++){
          crit_buffer[offset+j] = CritterPaths[i][j];
        }
      } else{
        for (auto j=0; j<crit_path_size_array[i]; j++){
          crit_buffer[offset+j] = int_int_double(0,0,0.);
        }
      }
      offset+=crit_path_size_array[i];
    }
    std::vector<int> block(2); std::vector<MPI_Aint> disp(2); std::vector<MPI_Datatype> type(2);
    MPI_Datatype MPI_INT_INT_DOUBLE;
    block[0] = 2;
    block[1] = 1;
    disp[0] = 0;
    disp[1] = 2*sizeof(int);
    type[0] = MPI_INT;
    type[1] = MPI_DOUBLE;
    MPI_Type_create_struct(2,&block[0],&disp[0],&type[0], &MPI_INT_INT_DOUBLE);
    MPI_Type_commit(&MPI_INT_INT_DOUBLE);
    MPI_Op op; MPI_Op_create((MPI_User_function*) add_critical_path_data,1,&op);
    PMPI_Allreduce(MPI_IN_PLACE,&crit_buffer[0],crit_length,MPI_INT_INT_DOUBLE,op,cm);
    MPI_Op_free(&op);
    // now copy into 7 different buffers and change their lengths (via some resize)
    offset=0;
    for (int i=0; i<7; i++){
      CritterPaths[i].resize(crit_path_size_array[i]);
      for (int j=0; j<crit_path_size_array[i]; j++){
        CritterPaths[i][j] = crit_buffer[offset+j];
      }
      offset+=crit_path_size_array[i];
    }*/
  }
}

// Note: this function should be called once per start/stop, else it will double count
void compute_volume_path(MPI_Comm cm){
  constexpr auto num_costs = 5*list_size;
  int rank; MPI_Comm_rank(cm,&rank);
  for (int i=0; i<7; i++){
    old_vp[i] = costs[7+i];
  }
  for (int i=0; i<list_size; i++){
    list[i]->get_my_data(&old_cs[5*i]);
  }
  PMPI_Allreduce(&old_cs[0], &new_cs[0], num_costs, MPI_DOUBLE, MPI_SUM, cm);
  PMPI_Allreduce(&old_vp[0], &new_vp[0], 7, MPI_DOUBLE, MPI_SUM, cm);
  for (int i=0; i<7; i++){
    costs[7+i] = new_vp[i];
  }
  for (int i=0; i<list_size; i++){
    list[i]->set_volume_data(&new_cs[5*i]);
  }
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
  Stream << "\tNumBytes\tCommunicationTime\tIdleTime\tEstimatedCommCost\tEstimatedSynchCost\tComputationTime\tRunTime";// critical path
  Stream << "\tNumBytes\tCommunicationTime\tIdleTime\tEstimatedCommCost\tEstimatedSynchCost\tComputationTime\tRunTime";// volume
  for (auto i=0;i<10;i++){
    for (auto& it : save_info){
     Stream << "\t" << it.first;
    }
  }
}

void record(std::ofstream& Stream){
  assert(internal_comm_info.size() == 0);
  auto np=0; MPI_Comm_size(MPI_COMM_WORLD,&np);
  if (is_world_root){
    auto inputs = parse_file_string();
    // Save the critter information before printing
    for (int i=0; i<list_size; i++){
      list[i]->set_costs();
    }
    if (is_first_iter){
      print_header(Stream,inputs.size());
      Stream << "\n";
    }
    print_inputs(Stream,np,inputs);
    Stream << "\t" << costs[0] << "\t" << costs[1] << "\t" << costs[2];
    Stream << "\t" << costs[3] << "\t" << costs[4] << "\t" << costs[5] << "\t" << costs[6];
    Stream << "\t" << costs[7] << "\t" << costs[8] << "\t" << costs[9];
    Stream << "\t" << costs[10] << "\t" << costs[11] << "\t" << costs[12] << "\t" << costs[13];
    for (auto& it : save_info){
      Stream << "\t" << std::get<0>(it.second);
    }
    for (auto& it : save_info){
      Stream << "\t" << std::get<1>(it.second);
    }
    for (auto& it : save_info){
      Stream << "\t" << std::get<2>(it.second);
    }
    for (auto& it : save_info){
      Stream << "\t" << std::get<3>(it.second);
    }
    for (auto& it : save_info){
      Stream << "\t" << std::get<4>(it.second);
    }
    for (auto& it : save_info){
      Stream << "\t" << std::get<5>(it.second);
    }
    for (auto& it : save_info){
      Stream << "\t" << std::get<6>(it.second);
    }
    for (auto& it : save_info){
      Stream << "\t" << std::get<7>(it.second);
    }
    for (auto& it : save_info){
      Stream << "\t" << std::get<8>(it.second);
    }
    for (auto& it : save_info){
      Stream << "\t" << std::get<9>(it.second);
    }
/*
    for (auto i=0; i<critical_paths.size(); i++){
      for (auto j=0; j<critical_paths[i].size(); j++){
        stream_track << critical_paths[i][j].first << " " << critical_paths[i][j].second << " " << critical_paths[i][j].third << std::endl;
      }
      stream_track << "0\n";
    }
*/
  }
}

void record(std::ostream& Stream){
  assert(internal_comm_info.size() == 0);
  if (is_world_root){
    // Save the critter information before printing
    for (int i=0; i<list_size; i++){
      list[i]->set_costs();
    }
    Stream << "\n\n";
    Stream << std::left << std::setw(25) << "Critical path:";
    Stream << std::left << std::setw(25) << "NumBytes";
    Stream << std::left << std::setw(25) << "CommTime";
    Stream << std::left << std::setw(25) << "IdleTime";
    Stream << std::left << std::setw(25) << "EstCommCost";
    Stream << std::left << std::setw(25) << "EstSynchCost";
    Stream << std::left << std::setw(25) << "CompTime";
    Stream << std::left << std::setw(25) << "RunTime";
    Stream << "\n";
    Stream << std::left << std::setw(25) << "                  ";
    Stream << std::left << std::setw(25) << costs[0];
    Stream << std::left << std::setw(25) << costs[1];
    Stream << std::left << std::setw(25) << costs[2];
    Stream << std::left << std::setw(25) << costs[3];
    Stream << std::left << std::setw(25) << costs[4];
    Stream << std::left << std::setw(25) << costs[5];
    Stream << std::left << std::setw(25) << costs[6] << "\n\n";

    Stream << std::left << std::setw(25) << "Volume:";
    Stream << std::left << std::setw(25) << "NumBytes";
    Stream << std::left << std::setw(25) << "CommTime";
    Stream << std::left << std::setw(25) << "IdleTime";
    Stream << std::left << std::setw(25) << "EstCommCost";
    Stream << std::left << std::setw(25) << "EstSynchCost";
    Stream << std::left << std::setw(25) << "CompTime";
    Stream << std::left << std::setw(25) << "RunTime";
    Stream << "\n";
    Stream << std::left << std::setw(25) << "                  ";
    Stream << std::left << std::setw(25) << costs[7];
    Stream << std::left << std::setw(25) << costs[8];
    Stream << std::left << std::setw(25) << costs[9];
    Stream << std::left << std::setw(25) << costs[10];
    Stream << std::left << std::setw(25) << costs[11];
    Stream << std::left << std::setw(25) << costs[12];
    Stream << std::left << std::setw(25) << costs[13] << "\n\n";

    Stream << std::left << std::setw(25) << "Critical path:";
    Stream << std::left << std::setw(25) << "NumBytes";
    Stream << std::left << std::setw(25) << "CommTime";
    Stream << std::left << std::setw(25) << "IdleTime";
    Stream << std::left << std::setw(25) << "EstSynchCost";
    Stream << std::left << std::setw(25) << "EstCommCost";
    for (auto& it : save_info){
      Stream << "\n";
      Stream << std::left << std::setw(25) << it.first;
      Stream << std::left << std::setw(25) << std::get<0>(it.second);
      Stream << std::left << std::setw(25) << std::get<1>(it.second);
      Stream << std::left << std::setw(25) << std::get<2>(it.second);
      Stream << std::left << std::setw(25) << std::get<3>(it.second);
      Stream << std::left << std::setw(25) << std::get<4>(it.second);
    }
    Stream << "\n\n";
    Stream << std::left << std::setw(25) << "Volume:";
    Stream << std::left << std::setw(25) << "NumBytes";
    Stream << std::left << std::setw(25) << "CommTime";
    Stream << std::left << std::setw(25) << "IdleTime";
    Stream << std::left << std::setw(25) << "EstSynchCost";
    Stream << std::left << std::setw(25) << "EstCommCost";
    for (auto& it : save_info){
      Stream << "\n";
      Stream << std::left << std::setw(25) << it.first;
      Stream << std::left << std::setw(25) << std::get<5>(it.second);
      Stream << std::left << std::setw(25) << std::get<6>(it.second);
      Stream << std::left << std::setw(25) << std::get<7>(it.second);
      Stream << std::left << std::setw(25) << std::get<8>(it.second);
      Stream << std::left << std::setw(25) << std::get<9>(it.second);
    }
    Stream << "\n";
  }
}
};

void print(size_t num_data, double* data){
  assert(internal::internal_comm_info.size() == 0);
  if (internal::need_new_line){
    if (internal::flag) {internal::stream << "\n";} else {std::cout << "\n";}
    internal::need_new_line=false;
  }
  if (internal::is_world_root){
    for (auto i=0; i<num_data; i++){
      if (internal::flag) {internal::stream << "\t" << data[i];} else {std::cout << "\t" << data[i];}
    }
  }
  internal::need_new_line=true;
}

void start(){
  assert(internal::internal_comm_info.size() == 0);
  internal::track=true;
  for (int i=0; i<internal::list_size; i++){
    internal::list[i]->init();
  }
  if (internal::is_world_root){
    if (!internal::is_first_iter){
      if (internal::flag) {internal::stream << "\n";} else {std::cout << "\n";}
    }
  }
  for (auto i=0; i<internal::costs.size(); i++){
    internal::costs[i]=0.;
  }
  /*Initiate new timer*/
  internal::computation_timer=MPI_Wtime();
}

void stop(){
  auto last_time = MPI_Wtime();
  assert(internal::internal_comm_info.size() == 0);
  internal::costs[5]+=(last_time-internal::computation_timer);	// update critical path computation time
  internal::costs[6]+=(last_time-internal::computation_timer);	// update critical path runtime
  internal::costs[12]+=(last_time-internal::computation_timer);	// update computation time volume
  internal::costs[13]+=(last_time-internal::computation_timer);	// update runtime volume
  internal::propagate_critical_path(MPI_COMM_WORLD,-1,-1);
  internal::compute_volume_path(MPI_COMM_WORLD);
  if (internal::flag) {internal::record(internal::stream);} else {internal::record(std::cout);}
  internal::is_first_iter = false;\
  internal::track=false;
  internal::save_info.clear();
  for (auto i=0; i<internal::costs.size(); i++){
    internal::costs[i]=0.;
  }
  internal::need_new_line=false;
/*
  for (auto i=0; i<internal::critical_paths.size(); i++){
    internal::critical_paths[i].clear();
  }
*/
}
};
