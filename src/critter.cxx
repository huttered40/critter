#include "mpi.h"
#include "critter.h"
#include <string.h>
#include <string>
#include <assert.h>
#include <stdio.h>

namespace critter{
namespace internal{

void add_critter_path_data(int_int_double* in, int_int_double* inout, int* len, MPI_Datatype* dtype){
  int_int_double* invec = in;
  int_int_double* inoutvec = inout;
  for (int i=0; i<*len; i++){
    inoutvec[i].first = std::max(inoutvec[i].first,invec[i].first);
    inoutvec[i].second = std::max(inoutvec[i].second,invec[i].second);
    inoutvec[i].third = std::max(inoutvec[i].third,invec[i].third);
  }
}

_critter MPI_Barrier_critter("MPI_Barrier",0, 
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),0.); 
                          }), 
        MPI_Bcast_critter("MPI_Bcast",1,
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }), 

        MPI_Reduce_critter("MPI_Reduce",2, 
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }), 
        MPI_Allreduce_critter("MPI_Allreduce",3,
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }), 
        MPI_Gather_critter("MPI_Gather",4,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Gatherv_critter("MPI_Gatherv",5,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Allgather_critter("MPI_Allgather",6,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Allgatherv_critter("MPI_Allgatherv",7,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Scatter_critter("MPI_Scatter",8,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Scatterv_critter("MPI_Scatterv",9,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Reduce_scatter_critter("MPI_Reduce_scatter",10,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Alltoall_critter("MPI_Alltoall",11,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n); 
                          }), 
        MPI_Alltoallv_critter("MPI_Alltoallv",12,
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n); 
                          }), 
        MPI_Send_critter("MPI_Send",13,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }), 
        MPI_Recv_critter("MPI_Recv",14,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }), 
        MPI_Isend_critter("MPI_Isend",15,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }), 
        MPI_Irecv_critter("MPI_Irecv",16,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }), 
        MPI_Sendrecv_critter("MPI_Sendrecv",17,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }), 
        MPI_Sendrecv_replace_critter("MPI_Sendrecv_replace",18,
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          });


_critter * critter_list[NumCritters] = {
        &MPI_Barrier_critter,
        &MPI_Bcast_critter,
        &MPI_Reduce_critter,
        &MPI_Allreduce_critter,
        &MPI_Gather_critter,
        &MPI_Gatherv_critter,
        &MPI_Allgather_critter,
        &MPI_Allgatherv_critter,
        &MPI_Scatter_critter,
        &MPI_Scatterv_critter,
        &MPI_Reduce_scatter_critter,
        &MPI_Alltoall_critter,
        &MPI_Alltoallv_critter,
        &MPI_Send_critter,
        &MPI_Recv_critter,
        &MPI_Irecv_critter,
        &MPI_Isend_critter,
        &MPI_Sendrecv_critter,
        &MPI_Sendrecv_replace_critter };
std::map<MPI_Request, std::tuple<_critter*,double,MPI_Comm,int,int,int64_t,int,double>> critter_req;
std::vector<std::pair<MPI_Request,typename std::map<MPI_Request,std::tuple<_critter*,double,MPI_Comm,int,int,int64_t,int,double>>::iterator>> request_save;

double ComputationTimer,OverlapTimer;
std::vector<std::vector<int_int_double>> CritterPaths(8);
std::array<double,16> CritterCostMetrics;	// NumBytes,CommTime,IdleTime,EstCommCost,EstSynchCost,CompTime,OverlapTime,RunTime
// Instead of printing out each Critter for each iteration individually, I will save them for each iteration, print out the iteration, and then clear before next iteration
std::map<std::string,std::tuple<double,double,double,double,double,double,double,double,double,double>> saveCritterInfo;
std::string StreamName,StreamTrackName,FileName;
std::ofstream Stream,StreamTrack;
bool track,flag,IsWorldRoot,IsFirstIter,NeedNewLine;

void _critter::init(){
  this->last_start_time = -1.;
  this->my_bytes        = 0.;
  this->my_comm_time    = 0.;
  this->my_bar_time     = 0.;
  this->my_msg     = 0.;
  this->my_wrd     = 0.;
  this->crit_bytes      = 0.;
  this->crit_comm_time  = 0.;
  this->crit_bar_time   = 0.;
  this->crit_msg        = 0.;
  this->crit_wrd        = 0.;
  this->save_comp_time   = 0.;
}

_critter::_critter(std::string name_, int tag, std::function< std::pair<double,double>(int64_t,int) > 
              cost_func_){
  this->cost_func = cost_func_;
  this->name = std::move(name_);
  this->tag = tag;
  this->init();
}

_critter::_critter(_critter const & t){
  this->cost_func = t.cost_func;
  this->name = t.name;
  this->tag = t.tag;
  this->init();
}

_critter::~_critter(){}

void _critter::start(int64_t nelem, MPI_Datatype t, MPI_Comm cm, int nbr_pe, int nbr_pe2){
  //assert(this->last_start_time == -1.); //assert timer was not started twice without first being stopped
  
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical-path yet or that will screw up calculation of overlap!
  volatile double curTime = MPI_Wtime();
  this->save_comp_time = curTime - ComputationTimer;

  this->last_cm = cm;
  this->last_nbr_pe = nbr_pe;
  this->last_nbr_pe2 = nbr_pe2;
  int el_size;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  this->my_bytes+=nbytes;
  this->crit_bytes+=nbytes;
  int p;
  MPI_Comm_size(cm, &p);
  std::pair<double,double> dcost = cost_func(nbytes, p);
  this->last_nbytes = nbytes;
  this->last_p = p;
  this->crit_msg += dcost.first;
  this->crit_wrd += dcost.second;
  this->my_msg += dcost.first;
  this->my_wrd += dcost.second;

  volatile double init_time = MPI_Wtime();
  // If routine is asynchronous (MPI_Isend/MPI_Irecv), no need to wait for other processes
  if (nbr_pe == -1)
    PMPI_Barrier(cm);
  else {
    double sbuf, rbuf;
    sbuf = 0.;
    PMPI_Sendrecv(&sbuf, 1, MPI_DOUBLE, nbr_pe, 1232137, &rbuf, 1, MPI_DOUBLE, nbr_pe2, 1232137, cm, MPI_STATUS_IGNORE);
  }
  this->last_start_time = MPI_Wtime();
  double localBarrierTime = this->last_start_time - init_time;
  // crit_bar_time is a process-local data value for now, will get crittered after the communication routine is over.
  this->crit_bar_time+=localBarrierTime;	// Will get updated after an AllReduce to find the current critical path
  this->my_bar_time+=localBarrierTime;

  CritterCostMetrics[0] += nbytes;
  CritterCostMetrics[2] += localBarrierTime;
  CritterCostMetrics[3] += dcost.second;
  CritterCostMetrics[4] += dcost.first;
  CritterCostMetrics[8] += nbytes;
  CritterCostMetrics[10] += localBarrierTime;
  CritterCostMetrics[11] += dcost.second;
  CritterCostMetrics[12] += dcost.first;

  // Mark the local synchronization point before exchanging with its neighbors in the communicator
  CritterPaths[0].emplace_back(int_int_double(this->tag,p,nbytes));
  CritterPaths[1].emplace_back(int_int_double(this->tag,p,nbytes));
  CritterPaths[2].emplace_back(int_int_double(this->tag,p,nbytes));
  CritterPaths[3].emplace_back(int_int_double(this->tag,p,nbytes));
  CritterPaths[4].emplace_back(int_int_double(this->tag,p,nbytes));
  CritterPaths[5].emplace_back(int_int_double(this->tag,p,nbytes));
  CritterPaths[6].emplace_back(int_int_double(this->tag,p,nbytes));
  CritterPaths[7].emplace_back(int_int_double(this->tag,p,nbytes));
  // start timer for communication routine
  this->last_start_time = MPI_Wtime();
}

void _critter::stop(){
  double dt = MPI_Wtime() - this->last_start_time;
  this->my_comm_time += dt;
  this->crit_comm_time += dt;	// Will get updated after an AllReduce to find the current critical path
  CritterCostMetrics[1] += dt;
  CritterCostMetrics[9] += dt;
  CritterCostMetrics[5] += this->save_comp_time;
  CritterCostMetrics[6] += 0;
  CritterCostMetrics[7] += this->save_comp_time+dt;
  CritterCostMetrics[13] += this->save_comp_time;
  CritterCostMetrics[14] += 0;
  CritterCostMetrics[15] += this->save_comp_time+dt;
  compute_all_crit(this->last_cm, this->last_nbr_pe, this->last_nbr_pe2);
  PMPI_Barrier(this->last_cm);
  this->last_start_time = MPI_Wtime();
  ComputationTimer = this->last_start_time;
}

void _critter::istart(int64_t nelem, MPI_Datatype t, MPI_Comm cm, int nbr_pe, int nbr_pe2){
  //assert(this->last_start_time == -1.); //assert timer was not started twice without first being stopped
  
  // Deal with computational cost at the beginning, but don't synchronize to find computation-critical-path yet or that will screw up calculation of overlap!
  volatile double curTime = MPI_Wtime();
  this->save_comp_time = curTime - ComputationTimer;

  this->last_cm = cm;
  this->last_nbr_pe = nbr_pe;
  this->last_nbr_pe2 = nbr_pe2;
  int el_size;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  int p;
  MPI_Comm_size(cm, &p);
  this->last_nbytes = nbytes;
  this->last_p = p;

  // Asynchronous routines (MPI_Isend/MPI_Irecv) don't need to wait for other processes
  // There is thus no notion of barrier/idle time with overlapping communication
  // start timer for communication routine
  this->last_start_time = MPI_Wtime();
}

void _critter::istop1(MPI_Request* req){
  double dt = MPI_Wtime() - this->last_start_time;
  critter_req[*req] = std::make_tuple(this,dt,this->last_cm,this->last_nbr_pe,this->last_nbr_pe2,this->last_nbytes,this->last_p,this->save_comp_time);
  // do not write yet to any of the CritterCostMetrics
  this->last_start_time = MPI_Wtime();
  ComputationTimer = this->last_start_time;
  if (critter_req.size() == 0){
    OverlapTimer = this->last_start_time;
  }
  // Note: overlapping progress w/r/t other requests won't be corrupted
  //       because for nonblocking communication we do not block to propogate critical path info.
}

void _critter::istop2(MPI_Request req){
  assert(critter_req.find(req) != critter_req.end());
  this->last_cm = std::get<2>(critter_req[req]);
  this->last_nbr_pe = std::get<3>(critter_req[req]);
  this->last_nbr_pe2 = std::get<4>(critter_req[req]);
  double dt = std::get<1>(critter_req[req]);
  double save_comp_time = std::get<7>(critter_req[req]);
  int64_t nbytes = std::get<5>(critter_req[req]);
  int p = std::get<6>(critter_req[req]);
  this->my_comm_time += dt;
  this->crit_comm_time += dt;
  this->my_bytes+=nbytes;
  this->crit_bytes+=nbytes;
  std::pair<double,double> dcost = cost_func(nbytes, p);
  this->crit_msg += dcost.first;
  this->crit_wrd += dcost.second;
  this->my_msg += dcost.first;
  this->my_wrd += dcost.second;
  CritterCostMetrics[0] += nbytes;
  CritterCostMetrics[1] += dt;
  CritterCostMetrics[2] += 0;
  CritterCostMetrics[3] += dcost.second;
  CritterCostMetrics[4] += dcost.first;
  CritterCostMetrics[5] += save_comp_time;
  CritterCostMetrics[6] += 0;// Overlap timer should be used here
  CritterCostMetrics[7] += save_comp_time+dt;
  CritterCostMetrics[8] += nbytes;
  CritterCostMetrics[9] += dt;
  CritterCostMetrics[10] += 0;
  CritterCostMetrics[11] += dcost.second;
  CritterCostMetrics[12] += dcost.first;
  CritterCostMetrics[13] += save_comp_time;
  CritterCostMetrics[14] += 0;// Overlap timer should be used here
  CritterCostMetrics[15] += save_comp_time+dt;
  CritterPaths[0].emplace_back(int_int_double(this->tag,p,nbytes));
  CritterPaths[1].emplace_back(int_int_double(this->tag,p,nbytes));
  CritterPaths[2].emplace_back(int_int_double(this->tag,p,nbytes));
  CritterPaths[3].emplace_back(int_int_double(this->tag,p,nbytes));
  CritterPaths[4].emplace_back(int_int_double(this->tag,p,nbytes));
  CritterPaths[5].emplace_back(int_int_double(this->tag,p,nbytes));
  CritterPaths[6].emplace_back(int_int_double(this->tag,p,nbytes));
  CritterPaths[7].emplace_back(int_int_double(this->tag,p,nbytes));
  compute_all_crit(this->last_cm, this->last_nbr_pe, this->last_nbr_pe2);
  this->last_start_time = MPI_Wtime();
  ComputationTimer = this->last_start_time;
}

void _critter::get_crit_data(double* Container){
  Container[0] = this->crit_bytes;
  Container[1] = this->crit_comm_time;
  Container[2] = this->crit_bar_time;
  Container[3] = this->crit_msg;
  Container[4] = this->crit_wrd;
}

void _critter::set_crit_data(double* Container){
  this->crit_bytes     = Container[0];
  this->crit_comm_time = Container[1];
  this->crit_bar_time  = Container[2];
  this->crit_msg       = Container[3];
  this->crit_wrd       = Container[4];
}

void _critter::get_avg_data(double* Container){
  Container[0] = this->my_bytes;
  Container[1] = this->my_comm_time;
  Container[2] = this->my_bar_time;
  Container[3] = this->my_msg;
  Container[4] = this->my_wrd;
}

void _critter::set_avg_data(double* Container, int CommSize){
  this->my_bytes     = Container[0]/CommSize;
  this->my_comm_time = Container[1]/CommSize;
  this->my_bar_time  = Container[2]/CommSize;
  this->my_msg       = Container[3]/CommSize;
  this->my_wrd       = Container[4]/CommSize;
}

void _critter::save_crit(){
  if (this->last_start_time != -1.){
    saveCritterInfo[this->name] = std::make_tuple(this->crit_bytes, this->crit_comm_time, this->crit_bar_time, this->crit_msg, this->crit_wrd,
                                                  this->my_bytes,this->my_comm_time,this->my_bar_time, this->my_msg, this->my_wrd);
  }
}

std::pair<double,double> _critter::get_crit_cost(){
  return std::pair<double,double>(crit_msg, crit_wrd); 
}

std::vector<std::string> parse_file_string(){
  std::vector<std::string> Inputs;
  auto prev=0;
  auto First=false;
  for (auto i=0; i<FileName.size(); i++){
    if (FileName[i]=='+'){
      if (First){
        Inputs.emplace_back(FileName.substr(prev,i-prev));
      }
      else{
        First=true;
      }
      prev=i+1;
    }
  }
  return Inputs;
}

void compute_all_crit(MPI_Comm cm, int nbr_pe, int nbr_pe2){
  // First exchange the tracked routine critical path data
  constexpr auto NumCritMetrics = 5*NumCritters;
  double old_cs[NumCritMetrics];
  double new_cs[NumCritMetrics];
  for (int i=0; i<NumCritters; i++){
    critter_list[i]->get_crit_data(&old_cs[5*i]);
  }
  if (nbr_pe == -1)
    PMPI_Allreduce(old_cs, new_cs, NumCritMetrics, MPI_DOUBLE, MPI_MAX, cm);
  else {
    PMPI_Sendrecv(&old_cs, NumCritMetrics, MPI_DOUBLE, nbr_pe, 123213, &new_cs, NumCritMetrics, MPI_DOUBLE, nbr_pe, 123213, cm, MPI_STATUS_IGNORE);
    for (int i=0; i<NumCritMetrics; i++){
      new_cs[i] = std::max(old_cs[i], new_cs[i]);
    }
    if (nbr_pe2 != -1 && nbr_pe2 != nbr_pe){
      PMPI_Sendrecv(&new_cs, NumCritMetrics, MPI_DOUBLE, nbr_pe2, 123214, &old_cs, NumCritMetrics, MPI_DOUBLE, nbr_pe2, 123214, cm, MPI_STATUS_IGNORE);
      for (int i=0; i<NumCritMetrics; i++){
        new_cs[i] = std::max(old_cs[i], new_cs[i]);
      }
    }
  }
  for (int i=0; i<NumCritters; i++){
    critter_list[i]->set_crit_data(&new_cs[5*i]);
  }

  if (internal::flag){
    // Next, exchange the critical path metric, together with tracking the rank of the process that determines each critical path
    int rank; MPI_Comm_rank(cm,&rank);
    double_int old_cp[8];
    double_int new_cp[8];
    int root_array[8];
    int crit_path_size_array[8];
    for (int i=0; i<8; i++){
      old_cp[i].first = CritterCostMetrics[i];
      old_cp[i].second = rank;
    }
    if (nbr_pe == -1)
      PMPI_Allreduce(old_cp, new_cp, 8, MPI_DOUBLE_INT, MPI_MAXLOC, cm);
    else {
      PMPI_Sendrecv(&old_cp, 8, MPI_DOUBLE_INT, nbr_pe, 123213, &new_cp, 8, MPI_DOUBLE_INT, nbr_pe, 123213, cm, MPI_STATUS_IGNORE);
      for (int i=0; i<8; i++){
        new_cp[i].first = std::max(old_cp[i].first, new_cp[i].first);
        if (old_cp[i].first<new_cp[i].first){new_cp[i].second = nbr_pe;}
      }
      if (nbr_pe2 != -1 && nbr_pe2 != nbr_pe){
        PMPI_Sendrecv(&new_cp, 8, MPI_DOUBLE_INT, nbr_pe2, 123214, &old_cp, 8, MPI_DOUBLE_INT, nbr_pe2, 123214, cm, MPI_STATUS_IGNORE);
        for (int i=0; i<8; i++){
          new_cp[i].first = std::max(old_cp[i].first, new_cp[i].first);
          if (old_cp[i].first<new_cp[i].first){new_cp[i].second = nbr_pe2;}
        }
      }
    }
    for (int i=0; i<8; i++){
      CritterCostMetrics[i] = new_cp[i].first;
      root_array[i] = new_cp[i].second;
    }
    /* TODO: Fix later: notice that the MPI routines below use the communicator, but with send/recv, this is the wrong thing to do
    for (int i=0; i<8; i++){
      if (rank==root_array[i]){
        crit_path_size_array[i] = CritterPaths[i].size();
      } else{ crit_path_size_array[i]=0; }
    }
    // Note that instead of using a MPI_Bcast, we can use an AllReduce and that will require fewer messages (only 1)
    // set up new vectors to handle whats about to come
    PMPI_Allreduce(MPI_IN_PLACE,&crit_path_size_array[0],8,MPI_INT,MPI_MAX,cm);
    int crit_length=0;
    for (int i=0; i<8; i++){
      crit_length+=crit_path_size_array[i];
    }
    std::vector<int_int_double> crit_buffer(crit_length);
    int offset=0;
    for (int i=0; i<8; i++){
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
    MPI_Op op; MPI_Op_create((MPI_User_function*) add_critter_path_data,1,&op);
    PMPI_Allreduce(MPI_IN_PLACE,&crit_buffer[0],crit_length,MPI_INT_INT_DOUBLE,op,cm);
    MPI_Op_free(&op);
    // now copy into 8 different buffers and change their lengths (via some resize)
    offset=0;
    for (int i=0; i<8; i++){
      CritterPaths[i].resize(crit_path_size_array[i]);
      for (int j=0; j<crit_path_size_array[i]; j++){
        CritterPaths[i][j] = crit_buffer[offset+j];
      }
      offset+=crit_path_size_array[i];
    }*/
  }
}

void compute_all_avg(MPI_Comm cm){
  int CommSize; MPI_Comm_size(cm,&CommSize);
  constexpr auto NumCritMetrics = 5*NumCritters+8;
  double old_cs[NumCritMetrics];
  double new_cs[NumCritMetrics];
  for (int i=0; i<8; i++){
    old_cs[i] = CritterCostMetrics[8+i];
  }
  for (int i=0; i<NumCritters; i++){
    critter_list[i]->get_avg_data(&old_cs[5*i+8]);
  }
  PMPI_Allreduce(old_cs, new_cs, NumCritMetrics, MPI_DOUBLE, MPI_SUM, cm);
  for (int i=0; i<8; i++){
    CritterCostMetrics[8+i] = new_cs[i] / CommSize;
  }
  for (int i=0; i<NumCritters; i++){
    critter_list[i]->set_avg_data(&new_cs[5*i+8],CommSize);
  }
}

template<typename StreamType>
void PrintInputs(StreamType& Stream, int NumPEs, std::vector<std::string>& Inputs){
  Stream << NumPEs;
  for (auto InputStr : Inputs){
    Stream << "\t" << InputStr;
  }
}

template<typename StreamType>
void PrintHeader(StreamType& Stream, size_t NumInputs){
  for (size_t idx = 0; idx < (NumInputs+1); idx++){
    if (idx != 0){
      Stream << "\t";
    }
    Stream << "Input";
  }
  Stream << "\tNumBytes\tCommunicationTime\tIdleTime\tEstimatedCommCost\tEstimatedSynchCost\tComputationTime\tOverlapPotentalTime\tRunTime";// critical path
  Stream << "\tNumBytes\tCommunicationTime\tIdleTime\tEstimatedCommCost\tEstimatedSynchCost\tComputationTime\tOverlapPotentalTime\tRunTime";// average (per-process)
  for (auto i=0;i<10;i++){
    for (auto& it : saveCritterInfo){
     Stream << "\t" << it.first;
    }
  }
}

void record(std::ofstream& Stream){
  assert(critter_req.size() == 0);
  auto NumPEs=0; MPI_Comm_size(MPI_COMM_WORLD,&NumPEs);
  if (IsWorldRoot){
    auto Inputs = parse_file_string();
    // Save the critter information before printing
    for (int i=0; i<NumCritters; i++){
      critter_list[i]->save_crit();
    }
    if (IsFirstIter){
      PrintHeader(Stream,Inputs.size());
      Stream << "\n";
    }
    PrintInputs(Stream,NumPEs,Inputs);
    Stream << "\t" << CritterCostMetrics[0] << "\t" << CritterCostMetrics[1] << "\t" << CritterCostMetrics[2];
    Stream << "\t" << CritterCostMetrics[3] << "\t" << CritterCostMetrics[4] << "\t" << CritterCostMetrics[5];
    Stream << "\t" << CritterCostMetrics[6] << "\t" << CritterCostMetrics[7];
    Stream << "\t" << CritterCostMetrics[8] << "\t" << CritterCostMetrics[9] << "\t" << CritterCostMetrics[10];
    Stream << "\t" << CritterCostMetrics[11] << "\t" << CritterCostMetrics[12] << "\t" << CritterCostMetrics[13];
    Stream << "\t" << CritterCostMetrics[14] << "\t" << CritterCostMetrics[15];
    for (auto& it : saveCritterInfo){
      Stream << "\t" << std::get<0>(it.second);
    }
    for (auto& it : saveCritterInfo){
      Stream << "\t" << std::get<1>(it.second);
    }
    for (auto& it : saveCritterInfo){
      Stream << "\t" << std::get<2>(it.second);
    }
    for (auto& it : saveCritterInfo){
      Stream << "\t" << std::get<3>(it.second);
    }
    for (auto& it : saveCritterInfo){
      Stream << "\t" << std::get<4>(it.second);
    }
    for (auto& it : saveCritterInfo){
      Stream << "\t" << std::get<5>(it.second);
    }
    for (auto& it : saveCritterInfo){
      Stream << "\t" << std::get<6>(it.second);
    }
    for (auto& it : saveCritterInfo){
      Stream << "\t" << std::get<7>(it.second);
    }
    for (auto& it : saveCritterInfo){
      Stream << "\t" << std::get<8>(it.second);
    }
    for (auto& it : saveCritterInfo){
      Stream << "\t" << std::get<9>(it.second);
    }

    for (auto i=0; i<CritterPaths.size(); i++){
      for (auto j=0; j<CritterPaths[i].size(); j++){
        StreamTrack << CritterPaths[i][j].first << " " << CritterPaths[i][j].second << " " << CritterPaths[i][j].third << std::endl;
      }
      StreamTrack << "0\n";
    }
  }
}

void record(std::ostream& Stream){
  assert(critter_req.size() == 0);
  auto NumPEs=0; MPI_Comm_size(MPI_COMM_WORLD,&NumPEs);
  if (IsWorldRoot){
    // Save the critter information before printing
    for (int i=0; i<NumCritters; i++){
      critter_list[i]->save_crit();
    }
    Stream << "\n\n";
    Stream << std::left << std::setw(25) << "Critical path:";
    Stream << std::left << std::setw(25) << "NumBytes";
    Stream << std::left << std::setw(25) << "CommTime";
    Stream << std::left << std::setw(25) << "IdleTime";
    Stream << std::left << std::setw(25) << "EstCommCost";
    Stream << std::left << std::setw(25) << "EstSynchCost";
    Stream << std::left << std::setw(25) << "CompTime";
    Stream << std::left << std::setw(25) << "OverlapTime";
    Stream << std::left << std::setw(25) << "RunTime";
    Stream << "\n";
    Stream << std::left << std::setw(25) << "                  ";
    Stream << std::left << std::setw(25) << CritterCostMetrics[0];
    Stream << std::left << std::setw(25) << CritterCostMetrics[1];
    Stream << std::left << std::setw(25) << CritterCostMetrics[2];
    Stream << std::left << std::setw(25) << CritterCostMetrics[3];
    Stream << std::left << std::setw(25) << CritterCostMetrics[4];
    Stream << std::left << std::setw(25) << CritterCostMetrics[5];
    Stream << std::left << std::setw(25) << CritterCostMetrics[6];
    Stream << std::left << std::setw(25) << CritterCostMetrics[7] << "\n\n";

    Stream << std::left << std::setw(25) << "Avg (per-process):";
    Stream << std::left << std::setw(25) << "NumBytes";
    Stream << std::left << std::setw(25) << "CommTime";
    Stream << std::left << std::setw(25) << "IdleTime";
    Stream << std::left << std::setw(25) << "EstCommCost";
    Stream << std::left << std::setw(25) << "EstSynchCost";
    Stream << std::left << std::setw(25) << "CompTime";
    Stream << std::left << std::setw(25) << "OverlapTime";
    Stream << std::left << std::setw(25) << "RunTime";
    Stream << "\n";
    Stream << std::left << std::setw(25) << "                  ";
    Stream << std::left << std::setw(25) << CritterCostMetrics[8];
    Stream << std::left << std::setw(25) << CritterCostMetrics[9];
    Stream << std::left << std::setw(25) << CritterCostMetrics[10];
    Stream << std::left << std::setw(25) << CritterCostMetrics[11];
    Stream << std::left << std::setw(25) << CritterCostMetrics[12];
    Stream << std::left << std::setw(25) << CritterCostMetrics[13];
    Stream << std::left << std::setw(25) << CritterCostMetrics[14];
    Stream << std::left << std::setw(25) << CritterCostMetrics[15] << "\n\n";

    Stream << std::left << std::setw(25) << "Critical path:";
    Stream << std::left << std::setw(25) << "NumBytes";
    Stream << std::left << std::setw(25) << "CommTime";
    Stream << std::left << std::setw(25) << "IdleTime";
    Stream << std::left << std::setw(25) << "EstSynchCost";
    Stream << std::left << std::setw(25) << "EstCommCost";
    for (auto& it : saveCritterInfo){
      Stream << "\n";
      Stream << std::left << std::setw(25) << it.first;
      Stream << std::left << std::setw(25) << std::get<0>(it.second);
      Stream << std::left << std::setw(25) << std::get<1>(it.second);
      Stream << std::left << std::setw(25) << std::get<2>(it.second);
      Stream << std::left << std::setw(25) << std::get<3>(it.second);
      Stream << std::left << std::setw(25) << std::get<4>(it.second);
    }
    Stream << "\n\n";
    Stream << std::left << std::setw(25) << "Avg (per-process):";
    Stream << std::left << std::setw(25) << "NumBytes";
    Stream << std::left << std::setw(25) << "CommTime";
    Stream << std::left << std::setw(25) << "IdleTime";
    Stream << std::left << std::setw(25) << "EstSynchCost";
    Stream << std::left << std::setw(25) << "EstCommCost";
    for (auto& it : saveCritterInfo){
      Stream << "\n";
      Stream << std::left << std::setw(25) << it.first;
      Stream << std::left << std::setw(25) << std::get<5>(it.second);
      Stream << std::left << std::setw(25) << std::get<6>(it.second);
      Stream << std::left << std::setw(25) << std::get<7>(it.second);
      Stream << std::left << std::setw(25) << std::get<8>(it.second);
      Stream << std::left << std::setw(25) << std::get<9>(it.second);
    }
  }
}
};

void print(size_t NumData, double* Data){
  assert(internal::critter_req.size() == 0);
  if (internal::NeedNewLine){
    if (internal::flag) {internal::Stream << "\n";} else {std::cout << "\n";}
    internal::NeedNewLine=false;
  }
  auto NumPEs=0; MPI_Comm_size(MPI_COMM_WORLD,&NumPEs);
  if (internal::IsWorldRoot){
    auto Inputs = internal::parse_file_string();
    for (auto i=0; i<NumData; i++){
      if (internal::flag) {internal::Stream << "\t" << Data[i];} else {std::cout << "\t" << Data[i];}
    }
  }
  internal::NeedNewLine=true;
}

void start(){
  assert(internal::critter_req.size() == 0);
  internal::track=true;
  for (int i=0; i<internal::NumCritters; i++){
    internal::critter_list[i]->init();
  }
  if (internal::IsWorldRoot){
    if (!internal::IsFirstIter){
      if (internal::flag) {internal::Stream << "\n";} else {std::cout << "\n";}
    }
  }
  for (auto i=0; i<internal::CritterCostMetrics.size(); i++){
    internal::CritterCostMetrics[i]=0.;
  }
  /*Initiate new timer*/
  internal::ComputationTimer=MPI_Wtime();
  internal::OverlapTimer=0.;
}

void stop(){
  assert(internal::critter_req.size() == 0);
  auto last_time = MPI_Wtime();
  internal::CritterCostMetrics[5]+=(last_time-internal::ComputationTimer);
  internal::CritterCostMetrics[13]+=(last_time-internal::ComputationTimer);
  internal::compute_all_crit(MPI_COMM_WORLD,-1,-1);
  internal::compute_all_avg(MPI_COMM_WORLD);
  if (internal::flag) {internal::record(internal::Stream);} else {internal::record(std::cout);}
  internal::IsFirstIter = false;\
  internal::track=false;
  internal::saveCritterInfo.clear();
  for (auto i=0; i<internal::CritterCostMetrics.size(); i++){
    internal::CritterCostMetrics[i]=0.;
  }
  internal::NeedNewLine=false;
  for (auto i=0; i<internal::CritterPaths.size(); i++){
    internal::CritterPaths[i].clear();
  }
}
};
