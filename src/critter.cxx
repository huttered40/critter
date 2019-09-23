#include "mpi.h"
#include "critter.h"
#include <string.h>
#include <string>
#include <assert.h>
#include <stdio.h>

namespace critter{

_critter MPI_Barrier_critter("MPI_Barrier", 
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),0.); 
                          }), 
        MPI_Bcast_critter("MPI_Bcast", 
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }), 

        MPI_Reduce_critter("MPI_Reduce", 
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }), 
        MPI_Allreduce_critter("MPI_Allreduce",
                          [](int64_t n, int p){
                            return std::pair<double,double>(2.*log2((double)p),2.*n); 
                          }), 
        MPI_Gather_critter("MPI_Gather",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Gatherv_critter("MPI_Gatherv",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Allgather_critter("MPI_Allgather",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Allgatherv_critter("MPI_Allgatherv",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Scatter_critter("MPI_Scatter",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Scatterv_critter("MPI_Scatterv",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Reduce_scatter_critter("MPI_Reduce_scatter",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),n); 
                          }), 
        MPI_Alltoall_critter("MPI_Alltoall",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n); 
                          }), 
        MPI_Alltoallv_critter("MPI_Alltoallv",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),log2((double)p)*n); 
                          }), 
        MPI_Send_critter("MPI_Send",
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }), 
        MPI_Recv_critter("MPI_Recv",
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }), 
        MPI_Isend_critter("MPI_Isend",
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }), 
        MPI_Irecv_critter("MPI_Irecv",
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }), 
        MPI_Sendrecv_critter("MPI_Sendrecv",
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }), 
        MPI_Comm_split_critter("MPI_Comm_split",
                          [](int64_t n, int p){
                            return std::pair<double,double>(log2((double)p),0); 
                          }), 
        MPI_Sendrecv_replace_critter("MPI_Sendrecv_replace",
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          });


_critter * critter_list[NumCritters] = {
        &MPI_Barrier_critter,
        &MPI_Bcast_critter,
        &MPI_Reduce_critter,
        &MPI_Allreduce_critter,
        &MPI_Scatter_critter,
        &MPI_Gather_critter,
        &MPI_Allgather_critter,
        &MPI_Allgatherv_critter,
        &MPI_Reduce_scatter_critter,
        &MPI_Alltoall_critter,
        &MPI_Alltoallv_critter,
        &MPI_Send_critter,
        &MPI_Recv_critter,
        &MPI_Irecv_critter,
        &MPI_Isend_critter,
        &MPI_Sendrecv_critter,
        &MPI_Comm_split_critter,
        &MPI_Sendrecv_replace_critter };
std::map<MPI_Request, _critter*> critter_req;

double ComputationTimer;
std::array<double,14> CritterCostMetrics;	// NumBytes,CommTime,IdleTime,EstCommCost,EstSynchCost,CompTime,OverlapTime
// Instead of printing out each Critter for each iteration individually, I will save them for each iteration, print out the iteration, and then clear before next iteration
std::map<std::string,std::tuple<double,double,double,double,double,double,double,double,double,double>> saveCritterInfo;
std::string StreamName,FileName;
std::ofstream Stream;
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

_critter::_critter(std::string name_, std::function< std::pair<double,double>(int64_t,int) > 
              cost_func_){
  this->cost_func = cost_func_;
  this->name = std::move(name_);
  this->init();
}

_critter::_critter(_critter const & t){
  this->cost_func = t.cost_func;
  this->name = t.name;
  this->init();
}

_critter::~_critter(){}

void _critter::start(int64_t nelem, MPI_Datatype t, MPI_Comm cm, int nbr_pe, int nbr_pe2, bool is_async){
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
  this->crit_msg += dcost.first;
  this->crit_wrd += dcost.second;
  this->my_msg += dcost.first;
  this->my_wrd += dcost.second;

  volatile double init_time = MPI_Wtime();
  if (!is_async){
    if (nbr_pe == -1)
      PMPI_Barrier(cm);
    else {
      double sbuf, rbuf;
      sbuf = 0.;
      PMPI_Sendrecv(&sbuf, 1, MPI_DOUBLE, nbr_pe, 1232137, &rbuf, 1, MPI_DOUBLE, nbr_pe2, 1232137, cm, MPI_STATUS_IGNORE);
    }
  }
  this->last_start_time = MPI_Wtime();
  double localBarrierTime = this->last_start_time - init_time;
  // crit_bar_time is a process-local data value for now, will get crittered after the communication routine is over.
  this->crit_bar_time+=localBarrierTime;	// Will get updated after an AllReduce to find the current critical path
  this->my_bar_time+=localBarrierTime;
  // start timer for communication routine
  this->last_start_time = MPI_Wtime();

  CritterCostMetrics[0] += nbytes;
  CritterCostMetrics[2] += localBarrierTime;
  CritterCostMetrics[3] += dcost.second;
  CritterCostMetrics[4] += dcost.first;
  CritterCostMetrics[7] += nbytes;
  CritterCostMetrics[9] += localBarrierTime;
  CritterCostMetrics[10] += dcost.second;
  CritterCostMetrics[11] += dcost.first;
}

void _critter::stop(){
  double dt = MPI_Wtime() - this->last_start_time;
  this->my_comm_time += dt;
  this->crit_comm_time += dt;	// Will get updated after an AllReduce to find the current critical path
  CritterCostMetrics[1] += dt;
  CritterCostMetrics[8] += dt;
  this->last_start_time = MPI_Wtime();
  CritterCostMetrics[5] += this->save_comp_time;
  CritterCostMetrics[6] += this->save_comp_time+dt;
  CritterCostMetrics[12] += this->save_comp_time;
  CritterCostMetrics[13] += this->save_comp_time+dt;
  compute_all_crit(this->last_cm, this->last_nbr_pe, this->last_nbr_pe2);
  // Just for sanity, lets have all processors start at same place
  PMPI_Barrier(MPI_COMM_WORLD);
  ComputationTimer = MPI_Wtime();		// reset this again
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

void _critter::set_avg_data(double* Container){
  this->my_bytes     = Container[0];
  this->my_comm_time = Container[1];
  this->my_bar_time  = Container[2];
  this->my_msg       = Container[3];
  this->my_wrd       = Container[4];
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
  constexpr auto NumCritMetrics = 5*NumCritters+7;
  double old_cs[NumCritMetrics];
  double new_cs[NumCritMetrics];
  for (int i=0; i<7; i++){
    old_cs[i] = CritterCostMetrics[i];
  }
  for (int i=0; i<NumCritters; i++){
    critter_list[i]->get_crit_data(&old_cs[5*i+7]);
  }
  if (nbr_pe == -1)
    PMPI_Allreduce(old_cs, new_cs, NumCritMetrics, MPI_DOUBLE, MPI_MAX, cm);
  else {
    PMPI_Sendrecv(&old_cs, NumCritMetrics, MPI_DOUBLE, nbr_pe, 123213, &new_cs, NumCritMetrics, MPI_DOUBLE, nbr_pe, 123213, cm, MPI_STATUS_IGNORE);
    for (int i=0; i<5; i++){
      new_cs[i] = std::max(old_cs[i], new_cs[i]);
    }
    if (nbr_pe2 != -1 && nbr_pe2 != nbr_pe){
      PMPI_Sendrecv(&new_cs, NumCritMetrics, MPI_DOUBLE, nbr_pe2, 123214, &old_cs, NumCritMetrics, MPI_DOUBLE, nbr_pe2, 123214, cm, MPI_STATUS_IGNORE);
      for (int i=0; i<5; i++){
        new_cs[i] = std::max(old_cs[i], new_cs[i]);
      }
    }
  }
  for (int i=0; i<7; i++){
    CritterCostMetrics[i] = new_cs[i];
  }
  for (int i=0; i<NumCritters; i++){
    critter_list[i]->set_crit_data(&new_cs[5*i+7]);
  }
}

void compute_all_avg(MPI_Comm cm){
  int CommSize; MPI_Comm_size(cm,&CommSize);
  constexpr auto NumCritMetrics = 5*NumCritters+7;
  double old_cs[NumCritMetrics];
  double new_cs[NumCritMetrics];
  for (int i=0; i<7; i++){
    old_cs[i] = CritterCostMetrics[7+i];
  }
  for (int i=0; i<NumCritters; i++){
    critter_list[i]->get_avg_data(&old_cs[5*i+7]);
  }
  PMPI_Allreduce(old_cs, new_cs, NumCritMetrics, MPI_DOUBLE, MPI_SUM, cm);
  for (int i=0; i<7; i++){
    CritterCostMetrics[7+i] = new_cs[i] / CommSize;
  }
  for (int i=0; i<NumCritters; i++){
    new_cs[5*i+7] /= CommSize;
    critter_list[i]->set_avg_data(&new_cs[5*i+7]);
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
  Stream << "\tNumBytes\tCommunicationTime\tIdleTime\tEstimatedCommCost\tEstimatedSynchCost\tComputationTime\tOverlapPotentalTime";// critical path
  Stream << "\tNumBytes\tCommunicationTime\tIdleTime\tEstimatedCommCost\tEstimatedSynchCost\tComputationTime\tOverlapPotentalTime";// average (per-process)
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
    Stream << "\t" << CritterCostMetrics[1]+CritterCostMetrics[5]-CritterCostMetrics[6];
    Stream << "\t" << CritterCostMetrics[7] << "\t" << CritterCostMetrics[8] << "\t" << CritterCostMetrics[9];
    Stream << "\t" << CritterCostMetrics[10] << "\t" << CritterCostMetrics[11] << "\t" << CritterCostMetrics[12];
    Stream << "\t" << CritterCostMetrics[8]+CritterCostMetrics[12]-CritterCostMetrics[13];
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
    if (IsFirstIter){
      Stream << "Critical path:\t\t\tNumBytes\t\t\tCommTime\t\t\tIdleTime\t\t\tEstCommCost\t\t\tEstSynchCost\t\t\tCompTime\t\t\tOverlapTime\n";
    }
    Stream << "\t\t\t" << CritterCostMetrics[0] << "\t\t\t" << CritterCostMetrics[1] << "\t\t\t" << CritterCostMetrics[2];
    Stream << "\t\t\t" << CritterCostMetrics[3] << "\t\t\t" << CritterCostMetrics[4] << "\t\t\t" << CritterCostMetrics[5];
    Stream << "\t\t\t" << CritterCostMetrics[1]+CritterCostMetrics[5]-CritterCostMetrics[6] << "\n\n";

    if (IsFirstIter){
      Stream << "Avg (per-process):\t\t\tNumBytes\t\t\tCommTime\t\t\tIdleTime\t\t\tEstCommCost\t\t\tEstSynchCost\t\t\tCompTime\t\t\tOverlapTime\n";
    }
    Stream << "\t\t\t" << CritterCostMetrics[7] << "\t\t\t" << CritterCostMetrics[8] << "\t\t\t" << CritterCostMetrics[9];
    Stream << "\t\t\t" << CritterCostMetrics[10] << "\t\t\t" << CritterCostMetrics[11] << "\t\t\t" << CritterCostMetrics[12];
    Stream << "\t\t\t" << CritterCostMetrics[8]+CritterCostMetrics[12]-CritterCostMetrics[13] << "\n\n";

    Stream << "Critical path:\t\t\tNumBytes\t\t\tCommTime\t\t\tIdleTime\t\t\tEstSynchCost\t\t\tEstCommCost\n";
    for (auto& it : saveCritterInfo){
      Stream << it.first;
      Stream << "\t\t\t" << std::get<0>(it.second);
      Stream << "\t\t\t" << std::get<1>(it.second);
      Stream << "\t\t\t" << std::get<2>(it.second);
      Stream << "\t\t\t" << std::get<3>(it.second);
      Stream << "\t\t\t" << std::get<4>(it.second) << "\n\n";
    }
    Stream << "Avg (per-process):\t\t\tNumBytes\t\t\tCommTime\t\t\tIdleTime\t\t\tEstSynchCost\t\t\tEstCommCost\n";
    for (auto& it : saveCritterInfo){
      Stream << it.first;
      Stream << "\t\t\t" << std::get<5>(it.second);
      Stream << "\t\t\t" << std::get<6>(it.second);
      Stream << "\t\t\t" << std::get<7>(it.second);
      Stream << "\t\t\t" << std::get<8>(it.second);
      Stream << "\t\t\t" << std::get<9>(it.second) << "\n\n";
    }
  }
}

void print(size_t NumData, double* Data){
  assert(critter_req.size() == 0);
  if (NeedNewLine){
    if (flag) {Stream << "\n";} else {std::cout << "\n";}
    NeedNewLine=false;
  }
  auto NumPEs=0; MPI_Comm_size(MPI_COMM_WORLD,&NumPEs);
  if (IsWorldRoot){
    auto Inputs = parse_file_string();
    for (auto i=0; i<NumData; i++){
      if (flag) {Stream << "\t" << Data[i];} else {std::cout << "\t" << Data[i];}
    }
  }
  NeedNewLine=true;
}

void start(){
  assert(critter_req.size() == 0);
  track=true;
  for (int i=0; i<NumCritters; i++){
    critter_list[i]->init();
  }
  if (IsWorldRoot){
    if (!IsFirstIter){
      if (flag) {Stream << "\n";} else {std::cout << "\n";}
    }
  }
  for (auto i=0; i<CritterCostMetrics.size(); i++){
    CritterCostMetrics[i]=0.;
  }
  /*Initiate new timer*/
  ComputationTimer=MPI_Wtime();
}

void stop(){
  assert(critter_req.size() == 0);
  CritterCostMetrics[5]+=(MPI_Wtime()-ComputationTimer);
  CritterCostMetrics[12]+=(MPI_Wtime()-ComputationTimer);
  compute_all_crit(MPI_COMM_WORLD,-1,-1);
  compute_all_avg(MPI_COMM_WORLD);
  if (flag) {record(Stream);} else {record(std::cout);}
  IsFirstIter = false;\
  track=false;
  saveCritterInfo.clear();
  for (auto i=0; i<CritterCostMetrics.size(); i++){
    CritterCostMetrics[i]=0.;
  }
  NeedNewLine=false;
}
};
