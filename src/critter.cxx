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
double CritterCostMetrics[6];	// NumBytes,CommTime,IdleTime,EstCommCost,EstSynchCost,CompTime,OverlapTime
// Instead of printing out each Critter for each iteration individually, I will save them for each iteration, print out the iteration, and then clear before next iteration
std::map<std::string,std::tuple<double,double,double,double,double,double,double,double>> saveCritterInfo;
std::map<std::string,std::vector<std::string>> AlgCritters;
std::string StreamName,FileName;
bool UseCritter;
std::ofstream Stream;
bool IsWorldRoot;

bool InAlgCritterList(std::string AlgName, std::string CritterName){
  for (auto Critter : AlgCritters[AlgName]){
    if (CritterName == Critter){
      return true;
    }
  }
  return false;
}

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

void _critter::initSums(){
/*
  this->my_bytesSum        = 0.;
  this->my_comm_timeSum    = 0.;
  this->my_bar_timeSum     = 0.;
  this->my_msg_Sum    = 0.;
  this->my_wrd_Sum     = 0.;
  this->crit_bytesSum      = 0.;
  this->crit_comm_timeSum  = 0.;
  this->crit_bar_timeSum   = 0.;
  this->crit_msgSum        = 0.;
  this->crit_wrdSum        = 0.;
*/
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
  this->my_bytes.push_back(nbytes);
  this->crit_bytes.push_back(nbytes);
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
  this->crit_bar_time.push_back(localBarrierTime);	// Will get updated after an AllReduce to find the current critical path
  this->my_bar_time.push_back(localBarrierTime);
  // start timer for communication routine
  this->last_start_time = MPI_Wtime();

  CritterCostMetrics[0] += nbytes;
  CritterCostMetrics[2] += localBarrierTime;
  CritterCostMetrics[3] += dcost.second;
  CritterCostMetrics[4] += dcost.first;
}

void _critter::stop(){
  double dt = MPI_Wtime() - this->last_start_time;
  this->my_comm_time += dt;
  CritterCostMetrics[1] += dt;
  this->crit_comm_time += dt;	// Will get updated after an AllReduce to find the current critical path
  this->last_start_time = MPI_Wtime();
  CritterCostMetrics[5] += this->save_comp_time;
  CritterCostMetrics[6] += this->save_comp_time+dt;
  compute_all_max_crit(this->last_cm, this->last_nbr_pe, this->last_nbr_pe2);
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

void _critter::compute_avg_crit_update(){
  /*
  // Contribute to running totals
  this->crit_bytesSum += this->crit_bytes;
  this->crit_comm_timeSum += this->crit_comm_time;
  this->crit_bar_timeSum += this->crit_bar_time;
  this->crit_msgSum += this->crit_msg;
  this->crit_wrdSum += this->crit_wrd;
  */
  // New meaning to this method: Compute the average of my_bytes, my_comm_time, and my_bar_time over all processes. Store them in same variables
  int WorldSize;
  MPI_Comm_size(MPI_COMM_WORLD, &WorldSize);
  double old_cs[3];
  double new_cs[3];
  old_cs[0] = this->my_bytes;
  old_cs[1] = this->my_comm_time;
  old_cs[2] = this->my_bar_time;
  PMPI_Allreduce(old_cs, new_cs, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  this->my_bytes     = new_cs[0] / WorldSize;
  this->my_comm_time = new_cs[1] / WorldSize;
  this->my_bar_time  = new_cs[2] / WorldSize;
}

void _critter::print_crit(std::ofstream& fptr, std::string name){
  if (this->last_start_time != -1.){
    // No real reason to add an iteration number column to the first print statement, as that will be in order in the file its written to.
    // Only needed when writing to the file that gnuplot will then read.
    //printf("%s\t %1.3E\t %1.3E\t %1.3E\t %1.3E\t %1.3E\n", this->name, this->crit_bytes, this->crit_comm_time, this->crit_bar_time, this->crit_msg, this->crit_wrd);
  }
  if ((this->last_start_time != -1.)){// || (InAlgCritterList(name,this->name))){
    // Instead of printing, as I did before (see below), I will save to a map and print out at the end of the iteration.
    //fptr << this->name << "\t" << this->crit_bytes << "\t" << this->crit_comm_time << "\t" << this->crit_bar_time << "\t" << this->crit_msg << "\t" << this->crit_wrd << std::endl;
    // Note: the last 3 variables (this->my_*) have been AllReduced and averaged already. They are giving the average.
    saveCritterInfo[this->name] = std::make_tuple(this->crit_bytes, this->crit_comm_time, this->crit_bar_time, this->crit_msg, this->crit_wrd,
                                                  this->my_bytes,this->my_comm_time,this->my_bar_time);
  }
}

void _critter::print_crit_avg(std::ofstream& fptr, int numIter){
/*
  if (this->last_start_time != -1.){
    fptr << this->name << "\t" << this->crit_bytesSum/numIter << "\t" << this->crit_comm_timeSum/numIter << "\t" << this->crit_bar_timeSum/numIter << "\t" << this->crit_msgSum/numIter << "\t" << this->crit_wrdSum/numIter << std::endl;
  }
*/
}

void _critter::print_local(){
  if (this->last_start_time != -1.){
    //printf("loc%s\t %1.3E\t %1.3E\t %1.3E\n", this->name, this->my_bytes, this->my_comm_time, this->my_bar_time);
    //printf("Critter %s: local_bytes %1.3E local_comm_time %lf local_bar_time %lf\n", this->name, this->my_bytes, this->my_comm_time, this->my_bar_time);
  }
}

std::pair<double,double> _critter::get_crit_cost(){
  return std::pair<double,double>(crit_msg, crit_wrd); 
}

void compute_all_max_crit(MPI_Comm cm, int nbr_pe, int nbr_pe2){
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

void compute_all_avg_crit_updates(){
  for (int i=0; i<NumCritters; i++){
    critter_list[i]->compute_avg_crit_update();
  }
}

void reset(){
  assert(critter_req.size() == 0);
  FillAlgCritterList();
  for (int i=0; i<NumCritters; i++){
    critter_list[i]->init();
  }
  CritterCostMetrics[0]=0.;
  CritterCostMetrics[1]=0.;
  CritterCostMetrics[2]=0.;
  CritterCostMetrics[3]=0.;
  CritterCostMetrics[4]=0.;
  CritterCostMetrics[5]=0.;
  CritterCostMetrics[6]=0.;
  /*Initiate new timer*/
  ComputationTimer=MPI_Wtime();
}

void PrintInputs(std::ofstream& Stream, int NumPEs, size_t NumInputs, const char** InputNames, size_t* Inputs){
  Stream << NumPEs;
  for (size_t idx = 0; idx < NumInputs; idx++){
    Stream << "\t" << Inputs[idx];
  }
}

void PrintHeader(std::ofstream& Stream, size_t NumInputs){
  for (size_t idx = 0; idx < (NumInputs+1); idx++){
    if (idx != 0){
      Stream << "\t";
    }
    Stream << "Input";
  }
  if (UseCritter){
    Stream << "\tNumBytes\tCommunicationTime\tIdleTime\tEstimatedCommCost\tEstimatedSynchCost\tComputationTime\tOverlapPotentalTime";
    for (auto i=0;i<6;i++){
      for (auto& it : saveCritterInfo){
       Stream << "\t" << it.first;
      }
    }
  }
}

void print(bool IsFirstIter, std::string AlgName, int NumPEs, size_t NumInputs, size_t* Inputs, const char** InputNames, size_t NumData, double* Data){
  assert(critter_req.size() == 0);
  .. will need to derive Inputs from FileName, via parsing it and saving it to string
  .. we need to find NumIter in the string, and then use that to divide into each critter's array of information.
  if (UseCritter){
    volatile double endTimer = MPI_Wtime();
    double timeDiff = endTimer - ComputationTimer;
    CritterCostMetrics[5]+=timeDiff;
    compute_all_max_crit(MPI_COMM_WORLD,-1,-1);
    compute_all_avg_crit_updates();
    if (IsWorldRoot){
      // Save the critter information before printing
      for (int i=0; i<NumCritters; i++){
        critter_list[i]->print_crit(Stream,AlgName);
      }
      if (IsFirstIter){
        PrintHeader(Stream,NumInputs);
        Stream << "\n";
      }
      PrintInputs(Stream,NumPEs,NumInputs,InputNames,Inputs);
      Stream << "\t" << CritterCostMetrics[0] << "\t" << CritterCostMetrics[1] << "\t" << CritterCostMetrics[2];
      Stream << "\t" << CritterCostMetrics[3] << "\t" << CritterCostMetrics[4] << "\t" << CritterCostMetrics[5];
      Stream << "\t" << CritterCostMetrics[1]+CritterCostMetrics[5]-CritterCostMetrics[6];
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
      Stream << "\n";
      saveCritterInfo.clear();
    }
  }
  else{
    if (IsWorldRoot){
      PrintHeader(Stream,NumInputs);
      Stream << "\n";
      PrintInputs(Stream,NumPEs,NumInputs,Inputs);
      for (auto i=0; i<NumData; i++){
        Stream << "\t" << Data[i];
      }
      Stream << "\n";
    }
  }
}

void start(){
  assert(critter_req.size() == 0);
  FillAlgCritterList();
  for (int i=0; i<NUM_CRITTERS; i++){
    critter_list[i]->init();
  }
  totalComputationTime=0;
  totalCommunicationTime=0;
  totalOverlapTime=0;
  /*Initiate new timer*/
  curComputationTimer=MPI_Wtime();
}

void stop(){
  assert(critter_req.size() == 0);
  record();	.. record will be used internally. print can be used by user.
}
};
