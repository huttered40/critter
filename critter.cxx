#include "mpi.h"
#include "critter.h"
#include <string.h>
#include <string>
#include <assert.h>
#include <stdio.h>

Critter MPI_Barrier_critter("MPI_Barrier", 
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
        MPI_Sendrecv_replace_critter("MPI_Sendrecv_replace",
                          [](int64_t n, int p){
                            return std::pair<double,double>(1,n); 
                          }); 


Critter * critter_list[NUM_CRITTERS] = {
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
        &MPI_Sendrecv_replace_critter };
std::map<MPI_Request, Critter*> critter_req;

std::map<std::string,std::tuple<double,double,double,double,double> > saveCritterInfo;

double totalCritComputationTime;
double curComputationTimer;
double totalCommunicationTime;
double totalIdleTime;

void Critter::init(){
  this->last_start_time = -1.;
  this->my_bytes        = 0.;
  this->my_comm_time    = 0.;
  this->my_bar_time     = 0.;
  this->crit_bytes      = 0.;
  this->crit_comm_time  = 0.;
  this->crit_bar_time   = 0.;
  this->crit_msg        = 0.;
  this->crit_wrd        = 0.;
}

void Critter::initSums(){
  this->my_bytesSum        = 0.;
  this->my_comm_timeSum    = 0.;
  this->my_bar_timeSum     = 0.;
  this->crit_bytesSum      = 0.;
  this->crit_comm_timeSum  = 0.;
  this->crit_bar_timeSum   = 0.;
  this->crit_msgSum        = 0.;
  this->crit_wrdSum        = 0.;
}

Critter::Critter(char const * name_, std::function< std::pair<double,double>(int64_t,int) > 
              cost_func_){
  this->cost_func = cost_func_;
  this->name = (char*)malloc(strlen(name_)+1);
  strcpy(this->name, name_);
  this->init();
  this->initSums();
}

Critter::Critter(Critter const & t){
  this->cost_func = t.cost_func;
  this->name = (char*)malloc(strlen(t.name)+1);
  strcpy(this->name, t.name);
  this->init();
}

Critter::~Critter(){
  free(this->name);
}

void Critter::start(int64_t nelem, MPI_Datatype t, MPI_Comm cm, int nbr_pe, int nbr_pe2, bool is_async){
  //assert(this->last_start_time == -1.); //assert timer was not started twice without first being stopped
  
  // Deal with computational cost at the beginning
  volatile double curTime = MPI_Wtime();
  double computationalTimeDiff = curTime - curComputationTimer;
  if (!is_async){
    if (nbr_pe == -1){
      double curMaxTime;
      PMPI_Allreduce(&computationalTimeDiff, &curMaxTime, 1, MPI_DOUBLE, MPI_MAX, cm);
      totalCritComputationTime += curMaxTime;
    }
    else {
      double partnerCompTime;
      PMPI_Sendrecv(&computationalTimeDiff, 1, MPI_DOUBLE, nbr_pe, 1232137, &partnerCompTime, 1, MPI_DOUBLE, nbr_pe2, 1232137, cm, MPI_STATUS_IGNORE);
      totalCritComputationTime += std::max(computationalTimeDiff,partnerCompTime);
    }
  }

  // debugging check
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
  {
    // Note: this is just a sanity check. This might not be along the critical path if we have a p2p routine such as MPI_Sendrecv_replace.
    //   In such a case, each pair of processes saves the max of themselves, and then this requires a final MPI_Allreduce over total
    //     computation time at the very end in Critter_Print function.
    printf("Updated total computational time before routine: %s : %g\n", this->name, totalCritComputationTime);
  }

  this->last_cm = cm;
  this->last_nbr_pe = nbr_pe;
  this->last_nbr_pe2 = nbr_pe2;
  int el_size;
  MPI_Type_size(t, &el_size);
  int64_t nbytes = el_size * nelem;
  this->my_bytes += nbytes;
  this->crit_bytes += nbytes;
  int p;
  MPI_Comm_size(cm, &p);
  std::pair<double,double> dcost = cost_func(nbytes, p);
  this->crit_msg += dcost.first;
  this->crit_wrd += dcost.second;

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
  double idleTime;
  if (!is_async){
    if (nbr_pe == -1){
      PMPI_Allreduce(&localBarrierTime, &idleTime, 1, MPI_DOUBLE, MPI_MIN, cm);
    }
    else {
      double partnerIdleTime;
      PMPI_Sendrecv(&localBarrierTime, 1, MPI_DOUBLE, nbr_pe, 1232137, &partnerIdleTime, 1, MPI_DOUBLE, nbr_pe2, 1232137, cm, MPI_STATUS_IGNORE);
      idleTime = std::max(localBarrierTime,partnerIdleTime);
    }
  }

  // Note: each process will essentially have the same values for this now, so the Reduction at the end is uneccessary.
  this->my_bar_time += idleTime;		// we actually don't even use this
  this->crit_bar_time += idleTime;
}

void Critter::stop(){
  double dt = MPI_Wtime() - this->last_start_time;
  this->my_comm_time += dt;
  this->crit_comm_time += dt;
  this->last_start_time = MPI_Wtime();
  compute_all_max_crit(this->last_cm, this->last_nbr_pe, this->last_nbr_pe2);

  curComputationTimer = MPI_Wtime();		// reset this again

}

void Critter::compute_max_crit(MPI_Comm cm, int nbr_pe, int nbr_pe2){
  double old_cs[5];
  double new_cs[5];
  old_cs[0] = this->crit_bytes;
  old_cs[1] = this->crit_comm_time;
  old_cs[2] = this->crit_bar_time;
  old_cs[3] = this->crit_msg;
  old_cs[4] = this->crit_wrd;

  if (nbr_pe == -1)
    PMPI_Allreduce(old_cs, new_cs, 5, MPI_DOUBLE, MPI_MAX, cm);
  else {
    PMPI_Sendrecv(&old_cs, 5, MPI_DOUBLE, nbr_pe, 123213, &new_cs, 5, MPI_DOUBLE, nbr_pe, 123213, cm, MPI_STATUS_IGNORE);
    for (int i=0; i<5; i++){
      new_cs[i] = std::max(old_cs[i], new_cs[i]);
    }
    if (nbr_pe2 != -1 && nbr_pe2 != nbr_pe){
      PMPI_Sendrecv(&new_cs, 5, MPI_DOUBLE, nbr_pe2, 123214, &old_cs, 5, MPI_DOUBLE, nbr_pe2, 123214, cm, MPI_STATUS_IGNORE);
      for (int i=0; i<5; i++){
        new_cs[i] = std::max(old_cs[i], new_cs[i]);
      }
    }
  }
  this->crit_bytes     = new_cs[0];
  this->crit_comm_time = new_cs[1];
  this->crit_bar_time  = new_cs[2];
  this->crit_msg       = new_cs[3];
  this->crit_wrd       = new_cs[4];
}

void Critter::compute_avg_crit_update(){
  // Contribute to running totals
  this->crit_bytesSum += this->crit_bytes;
  this->crit_comm_timeSum += this->crit_comm_time;
  this->crit_bar_timeSum += this->crit_bar_time;
  this->crit_msgSum += this->crit_msg;
  this->crit_wrdSum += this->crit_wrd;
}

void Critter::print_crit(std::ofstream& fptr){
  if (this->last_start_time != -1.)
  {
    // No real reason to add an iteration number column to the first print statement, as that will be in order in the file its written to.
    // Only needed when writing to the file that gnuplot will then read.
    printf("%s\t %1.3E\t %1.3E\t %1.3E\t %1.3E\t %1.3E\n", this->name, this->crit_bytes, this->crit_comm_time, this->crit_bar_time, this->crit_msg, this->crit_wrd);
    
    // Instead of printing, as I did before (see below), I will save to a map and print out at the end of the iteration.
    //fptr << this->name << "\t" << this->crit_bytes << "\t" << this->crit_comm_time << "\t" << this->crit_bar_time << "\t" << this->crit_msg << "\t" << this->crit_wrd << std::endl;
    saveCritterInfo[this->name] = std::make_tuple(this->crit_bytes, this->crit_comm_time, this->crit_bar_time, this->crit_msg, this->crit_wrd);
  }
}

void Critter::print_crit_avg(std::ofstream& fptr, int numIter){
  if (this->last_start_time != -1.)
  {
    fptr << this->name << "\t" << this->crit_bytesSum/numIter << "\t" << this->crit_comm_timeSum/numIter << "\t" << this->crit_bar_timeSum/numIter << "\t" << this->crit_msgSum/numIter << "\t" << this->crit_wrdSum/numIter << std::endl;
  }
}

void Critter::print_local(){
  if (this->last_start_time != -1.)
    printf("loc%s\t %1.3E\t %1.3E\t %1.3E\n", this->name, this->my_bytes, this->my_comm_time, this->my_bar_time);
    //printf("Critter %s: local_bytes %1.3E local_comm_time %lf local_bar_time %lf\n", this->name, this->my_bytes, this->my_comm_time, this->my_bar_time);
}

std::pair<double,double> Critter::get_crit_cost(){
  return std::pair<double,double>(crit_msg, crit_wrd); 
}

void compute_all_max_crit(MPI_Comm cm, int nbr_pe, int nbr_pe2){
  for (int i=0; i<NUM_CRITTERS; i++){
    critter_list[i]->compute_max_crit(cm, nbr_pe, nbr_pe2);
  }
}

void compute_all_avg_crit_updates(){
  for (int i=0; i<NUM_CRITTERS; i++){
    critter_list[i]->compute_avg_crit_update();
  }
}
