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
        &MPI_Sendrecv_critter };
std::map<MPI_Request, Critter*> critter_req;


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

Critter::Critter(char const * name_, std::function< std::pair<double,double>(int64_t,int) > 
              cost_func_){
  this->cost_func = cost_func_;
  this->name = (char*)malloc(strlen(name_)+1);
  strcpy(this->name, name_);
  this->init();
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

  double init_time = MPI_Wtime();
  if (!is_async){
    if (nbr_pe != -1)
      PMPI_Barrier(cm);
    else {
      double sbuf, rbuf;
      sbuf = 0.;
      PMPI_Sendrecv(&sbuf, 1, MPI_DOUBLE, nbr_pe, 1232137, &rbuf, 1, MPI_DOUBLE, nbr_pe2, 1232137, cm, MPI_STATUS_IGNORE);
    }
  }
  this->last_start_time = MPI_Wtime();
  this->my_bar_time += this->last_start_time - init_time;
  this->crit_bar_time += this->last_start_time - init_time;
}

void Critter::stop(){
  double dt = MPI_Wtime() - this->last_start_time;
  this->my_comm_time += dt;
  this->crit_comm_time += dt;
  this->last_start_time = MPI_Wtime();
  this->compute_max_crit(this->last_cm, this->last_nbr_pe, this->last_nbr_pe2);
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

void Critter::print_crit(){
  if (this->last_start_time != -1.)
    printf("Critter %s: crit_bytes %1.3E crit_comm_time %lf crit_bar_time %lf crit_msg_cost %1.3E crit_wrd_cost %1.3E\n", this->name, this->crit_bytes, this->crit_comm_time, this->crit_bar_time, this->crit_msg, this->crit_wrd);
}

void Critter::print_local(){
  if (this->last_start_time != -1.)
    printf("Critter %s: local_bytes %1.3E local_comm_time %lf local_bar_time %lf\n", this->name, this->my_bytes, this->my_comm_time, this->my_bar_time);
}

std::pair<double,double> Critter::get_crit_cost(){
  return std::pair<double,double>(crit_msg, crit_wrd); 
}
