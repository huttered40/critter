#include "mpi.h"
#include "critter.h"

void Critter::init(){
  this->last_start_time = -1.;
  this->my_bytes        = 0.;
  this->my_comm_time    = 0.;
  this->my_bar_time     = 0.;
  this->crit_bytes      = 0.;
  this->crit_comm_time  = 0.;
  this->crit_bar_time   = 0.;
}

Critter::Critter(char const * name_){
  this->name = (char*)malloc(strlen(name_)+1);
  strcpy(this->name, name_);
  this->init();
}

Critter::Critter(Critter const & t){
  this->name = (char*)malloc(strlen(t.name)+1);
  strcpy(this->name, t.name);
  this->init();
}

Critter::~Critter(){
  free(this->name);
}

void Critter::start(int64_t nbytes, MPI_Comm cm, int nbr_pe, int nbr_pe2){
  assert(this->last_start_time == -1.); //assert timer was not started twice without first being stopped
  this->last_cm = cm;
  this->last_nbr_pe = nbr_pe;
  this->last_bytes = nbytes;

  double init_time = MPI_Wtime();
  if (nbr_pe != -1)
    PMPI_Barrier(cm);
  else {
    double sbuf, rbuf;
    sbuf = 0.;
    PMPI_Sendrecv(&buf, 1, MPI_DOUBLE, nbr_pe, 1232137, &rbuf, 1, MPI_DOUBLE, nbr_pe2, 1232137, cm, MPI_STATUS_IGNORE);
  }
  this->last_start_time = MPI_Wtime();
  this->my_bar_time += this->last_start_time - init_time;
  this->crit_bar_time += this->last_start_time - init_time;
}

void Critter::stop(){
  assert(this->last_start_time != -1.); //assert timer was started 
  double dt = MPI_Wtime() - this->last_start_time;
  this->my_comm_time += dt;
  this->crit_comm_time += dt;
  this->my_bytes += this->last_bytes;
  this->crit_bytes += this->last_bytes;
}

void Critter::compute_max_crit(MPI_Comm cm){
  double old_crit_bytes = this->crit_bytes;
  PMPI_Allreduce(&old_crit_bytes, &this->crit_bytes, 1, MPI_DOUBLE, MPI_MAX, cm);
}
