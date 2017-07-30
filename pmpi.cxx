#include "pmpi.h"

void PMPI_Timer::init(){
  this->last_start_time = -1.;
  this->my_bytes        = 0.;
  this->my_comm_time    = 0.;
  this->my_bar_time     = 0.;
  this->crit_bytes      = 0.;
  this->crit_comm_time  = 0.;
  this->crit_bar_time   = 0.;
}

PMPI_Timer::PMPI_Timer(char const * name_){
  this->name = (char*)malloc(strlen(name_)+1);
  strcpy(this->name, name_);
  this->init();
}

PMPI_Timer::PMPI_Timer(PMPI_Timer const & t){
  this->name = (char*)malloc(strlen(t.name)+1);
  strcpy(this->name, t.name);
  this->init();
}

PMPI_Timer::~PMPI_Timer(){
  free(this->name);
}

void PMPI_Timer::start(int64_t nbytes, MPI_Comm cm, int nbr_pe){
  assert(this->last_start_time == -1.); //assert timer was not started twice without first being stopped
  this->last_start_time = MPI_Wtime();

}

void PMPI_Timer::stop(){

}

void PMPI_Timer::compute_max_crit(MPI_Comm cm){

}
