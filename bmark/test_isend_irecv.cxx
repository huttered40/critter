#include "critter.hpp"
#include <random>

int main(int argc, char ** argv){
  size_t msg_size = atoi(argv[1]);
  size_t jump = atoi(argv[2]);
  size_t num_iter = atoi(argv[3]);
  MPI_Init(&argc, &argv);
  int rank, p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Request request,request2;
  MPI_Status status;
  srand48(rank);
  double* buf = (double*)malloc(msg_size*sizeof(double));
  /*for (int j=0; j<msg_size; j++){
    buf[j] = drand48();
  }*/
  MPI_Status st;
  int mod_rank = rank % jump;
  int half_jump = jump>>1;
  for (auto i=0; i<num_iter; i++){
    MPI_Status st;
    critter::start();
    if (mod_rank<half_jump){
      MPI_Isend(buf, msg_size, MPI_DOUBLE, (rank+half_jump)%p, 0, MPI_COMM_WORLD, &request);
      MPI_Wait(&request,&status);
    } else{
      MPI_Irecv(buf, msg_size, MPI_DOUBLE, (rank-half_jump)%p, 0, MPI_COMM_WORLD, &request2);
      MPI_Wait(&request2,&status);
    }
    critter::stop();
  }
  free(buf);
  MPI_Finalize();
  return 0;
}
