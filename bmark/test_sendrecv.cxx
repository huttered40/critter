#include "../src/critter.h"
#include <random>

int main(int argc, char ** argv){
  size_t msg_size = atoi(argv[1]);
  size_t sub_comm_size = atoi(argv[2]);
  size_t num_iter = atoi(argv[3]);
  MPI_Init(&argc, &argv);
  int rank, p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  srand48(rank);
  double* buf_send = (double*)malloc(msg_size*sizeof(double));
  double* buf_recv = (double*)malloc(msg_size*sizeof(double));
  /*for (int j=0; j<msg_size; j++){
    buf[j] = drand48();
  }*/
  MPI_Comm sub_comm;
  MPI_Comm_split(MPI_COMM_WORLD, rank/sub_comm_size, rank, &sub_comm);
  MPI_Status st;
  int partner=rank%sub_comm_size;
  if ((rank%sub_comm_size) < (sub_comm_size>>1)){
    partner+=(sub_comm_size>>1);
  } else{
    partner-=(sub_comm_size>>1);
  }
  for (auto i=0; i<num_iter; i++){
    critter::start();
    MPI_Sendrecv(buf_send, msg_size, MPI_DOUBLE, partner, 0, buf_recv, msg_size, MPI_DOUBLE, partner, 0, sub_comm, &st);
    critter::stop();
  }
  MPI_Comm_free(&sub_comm);
  free(buf_send);
  free(buf_recv);
  MPI_Finalize();
  return 0;
}
