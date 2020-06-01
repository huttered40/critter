#include "critter.hpp"
#include <random>

int main(int argc, char ** argv){
  size_t msg_size             = atoi(argv[1]);
  size_t sub_comm_size_factor = atoi(argv[2]);
  size_t id                   = atoi(argv[3]);
  size_t num_iter             = atoi(argv[4]);
  MPI_Init(&argc, &argv);
  int rank, world_size, color;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  srand48(rank);
  double* buf = (double*)malloc(msg_size*sizeof(double));
  for (int j=0; j<msg_size; j++){ buf[j] = drand48(); }
  MPI_Comm sub_comm;
  if (id%2==0) color = rank/sub_comm_size_factor;  // "local" contiguous-rank collective. When sub_comm_size_factor==ppn, each node has its own sub-communicator;
  if (id%2==1) color = color = rank%sub_comm_size_factor;  // strided-rank collectives
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &sub_comm);
  int root = 0;
  for (auto i=0; i<num_iter; i++){
    critter::start();
    if (id/2==0) MPI_Bcast(buf, msg_size, MPI_DOUBLE, root, sub_comm);
    else{
      if (color==0) MPI_Bcast(buf, msg_size, MPI_DOUBLE, root, sub_comm);
    }
    critter::stop();
  }
  MPI_Comm_free(&sub_comm);
  free(buf);
  MPI_Finalize();
  return 0;
}
