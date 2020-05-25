#include "../src/critter.h"
#include <random>

int main(int argc, char ** argv){
  size_t msg_size             = atoi(argv[1]);
  size_t sub_comm_size_factor = atoi(argv[2]);
  size_t id                   = atoi(argv[3]);
  size_t num_iter             = atoi(argv[4]);
  MPI_Init(&argc, &argv);
  int rank, size, color;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  srand48(rank);
  double* buf = (double*)malloc(msg_size*sizeof(double));	// big enough for root process to scatter
  int* recvcounts = (int*)malloc(sub_comm_size_factor*sizeof(int));
  for (int i=0; i<sub_comm_size_factor; i++) { recvcounts[i] = msg_size/sub_comm_size_factor; }
  /*for (int j=0; j<msg_size; j++){ buf[j] = drand48(); }*/
  MPI_Comm sub_comm;
  if (id%2==0) color = rank/sub_comm_size_factor;  // "local" contiguous-rank collective. When sub_comm_size_factor==ppn, each node has its own sub-communicator;
  if (id%2==1) color = color = rank%sub_comm_size_factor;  // strided-rank collectives
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &sub_comm);
  for (auto i=0; i<num_iter; i++){
    critter::start();
    if (id/2==0) MPI_Reduce_scatter(MPI_IN_PLACE, buf, recvcounts, MPI_DOUBLE, MPI_SUM, sub_comm);
    else{
      if (color==0) MPI_Reduce_scatter(MPI_IN_PLACE, buf, recvcounts, MPI_DOUBLE, MPI_SUM, sub_comm);
    }
    critter::stop();
  }
  MPI_Comm_free(&sub_comm);
  free(buf); free(recvcounts);
  MPI_Finalize();
  return 0;
}