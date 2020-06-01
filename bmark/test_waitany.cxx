#include "critter.hpp"
#include <stdio.h>

int main(int argc, char **argv){
  int rank, size, i, index;

  size_t id = atoi(argv[1]);

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<MPI_Request> requests(size);
  std::vector<MPI_Status> statuses(size);
  std::vector<int> buffer(size*10);

  critter::start(id);
  if (rank == 0){
    for (i=0; i<size * 10; i++){ buffer[i] = i/10; }
    for (i=0; i<size-1; i++){
      MPI_Isend(&buffer[i*10], 10, MPI_INT, i+1, 123, MPI_COMM_WORLD, &requests[i]);
    }
    for (i=0; i<size-1; i++){
      MPI_Waitany(size-1, &requests[0], &index, &statuses[0]);
    }
  }
  else{
    MPI_Recv(&buffer[0], 10, MPI_INT, 0, 123, MPI_COMM_WORLD, &statuses[0]);
  }
  critter::stop(id);

  MPI_Finalize();
  return 0;
}
