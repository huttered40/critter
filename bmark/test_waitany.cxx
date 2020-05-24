#include "../src/critter.h"
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
    for (i=0; i<size * 10; i++)
      buffer[i] = i/10;
      for (i=0; i<size-1; i++){
        MPI_Isend(&buffer[i*10], 10, MPI_INT, i+1, 123, MPI_COMM_WORLD, &requests[i]);
      }
      for (i=0; i<size-1; i++){
        std::cout << "request before - " << requests[i] << std::endl;
        MPI_Waitany(size-1, &requests[0], &index, &statuses[0]);
        std::cout << "request after - " << requests[i] << std::endl;
        std::cout << "index - " << index << std::endl;
      }
    }
    else{
      MPI_Recv(&buffer[0], 10, MPI_INT, 0, 123, MPI_COMM_WORLD, &statuses[0]);
      //printf("%d: buffer[0] = %d\n", rank, buffer[0]);fflush(stdout);
    }
    critter::stop(id);

    MPI_Finalize();
    return 0;
}
