#include "../src/critter.h"
#include <random>

int main(int argc, char ** argv){
  MPI_Init(&argc, &argv);
  int rank, p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  srand48(rank);
  for (int msg_size = 1; msg_size<=(1<<20); msg_size<<=1){
    double * buf = (double*)malloc(msg_size*sizeof(double));
    for (int j=0; j<msg_size; j++){
      buf[j] = drand48();
    }
    auto pcount=1;
    while (pcount<p){
      MPI_Comm sub_comm;
      MPI_Comm_split(MPI_COMM_WORLD, rank<pcount, rank, &sub_comm);
      critter::start();
      MPI_Allgather(buf, msg_size, MPI_DOUBLE, MPI_IN_PLACE, msg_size, MPI_DOUBLE, sub_comm);
      critter::stop();
      MPI_Comm_free(&sub_comm);
      pcount*=2;
    }
    free(buf);
  }
  MPI_Finalize();
  return 0;
}
