#include "../critter.h"
#include <random>

int main(int argc, char ** argv){
  MPI_Init(&argc, &argv);

  int rank, p;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  srand48(rank);

  for (int msg_size = 1; msg_size<=(1<<20); msg_size<<=1){
    for (int i=1; i<p; i++){
      MPI_Comm sub_comm;
      MPI_Comm_split(MPI_COMM_WORLD, rank<i, rank, &sub_comm);
      if (rank<i){
        double * buf = (double*)malloc(msg_size*sizeof(double));
        for (int j=0; j<msg_size; j++){
          buf[j] = drand48();
        }
        MPI_Bcast(buf, msg_size, MPI_DOUBLE, i-1, sub_comm);
        free(buf);
      }
      MPI_Comm_free(&sub_comm);
    }
  }

  MPI_Finalize();

  return 0;
}
