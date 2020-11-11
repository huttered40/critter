#ifndef CRITTER__INTERCEPT__COMM_H_
#define CRITTER__INTERCEPT__COMM_H_

#include "mpi.h"

namespace critter{

void start(bool schedule_kernels_override = true, bool force_steady_statistical_data_overide = false);
void stop();
void record(int variantID=-1, int print_mode=1, double overhead_time=0.);
void clear(int mode=0, int tag_count=0, int* distribution_tags = nullptr);

void set_mechanism(int input_mechanism=-1);
void set_mode(int input_mode=-1);
void set_debug(int debug_mode);

namespace internal{

void init(int* argc, char*** argv);
void init_thread(int* argc, char*** argv, int required, int* provided);
void barrier(MPI_Comm comm);
void comm_split(MPI_Comm comm, int color, int key, MPI_Comm* newcomm);
void comm_dup(MPI_Comm comm, MPI_Comm* newcomm);
void comm_free(MPI_Comm* comm);
void comm_split(MPI_Comm comm, int color, int key, MPI_Comm* new_comm);
void bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
void reduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
void allreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
void gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
void allgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
void scatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
void reduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
void alltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
void gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int *recvcounts, const int *displs,
             MPI_Datatype recvtype, int root, MPI_Comm comm);
void allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int *recvcounts, const int *displs,
             MPI_Datatype recvtype, MPI_Comm comm);
void scatterv(const void *sendbuf, const int *sendcounts, const int *displs, MPI_Datatype sendtype,
              void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
void alltoallv(const void *sendbuf, const int *sendcounts, const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
               const int *recvcounts, const int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);
void sendrecv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void* recvbuf, int recvcount,
              MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status* status);
void sendrecv_replace(void* buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                      MPI_Comm comm, MPI_Status* status);
void ssend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
void send(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
void recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status);
void isend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request);
void irecv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request* request);
void ibcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request* request);
void iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                MPI_Request *request);
void ireduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request* request);
void igather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
             int root, MPI_Comm comm, MPI_Request* request);
void igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[],
              MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
void iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
                MPI_Comm comm, MPI_Request* request);
void iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
                MPI_Comm comm, MPI_Request* request);
void iallgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int recvcounts[], const int displs[],
                 MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
void iscatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,
              MPI_Comm comm, MPI_Request* request);
void iscatterv(const void* sendbuf, const int sendcounts[], const int displs[], MPI_Datatype sendtype, void* recvbuf, int recvcount,
               MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request* request);
void ireduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op,
                     MPI_Comm comm, MPI_Request* request);
void ialltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
               MPI_Comm comm, MPI_Request* request);
void ialltoallv(const void* sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void* recvbuf,
                const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
void wait(MPI_Request* request, MPI_Status* status);
void waitany(int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status);
void waitsome(int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]);
void waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);
void finalize();

}
}

#endif /*CRITTER__INTERCEPT__COMM_H_*/
