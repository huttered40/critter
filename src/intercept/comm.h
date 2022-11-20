#ifndef CRITTER__INTERCEPT__COMM_H_
#define CRITTER__INTERCEPT__COMM_H_

#include <mpi.h>

int critter_init(int* argc, char*** argv);
int critter_init_thread(int* argc, char*** argv, int required, int* provided);
int critter_finalize();

int critter_comm_split(MPI_Comm comm, int color, int key, MPI_Comm* newcomm);
int critter_comm_dup(MPI_Comm comm, MPI_Comm* newcomm);
int critter_comm_free(MPI_Comm* comm);

int critter_barrier(MPI_Comm comm);
int critter_bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
int critter_reduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
int critter_allreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int critter_gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int critter_allgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int critter_scatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int critter_reduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int critter_alltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int critter_gatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int *recvcounts, const int *displs,
             MPI_Datatype recvtype, int root, MPI_Comm comm);
int critter_allgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int *recvcounts, const int *displs,
             MPI_Datatype recvtype, MPI_Comm comm);
int critter_scatterv(const void* sendbuf, const int *sendcounts, const int *displs, MPI_Datatype sendtype,
              void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int critter_alltoallv(const void* sendbuf, const int *sendcounts, const int *sdispls, MPI_Datatype sendtype, void* recvbuf,
               const int *recvcounts, const int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);
int critter_sendrecv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void* recvbuf, int recvcount,
              MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status* status);
int critter_sendrecv_replace(void* buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                      MPI_Comm comm, MPI_Status* status);
int critter_ssend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int critter_send(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int critter_recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status);
int critter_isend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request);
int critter_irecv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request* request);
int critter_ibcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request* request);
int critter_iallreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                MPI_Request *request);
int critter_ireduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request* request);
int critter_igather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
             int root, MPI_Comm comm, MPI_Request* request);
int critter_igatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int recvcounts[], const int displs[],
              MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
int critter_iallgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
                MPI_Comm comm, MPI_Request* request);
int critter_iallgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
                MPI_Comm comm, MPI_Request* request);
int critter_iallgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int recvcounts[], const int displs[],
                 MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
int critter_iscatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,
              MPI_Comm comm, MPI_Request* request);
int critter_iscatterv(const void* sendbuf, const int sendcounts[], const int displs[], MPI_Datatype sendtype, void* recvbuf, int recvcount,
               MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request* request);
int critter_ireduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op,
                     MPI_Comm comm, MPI_Request* request);
int critter_ialltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
               MPI_Comm comm, MPI_Request* request);
int critter_ialltoallv(const void* sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void* recvbuf,
                const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
int critter_wait(MPI_Request* request, MPI_Status* status);
int critter_waitany(int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status);
int critter_waitsome(int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]);
int critter_waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);
int critter_testsome(int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]);
int critter_testall(int count, MPI_Request array_of_requests[], int* flag, MPI_Status array_of_statuses[]);

#endif /*CRITTER__INTERCEPT__COMM_H_*/
