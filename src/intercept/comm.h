#ifndef CRITTER__INTERCEPT__COMM_H_
#define CRITTER__INTERCEPT__COMM_H_

#include <mpi.h>
#include <vector>
#include <string>

namespace critter{

void init(std::vector<std::string>& symbols);
void start(bool schedule_kernels_override = true, bool force_steady_statistical_data_overide = false);
void stop();
void record(int variantID=-1, int print_mode=1, double overhead_time=0.);
void clear(int mode=0, int tag_count=0, int* distribution_tags = nullptr);

void set_mechanism(int input_mechanism=-1);
void set_mode(int input_mode=0);
void set_debug(int debug_mode);

int get_critical_path_costs();
void get_critical_path_costs(float* costs);
int get_max_per_process_costs();
void get_max_per_process_costs(float* costs);
int get_volumetric_costs();
void get_volumetric_costs(float* costs);

namespace internal{

int init(int* argc, char*** argv);
int init_thread(int* argc, char*** argv, int required, int* provided);
int finalize();

int comm_split(MPI_Comm comm, int color, int key, MPI_Comm* newcomm);
int comm_dup(MPI_Comm comm, MPI_Comm* newcomm);
int comm_free(MPI_Comm* comm);
int get_count(MPI_Status* status, MPI_Datatype, int* count);

int barrier(MPI_Comm comm);
int bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
int reduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
int allreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int allgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int scatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int reduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int alltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int gatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int *recvcounts, const int *displs,
             MPI_Datatype recvtype, int root, MPI_Comm comm);
int allgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int *recvcounts, const int *displs,
             MPI_Datatype recvtype, MPI_Comm comm);
int scatterv(const void* sendbuf, const int *sendcounts, const int *displs, MPI_Datatype sendtype,
              void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int alltoallv(const void* sendbuf, const int *sendcounts, const int *sdispls, MPI_Datatype sendtype, void* recvbuf,
               const int *recvcounts, const int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);
int sendrecv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void* recvbuf, int recvcount,
              MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status* status);
int sendrecv_replace(void* buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                      MPI_Comm comm, MPI_Status* status);
int ssend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int send(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status);
int isend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request);
int irecv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request* request);
int ibcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request* request);
int iallreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                MPI_Request *request);
int ireduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request* request);
int igather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
             int root, MPI_Comm comm, MPI_Request* request);
int igatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int recvcounts[], const int displs[],
              MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
int iallgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
                MPI_Comm comm, MPI_Request* request);
int iallgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
                MPI_Comm comm, MPI_Request* request);
int iallgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int recvcounts[], const int displs[],
                 MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
int iscatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,
              MPI_Comm comm, MPI_Request* request);
int iscatterv(const void* sendbuf, const int sendcounts[], const int displs[], MPI_Datatype sendtype, void* recvbuf, int recvcount,
               MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request* request);
int ireduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op,
                     MPI_Comm comm, MPI_Request* request);
int ialltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
               MPI_Comm comm, MPI_Request* request);
int ialltoallv(const void* sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void* recvbuf,
                const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request);
int wait(MPI_Request* request, MPI_Status* status);
int waitany(int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status);
int waitsome(int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]);
int waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);
int test(MPI_Request* request, int* flag, MPI_Status* status);
int probe(int source, int tag, MPI_Comm comm, MPI_Status* status);
int iprobe(int source, int tag, MPI_Comm comm, int* flag, MPI_Status* status);

}
}

#endif /*CRITTER__INTERCEPT__COMM_H_*/
