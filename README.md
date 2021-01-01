
# critter
Welcome! If you are looking for a lightweight tool to analyze the critical path costs of your distributed-memory MPI program, you have come to the right place. `critter` **seeks to understand the critical paths of your MPI program**, and decomposes critical paths defined by the following metrics:

1. execution time
2. computation time
3. communication time
4. computation cost
5. synchronization cost (in the alpha-beta or BSP model)
6. communication cost (in the alpha-beta or BSP model)

For example, the communication-time critical path is the schedule path that incurs the maximum communication time. This path will not necessarily incur the maximum execution time.

`critter` also provides both **maximum-per-process** and **volumetric** times and costs of the measures above.

`critter` decomposes parallel schedule paths into contributions from MPI routines and user-defined kernels. User-defined kernels are encapsulated within preprocessor directives CRITTER_START(kernel_name) and CRITTER_STOP(kernel_name), which must be added manually inside source code.

See the lists below for an accurate depiction of our current support.

## Build and use instructions
`configure` compiler and flags in `config/config.mk` (MPI installation and C++11 are required). Run `make` in the main directory to generate the library files `./lib/libcritter.a`. Include `critter.h` in all files that use MPI in your application (i.e. replace `include mpi.h`), and link to `./lib/libcritter.a`. Shared library `./lib/libcritter.so` is currently not generated.

`critter` provides a few routines to the user: `critter::start()`, `critter::stop()`, and `critter::record()`. These create the window within which all MPI routines and computation kernels are intercepted and tracked. The window-creation routines are not strictly needed, as one can set the environment variable `CRITTER_AUTO=1` to enable `critter` to start tracking immediately within `MPI_Init` or `MPI_Init_thread`. Use `critter::record()` to print `critter`'s analysis. See the other environment variables below for all customization options.

## Environment variables
|     Env variable        |   description   |   default value   |    
| ----------------------- | ----------- | ---------- |
| CRITTER_MODE | Switch to enable `critter`; set to 0 to instead use a primitive timer with no user code interception | 1 |
| CRITTER_MECHANISM | Switch to enable parallel schedule decomposition (0), autotuning (1), or modeling (2) | 0 |
| CRITTER_AUTO | Switch to activate `critter` inside MPI initialization; prevents need for manually inserting `critter::start()` and `critter::stop()` inside user code; set to 1 to activate | 0 |
| CRITTER_TRACK_BLAS1 | Switch to intercept and track BLAS-1 computation kernels; set to 1 to enable | 0 |
| CRITTER_TRACK_BLAS2 | Switch to intercept and track BLAS-2 computation kernels; set to 1 to enable | 0 |
| CRITTER_TRACK_BLAS3 | Switch to intercept and track BLAS-3 computation kernels; set to 0 to disable | 1 |
| CRITTER_TRACK_LAPACK | Switch to intercept and track LAPACK computation kernels; set to 0 to disable | 1 |
| CRITTER_TRACK_SYNCHRONIZATION | Switch to measure synchronization time of collective communication; set to 1 to enable | 0 |
| CRITTER_EAGER_LIMIT | Specify maximum message size (in bytes) that can utilize eager protocol | 32768 |
| CRITTER_INCLUDE_BARRIER_TIME | Switch to include idle time in volumetric measurements; set to 1 to enable | 0 |
| CRITTER_COST_MODEL | Specify cost model: BSP (0) or Alpha-Beta (1) | 0 |
| CRITTER_DELETE_COMM | Switch to allow interception of MPI_Comm_free(...); set to 0 if communicators are deleted within interception window | 1 |
| CRITTER_PATH_DECOMPOSITION | Specify how critical paths are decomposed: by MPI routine (1), user-defined kernels (2), or avoid decomposition altogether (0) | 0 |
| CRITTER_PATH_SELECT | Specify which critical paths are decomposed (via a 6-digit string according to the following order:  Synchronization cost, Communication cost, Computation cost, Communication time, Computation time, Execution time); as an example, specify 000001 to decompose the execution-time critical path | 000000 |
| CRITTER_PATH_MEASURE_SELECT | Specify which metrics to track along each critical path; if CRITTER_PATH_DECOMPOSITION=1, specify via a 4-digit string according to the following order: Synchronization cost, Communication cost, Communication time, Synchronization time; if CRITTER_PATH_DECOMPOSITION=2, specify via a 7-digit string according to the following order: Synchronization cost, Communication cost, Computation cost, Communication time, Synchronization time, Computation time, Execution time; as an example, specify 0010 to measure the communication time attributed to each MPI routine, or specify 0000001 to measure the execution time attributed to each user-defined kernel | 0000 or 0000000 |

## Current support
|     MPI routine         |   tracked   |   
| ----------------------- | ----------- |
| MPI_Barrier              |   yes      |
| MPI_Bcast                |   yes      |
| MPI_Reduce               |   yes      |
| MPI_Allreduce            |   yes      |
| MPI_Gather               |   yes      |
| MPI_Gatherv              |   yes      |
| MPI_Allgather            |   yes      |
| MPI_Allgatherv           |   yes      |
| MPI_Scatter              |   yes      |
| MPI_Scatterv             |   yes      |
| MPI_Reduce_Scatter       |   yes      |
| MPI_Alltoall             |   yes      |
| MPI_Alltoallv            |   yes      |
| MPI_Ibcast               |   yes      |
| MPI_Ireduce              |   yes      |
| MPI_Iallreduce           |   yes      |
| MPI_Igather              |   yes      |
| MPI_Igatherv             |   yes      |
| MPI_Iallgather           |   yes      |
| MPI_Iallgatherv          |   yes      |
| MPI_Iscatter             |   yes      |
| MPI_Iscatterv            |   yes      |
| MPI_Ireduce_Scatter      |   yes      |
| MPI_Ialltoall            |   yes      |
| MPI_Ialltoallv           |   yes      |
| MPI_Send                 |   yes      |
| MPI_Ssend                |   yes      |
| MPI_Bsend                |   yes      |
| MPI_Rsend                |   no       |
| MPI_Isend                |   yes      |
| MPI_Issend               |   no       |
| MPI_Ibsend               |   no       |
| MPI_Irsend               |   no       |
| MPI_Recv                 |   yes      |
| MPI_Irecv                |   yes      |
| MPI_Sendrecv             |   yes      |
| MPI_Sendrecv_replace     |   yes      |
| MPI_Test                 |   no       |
| MPI_Probe                |   no       |

## Warnings
3. `critter` cannot handle usage of `MPI_ANY_SOURCE`.
4. `critter` cannot track libraries that use `MPI_THREAD_MULTIPLE`.
