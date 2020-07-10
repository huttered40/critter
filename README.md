
# critter
Welcome! If you are looking for a lightweight tool to analyze the critical path costs of your distributed-memory MPI program, you have come to the right place. `critter` **seeks to understand the critical paths of your MPI program**, and decomposes critical paths defined by the following metrics:

1. execution time
2. computation time
3. synchronization time
4. communication time
5. estimated synchronization cost in the alpha-beta model
6. estimated synchronization cost in the BSP model
7. estimated communication cost in the alpha-beta model
8. estimated communication cost in the BSP model

For example, the communication-time critical path is the schedule path that incurs the maximum communication time. This path will not necessarily incur the maximum execution time.

`critter` also provides both **per-process** and **volumetric** times and costs of the measures above.

`critter` decomposes parallel schedule paths into contributions from MPI routines and user-defined kernels. User-defined kernels must be defined manually inside source code via CRITTER_START(kernel_name) and CRITTER_STOP(kernel_name).

See the lists below for an accurate depiction of our current support.

## Build and use instructions
`configure` compiler and flags in `config/config.mk` (MPI installation and C++11 are required). Run `make` in the main directory to generate the library files `./lib/libcritter.a`. Include `critter.h` in all files that use MPI in your application (i.e. replace `include mpi.h`), and link to `./lib/libcritter.a`. Shared library `./lib/libcritter.so` is currently not generated.

`critter` provides two routines to the user: `critter::start()` and `critter::stop()`. These create the window within which all MPI routines are intercepted and tracked. These routines are not strictly needed, as one can set the environment variable `CRITTER_AUTO=1` to enable `critter` to start tracking immediately within `MPI_Init` or `MPI_Init_thread`. See the other environment variables below for all customization options.

## Environment variables
|     Env variable        |   description   |   default value   |    
| ----------------------- | ----------- | ---------- |
| CRITTER_MODE            | serves as switch to enable `critter`; set to 1 to activate; set to 0 for simple timer with no user code interception          |   0       |
| CRITTER_AUTO            | activates `critter` inside MPI initialization; prevents need for manually inserting `critter::start()` and `critter::stop()` inside user code; set to 1 to activate          |   0       |
| CRITTER_SYMBOL_PATH_SELECT   | specifies which critical paths are decomposed by user-defined kernel; order: (estimated communication in BSP model, esimated communication in alpha-beta model, estimated synchronization in BSP model, estimated synchronization in alpha-beta model, communication time, synchronization time, computation time, execution time); as an example, specify 000000001 to decompose the execution-time critical path; specified string length must be 8          |   00000000       |
| CRITTER_COMM_PATH_SELECT   | specifies which critical paths are decomposed by MPI routines and computation/idle time; specify 000000001 to decompose the execution-time critical path; specified string length must be 8 |   00000000       |
| CRITTER_TRACK_COLLECTIVE   | intercepts blocking collective routines called within user library; set to 0 to disable          |   1       |
| CRITTER_TRACK_P2P   | intercepts p2p (blocking and nonblocking) routines called within user library; set to 0 to disable          |   1       |
| CRITTER_TRACK_P2P_IDLE   | enables idle time and synchronization time calculation for p2p communication; set to 0 to disable          |   1       |
| CRITTER_EAGER_P2P   | enforces buffered internal communication when propagating path data; set to 0 to enforce rendezvous protocol          |   0       |
| CRITTER_MAX_NUM_SYMBOLS   | max number of user-defined kernels set inside user library          |   15       |
| CRITTER_MAX_SYMBOL_LENGTH   | max length of any kernel name specified in user library          |   25       |

## Current support
|     MPI routine         |   tracked   |   tested   |    
| ----------------------- | ----------- | ---------- |
| MPI_Barrier              |   yes       |   yes       |
| MPI_Bcast                |   yes       |   yes      |
| MPI_Reduce               |   yes       |   yes      |
| MPI_Allreduce            |   yes       |   yes      |
| MPI_Gather               |   yes       |   yes       |
| MPI_Gatherv              |   yes       |   yes       |
| MPI_Allgather            |   yes       |   yes      |
| MPI_Allgatherv           |   yes       |   yes      |
| MPI_Scatter              |   yes       |   yes       |
| MPI_Scatterv             |   yes       |   yes       |
| MPI_Reduce_Scatter       |   yes       |   yes       |
| MPI_Alltoall             |   yes       |   yes       |
| MPI_Alltoallv            |   yes       |   yes       |
| MPI_Ibcast               |   yes       |   yes      |
| MPI_Ireduce              |   yes       |   no      |
| MPI_Iallreduce           |   yes       |   yes      |
| MPI_Igather              |   yes       |   no       |
| MPI_Igatherv             |   yes       |   no       |
| MPI_Iallgather           |   yes       |   no      |
| MPI_Iallgatherv          |   yes       |   no      |
| MPI_Iscatter             |   yes       |   no       |
| MPI_Iscatterv            |   yes       |   no       |
| MPI_Ireduce_Scatter      |   yes       |   no       |
| MPI_Ialltoall            |   yes       |   no       |
| MPI_Ialltoallv           |   yes       |   no       |
| MPI_Send                 |   yes       |   yes       |
| MPI_Ssend                 |   yes       |   yes       |
| MPI_Bsend                 |   yes       |   no       |
| MPI_Rsend                 |   no       |   no       |
| MPI_Isend                |   yes       |   yes       |
| MPI_Issend                 |   no       |   no       |
| MPI_Ibsend                 |   no       |   no       |
| MPI_Irsend                 |   no       |   no       |
| MPI_Recv                 |   yes       |   yes       |
| MPI_Irecv                |   yes       |   yes       |
| MPI_Sendrecv             |   yes       |   yes       |
| MPI_Sendrecv_replace     |   yes       |   yes       |

## Warnings
1. `critter` is currently not able to track user-defined kernels in any nonblocking collectives.
2. `critter` incurs large overhead when intercepting personalized collectives.
3. Any usage of `MPI_Waitany`, `MPI_Waitsome`, or `MPI_ANY_SOURCE` requires setting the environment variable `CRITTER_TRACK_P2P_IDLE=0`.
4. `critter` cannot track libraries that use `MPI_THREAD_MULTIPLE`.
