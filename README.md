
# critter
Welcome! If you are looking for a lightweight tool to analyze the critical path costs of your distributed-memory MPI program, you have come to the right place. `critter` **seeks to understand the critical paths of your MPI program**, and decomposes critical paths defined by the following metrics:

1. execution time
2. communication time
3. computation cost
4. synchronization cost (in the alpha-beta or Bulk-synchronous-parallel model)
5. communication cost (in the alpha-beta or Bulk-synchronous-parallel model)

For example, the communication-time critical path is the schedule path that incurs the maximum communication time. This path will not necessarily incur the maximum execution time.

`critter` also provides both **maximum-per-process** and **volumetric** times and costs of the measures above.

`critter` decomposes parallel schedule paths into contributions from MPI routines and user-defined kernels. User-defined kernels are encapsulated within preprocessor directives CRITTER_START(kernel_name) and CRITTER_STOP(kernel_name), which must be added manually inside source code.

See the lists below for an accurate depiction of our current support.

## Build and use instructions
Modify compiler and flags in `config/config.mk` (MPI installation and C++11 are required). Run `make` in the main directory to generate the library files `./lib/libcritter.a`. Include `critter.h` in all files that use MPI in your application (i.e. replace `include mpi.h`), and link to `./lib/libcritter.a`. Shared library `./lib/libcritter.so` is currently not generated.

`critter` provides the following C routines to the user:
1. `void critter_register_timer(const char* timer_name)`: register a user-defined timer for a kernel with name *timer_name* to fix ordering of intercepted kernels (each of which must be associated with a distinct name). If processes start and stop timers in more than one ordering, all intercepted and profiled kernels must register their timers using this routine
2. `void critter_start_timer(const char* timer_name, bool propagate_within = true, MPI_Comm cm = MPI_COMM_NULL)`: start timer for a kernel with name *timer_name*. Optionally, specify whether critical path information should be propagated within this kernel and/or whether to synchronize and propagate critical path information at the start of the kernel
3. `void critter_stop_timer(const char* timer_name, MPI_Comm cm = MPI_COMM_NULL)`: stop timer for a kernel with name *timer_name*. Optionally, specify whether to synchronize and propagate critical path information at the end of the kernel
4. `int critter_get_critical_path_costs()`: get the size of the critical path information (so that *critter_get_critical_path_costs(...)* can be invoked properly)
5. `void critter_get_critical_path_costs(float* costs)`: set critical path information to passed buffer *costs*
6. `void critter_start(MPI_Comm cm=MPI_COMM_WORLD)`: initiates the window within which all MPI routines and computation kernels are intercepted and profiled
7. `void critter_stop(MPI_Comm cm=MPI_COMM_WORLD)`: closes the window within which all MPI routines and computation kernels are intercepted and profiled
8. `void critter_record(int variantID=-1)`: print critter's analysis

Note that one can set the environment variable `CRITTER_AUTO_PROFILE=1` to enable *critter* to start profiling immediately following invocation `MPI_Init` or `MPI_Init_thread` and ending with invocation of `MPI_Finalize` (and thus avoid explicitly calling `critter_start(...)` and `critter_stop(...)`).
See the other environment variables below for all customization options.

## Environment variables
|     Env variable        |   description   |   default value   |    
| ----------------------- | ----------- | ---------- |
| CRITTER_MODE | Switch to enable `critter`; set to 0 to instead use a primitive timer with no user code interception | 1 |
| CRITTER_AUTO_PROFILE | Switch to activate `critter` inside MPI initialization; prevents need for manually inserting `critter::start()` and `critter::stop()` inside user code; set to 1 to activate | 0 |
| CRITTER_EAGER_LIMIT | Specify maximum message size (in bytes) that can utilize eager protocol | 32768 |
| CRITTER_COST_MODEL | Specify cost model: Bulk-Synchronous-Parallel (0) or Alpha-Beta (1) | 0 |
| CRITTER_PATH_PROFILE | Specify how critical paths are decomposed: by MPI routine (1), user-defined kernels (2), or avoid altogether (0) | 0 |
| CRITTER_PATH_SELECT | Specify which critical paths are decomposed (via a 5-digit string according to the following order:  Synchronization cost, Communication cost, Computation cost, Communication time, Execution time); as an example, specify 00001 to decompose the execution-time critical path | 00000 |
| CRITTER_PATH_MEASURE_SELECT | Specify which metrics to profile along each critical path (via a 5-digit string according to the following order:  Synchronization cost, Communication cost, Computation cost, Communication time, Execution time); as an example, specify 00010 to measure the communication time attributed to each MPI routine, or specify 00001 to measure the execution time attributed to each user-defined kernel | 00000 |
| CRITTER_PROFILE_EXCLUSIVE_TIME_ONLY | Specify whether to profile each kernel's exclusive time (1) and additionally inclusive time (0) | 0 |
| CRITTER_PROFILE_MAX_NUM_KERNELS | Specify maximum number of user-defined kernels to intercept, profile, and propagate during program runtime | 20 |
| CRITTER_PROFILE_P2P | Specify whether to profile point-to-point communications invoked during program runtime | 1 |
| CRITTER_PROFILE_COLLECTIVE | Specify whether to profile collective communications invoked during program runtime | 1 |
| CRITTER_PROPAGATE_P2P | Specify whether to propagate critical-path profiles during interception of point-to-point communications invoked during program runtime | 1 |
| CRITTER_PROPAGATE_COLLECTIVE | Specify whether to propagate critical-path profiles during interception of collective communications invoked during program runtime | 1 |
| CRITTER_EXECUTE_KERNELS | Specify whether to execute *intercepted* communication routines invoked during program runtime (subject to a maximum message size CRITTER_EXECUTE_KERNELS_MAX_MESSAGE_SIZE) | 1 |
| CRITTER_EXECUTE_KERNELS_MAX_MESSAGE_SIZE | Specify the maximum message size of an intercepted communication routine that should not be avoided if CRITTER_EXECUTE_KERNELS=1 | 1 |

## Current support (most of MPI-1, including non-blocking collectives)
|     MPI routine         |   profiled   |   
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
| MPI_Test                 |   yes      |
| MPI_Testany              |   yes      |
| MPI_Testsome             |   yes      |
| MPI_Testall              |   yes      |

## Warnings
1. `critter` cannot profile libraries that use any thread mechanism other than `MPI_THREAD_SINGLE`.
2. `critter` does not profile late receivers for point-to-point communication.
