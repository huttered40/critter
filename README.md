
# critter
Welcome! If you are looking for a lightweight tool to analyze the critical path costs of your distributed-memory MPI program, you have come to the right place. `critter` **seeks to understand the critical paths of your MPI program**, and tracks the critical path of each of the following measurements:

1. total runtime
2. communication time
3. computation time
4. number of bytes communicated
5. estimated communication cost
6. estimated synchronization cost

`critter` also provides both **per-process** and **volumetric** times and costs of the measures above.

Critical paths through parallel schedules incur contributions from many functions; `critter` can track the contributions (in terms of the measures listed above) of each function and/or any user-defined symbol.

In addition, `critter` can break down each of the critical path measures into the contributions from each MPI routine.

See the lists below for an accurate depiction of our current support.

## Build and use instructions
`configure` compiler and flags in `config/config.mk` (MPI installation and C++11 are required). Run `make` in the main directory to generate the library files `./lib/libcritter.a` and `./lib/libcritter.so`. Include `critter.h` in all files that use MPI in your application (i.e. replace `include mpi.h`), and link to `./lib/libcritter.a`.

`critter` provides a few variables to the user inside the `critter` namespace in `src/critter.h`. Along with a few tags `internal_tag*` to prevent critter's internal MPI communication from conflicting with user communication, `critter_breakdown` specifies which, if any, critical path measurement is broken down into contributions from individual MPI routines. `critter_breakdown_size` must match `critter_breakdown.count()`. `max_timer_name_length` and `max_num_symbols` additionally serve as modifiable compile-time variables that must be large enough to support the number of tracked symbols. Understand that increased profiling information comes at a cost of increased internal data transfer necessary to propogate critical path information.

`critter` borrows from the `tau` syntax in that all instances of `TAU_START(...)` and `TAU_STOP(...)` are intercepted and tracked.

`critter` provides two routines to the user: `critter::start(int mode)` and `critter::stop(int mode)`. These create the window within which all MPI routines are intercepted and tracked. The argument `mode` can be 0, 1, or 2. Each signifies a distinct set of profiling measurements:

0. Mode 0 acts as a simple timer with no profiling output apart from a walltime.
1. Mode 1 tracks the critical path, per-process, and volumetric measurements for all six measurements listed above along any of the six critical paths. Additionally, `critter` profiles the contributions from each MPI routine.
2. Mode 2 supersedes mode 1 by tracking the contributions (in terms of the measures listed above) of each function and/or any user-defined symbol. Only symbols enclosed in `TAU_START` and `TAU_STOP` are tracked.

## Current support
|     MPI routine         |   tracked   |   tested   |    benchmarks   |     
| ----------------------- | ----------- | ---------- | --------------- |
| MPI_Barrier              |   yes       |   yes       |   no            |
| MPI_Bcast                |   yes       |   yes      |   yes           |
| MPI_Reduce               |   yes       |   yes      |   yes           |
| MPI_Allreduce            |   yes       |   yes      |   yes           |
| MPI_Gather               |   yes       |   yes       |   no            |
| MPI_Gatherv              |   yes       |   yes       |   no            |
| MPI_Allgather            |   yes       |   yes      |   yes           |
| MPI_Allgatherv           |   yes       |   yes      |   no           |
| MPI_Scatter              |   yes       |   yes       |   no            |
| MPI_Scatterv             |   yes       |   yes       |   no            |
| MPI_Reduce_Scatter       |   yes       |   yes       |   no            |
| MPI_Alltoall             |   yes       |   yes       |   no            |
| MPI_Alltoallv            |   yes       |   yes       |   no            |
| MPI_Ibcast               |   yes       |   no      |   no           |
| MPI_Ireduce              |   yes       |   no      |   no           |
| MPI_Iallreduce           |   yes       |   no      |   no           |
| MPI_Igather              |   yes       |   no       |   no            |
| MPI_Igatherv             |   yes       |   no       |   no            |
| MPI_Iallgather           |   yes       |   no      |   no           |
| MPI_Iallgatherv          |   yes       |   no      |   no           |
| MPI_Iscatter             |   yes       |   no       |   no            |
| MPI_Iscatterv            |   yes       |   no       |   no            |
| MPI_Ireduce_Scatter      |   yes       |   no       |   no            |
| MPI_Ialltoall            |   yes       |   no       |   no            |
| MPI_Ialltoallv           |   yes       |   no       |   no            |
| MPI_Send                 |   yes       |   yes       |   yes            |
| MPI_Ssend                 |   yes       |   yes       |   no            |
| MPI_Bsend                 |   no       |   no       |   no            |
| MPI_Rsend                 |   no       |   no       |   no            |
| MPI_Isend                |   yes       |   yes       |   no            |
| MPI_Issend                 |   no       |   no       |   no            |
| MPI_Ibsend                 |   no       |   no       |   no            |
| MPI_Irsend                 |   no       |   no       |   no            |
| MPI_Recv                 |   yes       |   yes       |   yes            |
| MPI_Irecv                |   yes       |   yes       |   no            |
| MPI_Sendrecv             |   yes       |   yes       |   yes            |
| MPI_Sendrecv_replace     |   yes       |   yes       |   yes            |

`critter` is currently not able to track user-defined symbols in any nonblocking collectives, including p2p nonblocking communications requiring `MPI_Waitany` or `MPI_Waitsome`. In addition, `critter` does not track sendrecv blocking communications in which a process sends and receives from distinct processes.
