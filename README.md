
# critter
Welcome! If you are looking for a lightweight tool to analyze the critical path costs of your distributed-memory MPI program, you have come to the right place. `critter` **seeks to understand the critical paths of your MPI program**, and tracks the critical path of each of the following measurements:

1. total runtime
2. communication time
3. computation time
4. number of bytes communicated
5. estimated communication cost
6. estimated synchronization cost

`critter` also provides both **per-process** and **volumetric** costs of the measures above.

In addition, `critter` can break down each of the critical path measures into the contributions from each MPI routine.

Critical paths through parallel schedules incur contributions from many functions; thus we are currently engineering a way to track the contributions (in terms of the measures listed above) of each function.

See the lists below for an accurate depiction of our current support.

## Build and use instructions
`configure` compiler and flags in `config/config.mk` (MPI installation and C++11 are required). Run `make` in the main directory to generate the library file `./lib/libcritter.a`. Include `critter.h` in all files that use MPI in your application (i.e. replace `include mpi.h`), and link to `./lib/libcritter.a`.

`critter` provides two routines to the user: `critter::start(...)` and `critter::stop(...)`. These create the window within which all MPI routines are intercepted and tracked.

`critter` provies a few variables to the user inside the `critter` namespace in `src/critter.h`. Along with a tag `critter::internal_tag` to prevent critter's internal MPI communication from conflicting with user communication, `critter_breakdown` specifies which, if any, critical path measurement is broken down into contributions from individual MPI routines. Note that this extra information comes at a cost of increased internal data transfer necessary to propogate critical path information. `critter::critter_breakdown_size` must match `critter::critter_breakdown.count()`.

## Design decisions
Critter will never assume a communication protocol more limiting than what is specified in the MPI Standard to propogate critical paths. As such, it will never break the semantics of your parallel program and limit its forward progress. If and only if the MPI implementation performs a communication protocol more limiting than is necessary (i.e. requiring synchronous handshake between processes calling MPI_Send and MPI_Receive), it may report an erroneous critical path measure.
1. Critter assumes any blocking collective communication or synchronous p2p communication performs a synchronization at the start, and thereby does not force a limiting communication protocol in using a blocking collective or synchronous p2p routine to propogate the critical paths before the intercepted routine takes place. It allows processes to jump back into the user code without further synchronization once the intercepted collective routine completes, and therefore does not synchronize for the limiting costs of the intercepted routine.
2. Critter assumes that in any asynchronous p2p communication (including nonblocking and blocking protocols), the sender does not wait on the receiver, and thus needs not accumulate the receiver's critical path information. Unlike critter's support for synchronous protocol, if the intercepted communication is blocking, the receiver will update its critical path after the intercepted routine completes, thus taking into account the limiting costs of the intercepted routine into the propogated critical paths. If the MPI implementation requires a handshake between sender and receiver in a blocking p2p communication routine, critter will accurately 
3. All nonblocking communication routines, including p2p and collectives, use nonblocking communication to propogate critical paths. The sending process sends a nonblocking message to the receiver after posting its intercepted routine, to be completed directly after the completion of the intercepted routine. The receiver will also post a nonblocking message to retrieve the sender's critical path information, and will compare that against its critical path data after the intercepted routine completes. Such critical path propogation will not use stale sending-process data, as the receiver is not dependent on the sender until exactly when the intercepted routine completion takes place.

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
| MPI_Ssend                 |   yes       |   no       |   no            |
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
