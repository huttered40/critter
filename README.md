
# critter
Welcome! If you are looking for a lightweight tool to analyze the critical path costs of your distributed-memory MPI program, you have come to the right place. `critter` seeks to understand the critical paths of your MPI program, and counts the following metrics for each:

1. number of bytes communicated
2. communication time
3. idle time
4. estimated communication cost
5. estimated synchronization cost
6. computation time
7. total runtime

`critter` also provides volumetric costs for each of the measures above, and tracks these measures for each MPI routine.

`critter` provides three levels of support:
1. counting critical/per-process path measures and outputting data to the standard output file
2. both counting and tracking the critical/per-process path measures and critical path itself, as identified by local synchronization points inherent in your program (any intercepted MPI routine)
3. micro-benchmarks for predicting based on critical path from level2

Using both the micro-benchmarks and the tracked critical path, `critter` can re-assemble the communication costs and provide an accurate assessment for strong scaling and weak scaling tests for any node count.

See the lists below for an accurate depiction of our current support.

## Build and use instructions
`configure` compiler and flags in `config/config.mk` (MPI installation and C++11 are required). Run `make` in the main directory to generate the library file `./lib/libcritter.a`. Include `critter.h` in all files that use MPI in your application, remove all `include mpi.h`, and link to `./lib/libcritter.a`. Make sure to run `export CRITTER_VIZ=ON` before building your application for levels 2 and 3 support.

`critter` provides three routines to the user: `critter::start()`, `critter::stop()`, and `critter::print(int size, double* data)`. The first two create the window within which all MPI routines are intercepted and tracked. The latter is optional and is provided to allow the user to include diagnostic data for the test, including information such as residual or raw performance.

Each of `critter`'s three levels of support require increasing user input. If you want level 1 support, simply follow the directions above, create your own scripts, and look for `critter`'s output to stdout.

Levels 2 and 3 require more user input, but provide much stronger tools for analysis and prediction, and save all of the hassle in writing script files and plotting data.

Depending on the machine, the user should run `python setup.py develop --user` to register the python library. We have provided a few sample instruction files under `experiments/` that should be modified for your needs. Each invocation of the instructions file will generate a test directory under both `Tests/` and `$SCRATCH/`, the latter of which should contain all relevant data.

The user should build his/her program separately, and then inject the corresponding binary into `critter`'s `bin/` directory. All binary names should correspond to the `algorithm` class instance tag specified in the instructions file.

In addition to writing the instructions file, the user need also write a file under `experiments/machines` to specify exactly how the machine handles a few key details (such as providing the launch and scheduling information). See one of the provided example files.

## Advice
Large-scale tests can be expensive, and if you are using `critter`'s levels 2/3 support, we want to make sure you agree with the generated tests.
1. always run `generate` without `launch` first and inspect the generated script files in the test directory. Once satisfied, erase the generated test directories (optional but convenient) and add the `launch` method directly after the `generate` method
2. remove all binaries in the `bin/` (in test directory)  before applying `tar -cvf`
3. plot names can be specified in `instructions.py` for each test
4. due to scaplot assumptions, users can only specify node counts that range in a multiplicative manner
5. don't have any `'+'` characters in the file path to `critter`, nor in any of the arguments. `critter` uses this character to parse file strings
6. make sure your binaries are not named `test_allreduce` or similar for other benchmarks
7. make sure any communication requiring a tag doesn't use a tag of the same value as critter' internal tag critter::internal::tag. If there is a conflict, simply modify critter's internal tag and rebuild instead of modifying your source code.

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
| MPI_Ibcast               |   yes       |   yno      |   no           |
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
