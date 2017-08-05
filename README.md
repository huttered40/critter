# critter
Counter and timer for critical path costs of MPI programs, uses PMPI.

## Build and use instructions
Configure compiler and flags in config.mk, MPI installation and C++11 are required.

Include critter.h in all files that use MPI in your application and link to ./lib/libcritter.a, when you run, timers and costs will be displayed when MPI_Finalize() is called.

The current list of MPI functions being tracked: MPI_Bcast
