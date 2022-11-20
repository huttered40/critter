#ifndef CRITTER__INTERFACE_H_
#define CRITTER__INTERFACE_H_

#include <mpi.h>

void _init(int* argc, char*** argv);

// General purpose routines
/* start: establishes the start of a window within which critter operates */
void critter_start(MPI_Comm cm=MPI_COMM_WORLD);
/* stop: establishes the end of a window within which critter operates */
void critter_stop(MPI_Comm cm=MPI_COMM_WORLD);
/* record: prints statistics dependent on mechanism */
void critter_record(int variantID=-1);

// Profiling interface routines
/* init: establishes a fixed order of kernels to keep track of */
void critter_register_timer(const char* timer_name);
/* start_time: starts timing a particular kernel */
void critter_start_timer(const char* timer_name, bool propagate_within = true, MPI_Comm cm = MPI_COMM_NULL);
/* stop_time: stops timing a particular kernel */
void critter_stop_timer(const char* timer_name, MPI_Comm cm = MPI_COMM_NULL);
/* get_critical_path_costs: returns number of critical-path costs (dependent on mechanism) */
int critter_get_critical_path_costs();
/* get_critical_path_costs: returns critical-path costs (dependent on mechanism) */
void critter_get_critical_path_costs(float* costs);
/* get_critical_path_costs: returns number of max-per-process costs (dependent on mechanism) */
int critter_get_max_per_process_costs();
/* get_critical_path_costs: returns max-per-process costs (dependent on mechanism) */
void critter_get_max_per_process_costs(float* costs);
/* get_critical_path_costs: returns number of volumetric costs (dependent on mechanism) */
int critter_get_volumetric_costs();
/* get_critical_path_costs: returns volumetric costs (dependent on mechanism) */
void critter_get_volumetric_costs(float* costs);

// Auxiliary routines
/* set_mode: resets the mechanism (possibly different than environment variable CRITTER_MECHANISM */
void critter_set_mode(int input_mode=-1);

#endif /*CRITTER__INTERFACE_H_*/
