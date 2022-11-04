#include "volumetric.h"
#include "../util/util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace profile{

void volumetric::collect(MPI_Comm cm){
  if (mode==0) return;
  // First compute per-process max
  int rank; MPI_Comm_rank(cm,&rank);
  float_int buffer[num_pp_measures];
  for (size_t i=0; i<num_pp_measures; i++){
    max_pp_costs[i] = vol_costs[i];
    buffer[i].first = vol_costs[i];
    buffer[i].second = rank;
  }
  PMPI_Allreduce(MPI_IN_PLACE, &max_pp_costs[0], num_pp_measures, MPI_FLOAT, MPI_MAX, cm);
  PMPI_Allreduce(MPI_IN_PLACE, &buffer[0], num_pp_measures, MPI_FLOAT_INT, MPI_MAXLOC, cm);
  if (path_decomposition == 2){
    // First 'num_pp_measures' will be overwritten by the same data
    std::memcpy(&vol_costs[num_vol_measures], &max_pp_costs[num_pp_measures], (max_pp_costs.size()-num_pp_measures)*sizeof(float));
  }
  if (path_decomposition == 1 && path_count>0){
    size_t save=0;
    for (size_t i=0; i<path_index.size(); i++){
      size_t z = path_index[i];
      // If I am the processor that incurs the largest per-process time along any one path, then
      //   save the data and communicate it via reduction.
      if (rank == buffer[z].second){
        for (size_t j=0; j<num_decomp_pp_measures*list_size; j++){
          // The magic number 2 is for tracking the CompCost and ExecutionTime along each path
          max_pp_costs[num_pp_measures+save*(num_decomp_pp_measures*list_size+2)+j] = vol_costs[num_vol_measures+j];
        }
        // Set the two final measurements
        max_pp_costs[num_pp_measures+(save+1)*(num_decomp_pp_measures*list_size+2)-2] = vol_costs[num_vol_measures-3];// CompCost
        max_pp_costs[num_pp_measures+(save+1)*(num_decomp_pp_measures*list_size+2)-1] = vol_costs[num_vol_measures-1];// ExecTime
      }
      else{
        for (size_t j=0; j<num_decomp_pp_measures*list_size+2; j++){
          max_pp_costs[num_pp_measures+save*(num_decomp_pp_measures*list_size+2)+j] = 0.;
        }
      }
      PMPI_Allreduce(MPI_IN_PLACE, &max_pp_costs[num_pp_measures+save*(num_decomp_pp_measures*list_size+2)], num_decomp_pp_measures*list_size+2, MPI_FLOAT, MPI_MAX, cm);
      save++;
    }
  }
  else if (path_decomposition == 2 && path_count>0){
    size_t path_select_offset = num_kernel_ds*num_decomp_pp_measures+1;
    for (size_t i=0; i<path_index.size(); i++){
      size_t z = path_index[i];
      if (rank == buffer[z].second){
        for (size_t j=0; j<symbol_timers.size(); j++){
          std::memcpy(&max_pp_costs[num_pp_measures+j*path_select_offset*path_count+i*path_select_offset],
                      &vol_costs[num_pp_measures+j*path_select_offset*path_count+i*path_select_offset],
                      path_select_offset*sizeof(float));
        }
      }
      else{
        for (size_t j=0; j<symbol_timers.size(); j++){
          std::memset(&max_pp_costs[num_pp_measures+j*path_select_offset*path_count+i*path_select_offset],
                      (float)0.,
                      path_select_offset*sizeof(float));
        }
      }
    }
    PMPI_Allreduce(MPI_IN_PLACE, &max_pp_costs[0], (max_pp_costs.size()-num_pp_measures), MPI_FLOAT, MPI_MAX, cm);
  }
  // Now compute volumetric average
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  PMPI_Allreduce(MPI_IN_PLACE, &vol_costs[0], vol_costs.size(), MPI_FLOAT, MPI_SUM, cm);
  for (int i=0; i<vol_costs.size(); i++){ vol_costs[i] /= world_size; }
/*
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  if (mode && symbol_path_select_size>0){
    size_t active_size = world_size;
    size_t active_rank = world_rank;
    size_t active_mult = 1;
    while (active_size>1){
      if (active_rank % 2 == 1){
        int partner = (active_rank-1)*active_mult;
        PMPI_Send(...)
        break;
      }
      else if ((active_rank % 2 == 0) && (active_rank < (active_size-1))){
        int partner = (active_rank+1)*active_mult;
        PMPI_Recv(...)
      }
      active_size = active_size/2 + active_size%2;
      active_rank /= 2;
      active_mult *= 2;
    }
  }
*/
}

}
}
}
