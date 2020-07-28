#include "volumetric.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace decomposition{

void volumetric::collect(MPI_Comm cm){
  // First compute per-process max
  int rank; MPI_Comm_rank(cm,&rank);
  double_int buffer[num_per_process_measures];
  for (size_t i=0; i<num_per_process_measures; i++){
    max_per_process_costs[i] = volume_costs[i];
    buffer[i].first          = volume_costs[i];
    buffer[i].second         = rank;
  }
  PMPI_Allreduce(MPI_IN_PLACE, &max_per_process_costs[0], num_per_process_measures, MPI_DOUBLE, MPI_MAX, cm);
  PMPI_Allreduce(MPI_IN_PLACE, &buffer[0], num_per_process_measures, MPI_DOUBLE_INT, MPI_MAXLOC, cm);
  size_t save=0;
  for (size_t i=0; i<comm_path_select.size(); i++){// don't consider idle time an option
    if (comm_path_select[i] == '0') continue;
    size_t z = i<(2*cost_model_size) ? i : i+1;	// careful indexing to avoid idle time
    if (rank == buffer[z].second){
      for (size_t j=0; j<num_tracker_per_process_measures*list_size; j++){
        max_per_process_costs[num_per_process_measures+save*(num_tracker_per_process_measures*list_size+2)+j] = volume_costs[num_volume_measures+j];
      }
      max_per_process_costs[num_per_process_measures+(save+1)*(num_tracker_per_process_measures*list_size+2)-3] = volume_costs[num_volume_measures-6];
      max_per_process_costs[num_per_process_measures+(save+1)*(num_tracker_per_process_measures*list_size+2)-2] = volume_costs[num_volume_measures-2];
      max_per_process_costs[num_per_process_measures+(save+1)*(num_tracker_per_process_measures*list_size+2)-1] = volume_costs[num_volume_measures-5];
    }
    else{
      for (size_t j=0; j<num_tracker_per_process_measures*list_size+2; j++){
        max_per_process_costs[num_per_process_measures+save*(num_tracker_per_process_measures*list_size+2)+j] = 0.;
      }
    }
    PMPI_Allreduce(MPI_IN_PLACE, &max_per_process_costs[num_per_process_measures+save*(num_tracker_per_process_measures*list_size+2)], num_tracker_per_process_measures*list_size+2, MPI_DOUBLE, MPI_MAX, cm);
    save++;
  }
  // For now, buffer[num_per_process_measures-1].second holds the rank of the process with the max per-process runtime
  if (mode && symbol_path_select_size>0){
    // copy data to volume buffers to avoid corruption
    std::memcpy(&symbol_timer_pad_local_vol[0], &symbol_timer_pad_local_pp[0], (vol_symbol_class_count*num_volume_measures+1)*max_num_symbols*sizeof(double));

    int per_process_runtime_root_rank = buffer[num_per_process_measures-1].second;
    // We consider only critical path runtime
    int ftimer_size = 0;
    if (rank==per_process_runtime_root_rank){
      ftimer_size = symbol_timers.size();
    }
    PMPI_Allreduce(MPI_IN_PLACE,&ftimer_size,1,MPI_INT,MPI_SUM,cm);

    for (auto i=0; i<symbol_len_pad_cp.size(); i++){ symbol_len_pad_cp[i]=0.; }
    if (rank==per_process_runtime_root_rank){
      int symbol_offset = 0;
      for (auto i=0; i<symbol_timers.size(); i++){
        symbol_len_pad_cp[i] = symbol_order[i].size();
        for (auto j=0; j<symbol_len_pad_cp[i]; j++){
          symbol_pad_cp[symbol_offset+j] = symbol_order[i][j];
        }
        symbol_offset += symbol_len_pad_cp[i];
      }
    }
    PMPI_Allreduce(MPI_IN_PLACE,&symbol_len_pad_cp[0],ftimer_size,MPI_INT,MPI_SUM,cm);

    int num_chars = 0;
    for (auto i=0; i<ftimer_size; i++){
      num_chars += symbol_len_pad_cp[i];
    }
    if (rank == per_process_runtime_root_rank){
      PMPI_Bcast(&symbol_timer_pad_local_pp[0],(pp_symbol_class_count*num_per_process_measures+1)*ftimer_size,MPI_DOUBLE,rank,cm);
      PMPI_Bcast(&symbol_pad_cp[0],num_chars,MPI_CHAR,rank,cm);
    }
    else{
      PMPI_Bcast(&symbol_timer_pad_global_pp[0],(pp_symbol_class_count*num_per_process_measures+1)*ftimer_size,MPI_DOUBLE,per_process_runtime_root_rank,cm);
      PMPI_Bcast(&symbol_pad_cp[0],num_chars,MPI_CHAR,per_process_runtime_root_rank,cm);
      int symbol_offset = 0;
      for (int i=0; i<ftimer_size; i++){
        auto reconstructed_symbol = std::string(symbol_pad_cp.begin()+symbol_offset,symbol_pad_cp.begin()+symbol_offset+symbol_len_pad_cp[i]);
        if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
          symbol_timers[reconstructed_symbol] = symbol_tracker(reconstructed_symbol);
          symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
        }
        *symbol_timers[reconstructed_symbol].pp_numcalls = symbol_timer_pad_global_pp[(pp_symbol_class_count*num_per_process_measures+1)*i];
        for (int j=0; j<num_per_process_measures; j++){
          symbol_timers[reconstructed_symbol].pp_incl_measure[j] = symbol_timer_pad_global_pp[(pp_symbol_class_count*num_per_process_measures+1)*i+j+1];
          symbol_timers[reconstructed_symbol].pp_excl_measure[j] = symbol_timer_pad_global_pp[(pp_symbol_class_count*num_per_process_measures+1)*i+num_per_process_measures+j+1];
        }
        symbol_timers[reconstructed_symbol].has_been_processed = true;
        symbol_offset += symbol_len_pad_cp[i];
      }
      // Now cycle through and find the symbols that were not processed and set their accumulated measures to 0
      for (auto& it : symbol_timers){
        if (it.second.has_been_processed){ it.second.has_been_processed = false; }
        else{
          *it.second.pp_numcalls = 0;
          for (int j=0; j<num_per_process_measures; j++){
            it.second.pp_incl_measure[j] = 0;
            it.second.pp_excl_measure[j] = 0;
          }
        }
      }
    }
  }

  // Now compute volumetric average
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  if (mode){
    PMPI_Allreduce(MPI_IN_PLACE, &volume_costs[0], volume_costs.size(), MPI_DOUBLE, MPI_SUM, cm);
    for (int i=0; i<volume_costs.size(); i++){ volume_costs[i] /= (1.*world_size); }
  }
  if (mode && symbol_path_select_size>0){
    size_t active_size = world_size;
    size_t active_rank = world_rank;
    size_t active_mult = 1;
    while (active_size>1){
      if (active_rank % 2 == 1){
        int partner = (active_rank-1)*active_mult;
        int ftimer_size = symbol_timers.size();
        PMPI_Send(&ftimer_size, 1, MPI_INT, partner, internal_tag, cm);
        for (auto i=0; i<symbol_len_pad_cp.size(); i++){ symbol_len_pad_cp[i]=0.; }
        int symbol_offset = 0;
        for (auto i=0; i<symbol_timers.size(); i++){
          symbol_len_pad_cp[i] = symbol_order[i].size();
          for (auto j=0; j<symbol_len_pad_cp[i]; j++){
            symbol_pad_cp[symbol_offset+j] = symbol_order[i][j];
          }
          symbol_offset += symbol_len_pad_cp[i];
        }
        PMPI_Send(&symbol_len_pad_cp[0],ftimer_size,MPI_INT,partner,internal_tag1,cm);
        int num_chars = 0;
        for (auto i=0; i<ftimer_size; i++){
          num_chars += symbol_len_pad_cp[i];
        }
        PMPI_Send(&symbol_timer_pad_local_vol[0],(vol_symbol_class_count*num_volume_measures+1)*ftimer_size,MPI_DOUBLE,partner,internal_tag2,cm);
        PMPI_Send(&symbol_pad_cp[0],num_chars,MPI_CHAR,partner,internal_tag3,cm);
        break;
      }
      else if ((active_rank % 2 == 0) && (active_rank < (active_size-1))){
        int partner = (active_rank+1)*active_mult;
        int ftimer_size_foreign;
        PMPI_Recv(&ftimer_size_foreign, 1, MPI_INT, partner, internal_tag, cm, MPI_STATUS_IGNORE);
        PMPI_Recv(&symbol_len_pad_cp[0],ftimer_size_foreign,MPI_INT,partner, internal_tag1, cm, MPI_STATUS_IGNORE);
        int num_chars = 0;
        for (auto i=0; i<ftimer_size_foreign; i++){
          num_chars += symbol_len_pad_cp[i];
        }
        PMPI_Recv(&symbol_timer_pad_global_vol[0],(vol_symbol_class_count*num_volume_measures+1)*ftimer_size_foreign,MPI_DOUBLE,partner,internal_tag2,cm, MPI_STATUS_IGNORE);
        PMPI_Recv(&symbol_pad_cp[0],num_chars,MPI_CHAR,partner, internal_tag3,cm, MPI_STATUS_IGNORE);
        int symbol_offset = 0;
        for (int i=0; i<ftimer_size_foreign; i++){
          auto reconstructed_symbol = std::string(symbol_pad_cp.begin()+symbol_offset,symbol_pad_cp.begin()+symbol_offset+symbol_len_pad_cp[i]);
          if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
            symbol_timers[reconstructed_symbol] = symbol_tracker(reconstructed_symbol);
            symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
          }
          *symbol_timers[reconstructed_symbol].vol_numcalls += symbol_timer_pad_global_vol[(vol_symbol_class_count*num_per_process_measures+1)*i];
          for (int j=0; j<num_volume_measures; j++){
            symbol_timers[reconstructed_symbol].vol_incl_measure[j] += symbol_timer_pad_global_vol[(vol_symbol_class_count*num_volume_measures+1)*i+j+1];
            symbol_timers[reconstructed_symbol].vol_excl_measure[j] += symbol_timer_pad_global_vol[(vol_symbol_class_count*num_volume_measures+1)*i+num_volume_measures+j+1];
          }
          symbol_timers[reconstructed_symbol].has_been_processed = true;
          symbol_offset += symbol_len_pad_cp[i];
        }
      }
      active_size = active_size/2 + active_size%2;
      active_rank /= 2;
      active_mult *= 2;
    }
  }
}

}
}
}
