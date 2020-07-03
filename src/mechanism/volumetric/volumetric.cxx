#include "volumetric.h"
#include "../../container/comm_tracker.h"
#include "../../container/symbol_tracker.h"

namespace critter{
namespace internal{

void volumetric::collect(MPI_Comm cm){
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  if (mode>=1){
    PMPI_Allreduce(MPI_IN_PLACE, &volume_costs[0], volume_costs.size(), MPI_DOUBLE, MPI_SUM, cm);
    for (int i=0; i<volume_costs.size(); i++){ volume_costs[i] /= (1.*world_size); }
  }
  if (mode>=2){
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
        PMPI_Send(&symbol_timer_pad_local_vol[0],(symbol_class_count*num_volume_measures+1)*ftimer_size,MPI_DOUBLE,partner,internal_tag2,cm);
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
        PMPI_Recv(&symbol_timer_pad_global_vol[0],(symbol_class_count*num_volume_measures+1)*ftimer_size_foreign,MPI_DOUBLE,partner,internal_tag2,cm, MPI_STATUS_IGNORE);
        PMPI_Recv(&symbol_pad_cp[0],num_chars,MPI_CHAR,partner, internal_tag3,cm, MPI_STATUS_IGNORE);
        int symbol_offset = 0;
        for (int i=0; i<ftimer_size_foreign; i++){
          auto reconstructed_symbol = std::string(symbol_pad_cp.begin()+symbol_offset,symbol_pad_cp.begin()+symbol_offset+symbol_len_pad_cp[i]);
          if (symbol_timers.find(reconstructed_symbol) == symbol_timers.end()){
            symbol_timers[reconstructed_symbol] = symbol_tracker(reconstructed_symbol);
            symbol_order[(symbol_timers.size()-1)] = reconstructed_symbol;
          }
          *symbol_timers[reconstructed_symbol].vol_numcalls += symbol_timer_pad_global_vol[(symbol_class_count*num_per_process_measures+1)*i];
          for (int j=0; j<num_volume_measures; j++){
            *symbol_timers[reconstructed_symbol].vol_incl_measure[j] += symbol_timer_pad_global_vol[(symbol_class_count*num_volume_measures+1)*i+2*j+1];
            *symbol_timers[reconstructed_symbol].vol_excl_measure[j] += symbol_timer_pad_global_vol[(symbol_class_count*num_volume_measures+1)*i+2*(j+1)];
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
