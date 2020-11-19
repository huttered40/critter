#include "comp.h"
#include "../util/util.h"
#include "../dispatch/dispatch.h"

namespace critter{
namespace internal{

void _blk_to_cyc_rect_(double* blocked, double* cyclic, int num_rows_local, int num_columns_local, int sliceDim){
  if (mode){
    volatile double curtime = MPI_Wtime();
    double flops = 0;
    bool schedule_decision = initiate_comp(_CAPITAL_blktocyc__id,curtime,flops,num_rows_local,num_columns_local,sliceDim);
    if (schedule_decision){
      int write_idx = 0; int read_idx = 0;
      int offset = num_rows_local*num_columns_local;
      int num_rows_global = num_rows_local*sliceDim;
      int num_columns_global = num_columns_local*sliceDim;
      for (int i=0; i<num_columns_local; i++){
        for (int j=0; j<sliceDim; j++){
          for (int k=0; k<num_rows_local; k++){
            for (int z=0; z<sliceDim; z++){
              read_idx = z*offset*sliceDim + k + j*offset + i*num_rows_local;
              cyclic[write_idx++] = blocked[read_idx];
            }
          }
        }
      }
      // Remove scalars that should be zeros
      for (int i=0; i<num_columns_global; i++){
        for (int j=i+1; j<num_rows_global; j++){
          cyclic[i*num_rows_global+j]=0.;
        }
      }
    }
    complete_comp(0,_CAPITAL_blktocyc__id,flops,num_rows_local,num_columns_local,sliceDim);
  } else{
    int write_idx = 0; int read_idx = 0;
    int offset = num_rows_local*num_columns_local;
    int num_rows_global = num_rows_local*sliceDim;
    int num_columns_global = num_columns_local*sliceDim;
    for (int i=0; i<num_columns_local; i++){
      for (int j=0; j<sliceDim; j++){
        for (int k=0; k<num_rows_local; k++){
          for (int z=0; z<sliceDim; z++){
            read_idx = z*offset*sliceDim + k + j*offset + i*num_rows_local;
            cyclic[write_idx++] = blocked[read_idx];
          }
        }
      }
    }
    // Remove scalars that should be zeros
    for (int i=0; i<num_columns_global; i++){
      for (int j=i+1; j<num_rows_global; j++){
        cyclic[i*num_rows_global+j]=0.;
      }
    }
  }
}
void _cyc_to_blk_rect_(double* blocked, double* cyclic, int num_rows_local, int num_columns_local, int sliceDim){
  if (mode){
    volatile double curtime = MPI_Wtime();
    double flops = 0;
    bool schedule_decision = initiate_comp(_CAPITAL_cyctoblk__id,curtime,flops,num_rows_local,num_columns_local,sliceDim);
    if (schedule_decision){
      int write_idx = 0; int read_idx = 0; int offset = num_rows_local*num_columns_local;
      for (int i=0; i<num_columns_local; i++){
        for (int j=0; j<sliceDim; j++){
          for (int k=0; k<num_rows_local; k++){
            for (int z=0; z<sliceDim; z++){
              write_idx = z*offset*sliceDim + k + j*offset + i*num_rows_local;
              blocked[write_idx] = cyclic[read_idx++];
            }
          }
        }
      }
    }
    complete_comp(0,_CAPITAL_cyctoblk__id,flops,num_rows_local,num_columns_local,sliceDim);
  } else{
    int write_idx = 0; int read_idx = 0; int offset = num_rows_local*num_columns_local;
    for (int i=0; i<num_columns_local; i++){
      for (int j=0; j<sliceDim; j++){
        for (int k=0; k<num_rows_local; k++){
          for (int z=0; z<sliceDim; z++){
            write_idx = z*offset*sliceDim + k + j*offset + i*num_rows_local;
            blocked[write_idx] = cyclic[read_idx++];
          }
        }
      }
    }
  }
}

}
}
