#ifndef CRITTER__INTERCEPT__CUSTOM_H_
#define CRITTER__INTERCEPT__CUSTOM_H_

namespace critter{
namespace internal{

void _blk_to_cyc_rect_(double* blocked, double* cyclic, int num_rows_local, int num_columns_local, int sliceDim);
void _cyc_to_blk_rect_(double* blocked, double* cyclic, int num_rows_local, int num_columns_local, int sliceDim);
}
}

#endif /*CRITTER__INTERCEPT__COMP_H_*/
