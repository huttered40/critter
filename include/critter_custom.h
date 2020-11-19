#ifndef __CRITTER_CUSTOM_H__
#define __CRITTER_CUSTOM_H__

#include "../src/intercept/custom.h"

// Capital
#define blk_to_cyc_rect(blocked,cyclic,num_rows_local,num_columns_local,sliceDim)\
    critter::internal::_blk_to_cyc_rect_(blocked,cyclic,num_rows_local,num_columns_local,sliceDim)

#define cyc_to_blk_rect(blocked,cyclic,num_rows_local,num_columns_local,sliceDim)\
    critter::internal::_cyc_to_blk_rect_(blocked,cyclic,num_rows_local,num_columns_local,sliceDim)

#endif /*CRITTER_CUSTOM_H_*/
