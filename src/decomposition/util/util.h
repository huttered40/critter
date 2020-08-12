#ifndef CRITTER__DECOMPOSITION__UTIL__UTIL_H_
#define CRITTER__DECOMPOSITION__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace decomposition{

void allocate(MPI_Comm comm);
void reset(bool track_statistical_data_override, bool clear_statistical_data, bool schedule_kernels_override, bool propagate_statistical_data_overide);
void open_symbol(const char* symbol, double curtime);
void close_symbol(const char* symbol, double curtime);
void final_accumulate(double last_time);
void clear();

}
}
}

#endif /*CRITTER__DECOMPOSITION__UTIL__UTIL_H_*/
