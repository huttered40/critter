#ifndef CRITTER__DISCRETIZATION__UTIL__UTIL_H_
#define CRITTER__DISCRETIZATION__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

double get_arithmetic_mean(const pattern_key_id& index);
double get_variance(const pattern_key_id& index);
double get_std_dev(const pattern_key_id& index);
double get_std_error(const pattern_key_id& index);
double get_confidence_interval(const pattern_key_id& index, double level = .95);
/*
double get_skewness(pattern_key_id index);
double get_kurtosis(pattern_key_id index);
double get_jacque_barra(pattern_key_id index);
*/
void error_test(const pattern_key_id& index);
int should_schedule(const pattern_key_id& index);
int should_schedule_global(const pattern_key_id& index);
void set_schedule(const pattern_key_id& index, bool schedule_decision);
void update(const pattern_key_id& index, volatile double comp_time, double flop_count);

void allocate(MPI_Comm comm);
void open_symbol(const char* symbol, double curtime);
void close_symbol(const char* symbol, double curtime);
void final_accumulate(double last_time);
void reset();
void clear();

}
}
}

#endif /*CRITTER__DISCRETIZATION__UTIL__UTIL_H_*/
