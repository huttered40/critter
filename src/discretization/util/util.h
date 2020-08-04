#ifndef CRITTER__DISCRETIZATION__UTIL__UTIL_H_
#define CRITTER__DISCRETIZATION__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

double get_arithmetic_mean(pattern_key_id index);
double get_variance(pattern_key_id index);
double get_std_dev(pattern_key_id index);
double get_std_error(pattern_key_id index);
double get_confidence_interval(pattern_key_id index, double level = .95);
/*
double get_skewness(pattern_key_id index);
double get_kurtosis(pattern_key_id index);
double get_jacque_barra(pattern_key_id index);
*/
void error_test(pattern_key_id index);
bool should_schedule(pattern_key_id index);
void update(pattern_key_id index, volatile double comp_time, double flop_count);

}
}
}

#endif /*CRITTER__DISCRETIZATION__UTIL__UTIL_H_*/
