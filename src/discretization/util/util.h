#ifndef CRITTER__DISCRETIZATION__UTIL__UTIL_H_
#define CRITTER__DISCRETIZATION__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

double get_arithmetic_mean(std::pair<bool,int> index);
double get_variance(std::pair<bool,int> index);
double get_std_dev(std::pair<bool,int> index);
double get_std_error(std::pair<bool,int> index);
double get_confidence_interval(std::pair<bool,int> index, double level = .95);
/*
double get_skewness(std::pair<bool,int> index);
double get_kurtosis(std::pair<bool,int> index);
double get_jacque_barra(std::pair<bool,int> index);
*/
void error_test(std::pair<bool,int> index);
bool should_schedule(std::pair<bool,int> index);
void update(std::pair<bool,int> index, volatile double comp_time, double flop_count);

}
}
}

#endif /*CRITTER__DISCRETIZATION__UTIL__UTIL_H_*/
