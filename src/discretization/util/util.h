#ifndef CRITTER__DISCRETIZATION__UTIL__UTIL_H_
#define CRITTER__DISCRETIZATION__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

extern int autotuning_mode;
extern int autotuning_propagate;
extern int schedule_kernels;
extern MPI_Datatype comm_pattern_key_type;
extern MPI_Datatype comp_pattern_key_type;
extern MPI_Datatype pattern_type;
extern size_t pattern_count_limit;
extern double pattern_time_limit;
extern double pattern_error_limit;
extern std::map<MPI_Comm,std::pair<int,int>> communicator_map;
extern std::map<comm_pattern_key,pattern_key_id> comm_pattern_map;
extern std::map<comp_pattern_key,pattern_key_id> comp_pattern_map;
extern std::vector<comm_pattern_key> steady_state_comm_pattern_keys;
extern std::vector<comm_pattern_key> active_comm_pattern_keys;
extern std::vector<comp_pattern_key> steady_state_comp_pattern_keys;
extern std::vector<comp_pattern_key> active_comp_pattern_keys;
extern std::vector<pattern> steady_state_patterns;
extern std::vector<pattern> active_patterns;

double get_estimate(const pattern_key_id& index, int pattern_param, double unit_count=1.);
double get_arithmetic_mean(const pattern_key_id& index);
double get_harmonic_mean(const pattern_key_id& index);
double get_variance(const pattern_key_id& index, int pattern_param);
double get_std_dev(const pattern_key_id& index, int pattern_param);
double get_std_error(const pattern_key_id& index, int pattern_param);
double get_confidence_interval(const pattern_key_id& index, int pattern_param, double level = .95);
/*
double get_skewness(pattern_key_id index);
double get_kurtosis(pattern_key_id index);
double get_jacque_barra(pattern_key_id index);
*/
void error_test(const pattern_key_id& index, int pattern_param);
void update(const pattern_key_id& index, int pattern_param, volatile double comp_time, double flop_count);

int should_schedule(const pattern_key_id& index);
int should_schedule_global(const pattern_key_id& index);
void set_schedule(const pattern_key_id& index, bool schedule_decision);

void allocate(MPI_Comm comm);
void open_symbol(const char* symbol, double curtime);
void close_symbol(const char* symbol, double curtime);
void final_accumulate(double last_time);
void reset(bool track_statistical_data_override, bool clear_statistical_data, bool schedule_kernels_override, bool propagate_statistical_data_overide);
void clear();

}
}
}

#endif /*CRITTER__DISCRETIZATION__UTIL__UTIL_H_*/
