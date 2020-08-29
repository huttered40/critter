#ifndef CRITTER__DISCRETIZATION__UTIL__UTIL_H_
#define CRITTER__DISCRETIZATION__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

extern int analysis_mode;
extern int is_optimized;
extern int autotuning_propagate;
extern int schedule_kernels;
extern int update_analysis;
extern int comm_sample_include_idle;
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
extern comm_pattern_key previous_comm_key;
extern std::map<std::pair<comm_pattern_key,comm_pattern_key>,idle_pattern> comm_pattern_pair_map;

bool is_key_skipable(const comm_pattern_key& key);
bool is_key_skipable(const comp_pattern_key& key);

double get_estimate(const pattern& p, int analysis_param, double unit_count=1.);
double get_estimate(const pattern_key_id& index, int analysis_param, double unit_count=1.);

double get_arithmetic_mean(const pattern& p);
double get_arithmetic_mean(const pattern_key_id& index);

double get_harmonic_mean(const pattern& p);
double get_harmonic_mean(const pattern_key_id& index);

double get_variance(const pattern& p, int analysis_param);
double get_variance(const pattern_key_id& index, int analysis_param);

double get_std_dev(const pattern& p, int analysis_param);
double get_std_dev(const pattern_key_id& index, int analysis_param);

double get_std_error(const pattern& p, int analysis_param);
double get_std_error(const pattern_key_id& index, int analysis_param);

double get_confidence_interval(const pattern& p, int analysis_param, double level = .95);
double get_confidence_interval(const pattern_key_id& index, int analysis_param, double level = .95);

bool is_steady(const pattern& p, int analysis_param);
bool is_steady(const pattern_key_id& index, int analysis_param);

bool steady_test(const comm_pattern_key& key, const pattern& p, int analysis_param);
bool steady_test(const comm_pattern_key& key, const pattern_key_id& index, int analysis_param);
bool steady_test(const comp_pattern_key& key, const pattern& p, int analysis_param);
bool steady_test(const comp_pattern_key& key, const pattern_key_id& index, int analysis_param);

void update_kernel_stats(pattern& p, int analysis_param, volatile double exec_time, double unit_count);
void update_kernel_stats(const pattern_key_id& index, int analysis_param, volatile double exec_time, double unit_count);
void update_kernel_stats(pattern& dest, const pattern& src, int analysis_param);
void update_kernel_stats(idle_pattern& p, bool is_global_steady_state, volatile double exec_time);

int should_schedule(const pattern& p);
int should_schedule(const pattern_key_id& index);

int should_schedule_global(const pattern& p);
int should_schedule_global(const pattern_key_id& index);

void set_kernel_state(pattern& p, bool schedule_decision);
void set_kernel_state(const pattern_key_id& index, bool schedule_decision);

void set_kernel_state_global(pattern& p, bool schedule_decision);
void set_kernel_state_global(const pattern_key_id& index, bool schedule_decision);

void allocate(MPI_Comm comm);
void open_symbol(const char* symbol, double curtime);
void close_symbol(const char* symbol, double curtime);
void final_accumulate(MPI_Comm comm, double last_time);
void reset(bool track_statistical_data_override, bool schedule_kernels_override, bool force_steady_statistical_data_overide, bool update_statistical_data_overide);
void clear();

}
}
}

#endif /*CRITTER__DISCRETIZATION__UTIL__UTIL_H_*/
