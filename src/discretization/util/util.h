#ifndef CRITTER__DISCRETIZATION__UTIL__UTIL_H_
#define CRITTER__DISCRETIZATION__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

struct kernel_batch;
struct kernel_propagate;


struct sample_propagation_tree{
 solo_channel* root;
};

struct sample_propagation_forest{
  sample_propagation_forest();
  ~sample_propagation_forest();

  //void generate_span(channel* node, std::vector<std::pair<int,int>>& perm_tuples);
  void insert_node(solo_channel* node);
  void clear_info();
  void fill_ancestors(solo_channel* node, kernel_batch& batch);
  void fill_descendants(solo_channel* node, kernel_batch& batch);

  sample_propagation_tree* tree;
private:
  void delete_tree(solo_channel*& tree_root);
  void clear_tree_info(solo_channel* tree_root);
/*
  void generate_sibling_perm(std::vector<std::pair<int,int>>& static_info, std::vector<std::pair<int,int>>& gen_info, std::vector<std::pair<int,int>>& save_info, int level, bool& valid_siblings);
  void generate_partition_perm(std::vector<std::pair<int,int>>& static_info, std::vector<std::pair<int,int>>& gen_info, int level, bool& valid_partition,
                               int parent_max_span, int parent_min_stride);
*/
  //bool partition_test(channel* parent, int subtree_idx);
  void find_parent(solo_channel* tree_root, solo_channel* tree_node, solo_channel*& parent);
};

// ****************************************************************************************************************************************************
// Note: 'kernel_batch' is to be used for volumetric sampling modes only.
// Note: batches will never exist across variants, as they will be liquidated into the pathset during a variant, or else in 'final_accumulate'
struct kernel_batch{
  // If I leverage the kurtosis, I will have to utilize the arithmetic mean.
  //   Note that I'd rather utilize the geometric mean, but I'm not sure how to convert this algorithm
  //     to handle that.
  kernel_batch(channel* node = nullptr);
  kernel_batch(const kernel_batch& _copy);
  kernel_batch& operator=(const kernel_batch& _copy);

  int hash_id;
  int channel_count;
  int num_schedules;
  int num_local_schedules;
  double num_scheduled_units;
  double num_local_scheduled_units;
  double total_exec_time;
  double total_local_exec_time;
  double M1,M2;
  std::set<channel*> registered_channels;
};

// ****************************************************************************************************************************************************
// Note: 'kernel_batch_propagate' is to be used for volumetric sampling modes only.
struct kernel_batch_propagate{
  // If I leverage the kurtosis, I will have to utilize the arithmetic mean.
  //   Note that I'd rather utilize the geometric mean, but I'm not sure how to convert this algorithm
  //     to handle that.
  kernel_batch_propagate(){};
  kernel_batch_propagate(const kernel_batch& _copy);
  kernel_batch_propagate(const kernel_batch_propagate& _copy);
  kernel_batch_propagate& operator=(const kernel_batch_propagate& _copy);

  int hash_id;
  int channel_count;
  int num_schedules;
  double num_scheduled_units;
  double total_exec_time;
  double M1,M2;
};

// ****************************************************************************************************************************************************
struct kernel{
  // If I leverage the kurtosis, I will have to utilize the arithmetic mean.
  //   Note that I'd rather utilize the geometric mean, but I'm not sure how to convert this algorithm
  //     to handle that.
  kernel(channel* node = nullptr);
  kernel(const kernel_batch& _copy);
  kernel(const kernel& _copy);
  kernel& operator=(const kernel& _copy);
  kernel(const kernel_propagate& _copy);
  kernel& operator=(const kernel_propagate& _copy);
  void reset();
  void clear_distribution();

  int hash_id;
  int steady_state;
  int global_steady_state;
  int num_schedules;
  int num_non_schedules;
  int num_local_schedules;
  double num_scheduled_units;
  double num_non_scheduled_units;
  double num_local_scheduled_units;
  double M1,M2;
  double total_exec_time;
  double total_local_exec_time;
  std::set<channel*> registered_channels;
};

// ****************************************************************************************************************************************************
struct kernel_propagate{
  // If I leverage the kurtosis, I will have to utilize the arithmetic mean.
  //   Note that I'd rather utilize the geometric mean, but I'm not sure how to convert this algorithm
  //     to handle that.
  kernel_propagate(){};
  kernel_propagate(const kernel& _copy);
  kernel_propagate(const kernel_propagate& _copy);
  kernel_propagate& operator=(const kernel_propagate& _copy);

  int hash_id;
  int num_schedules;
  int num_local_schedules;
  double num_scheduled_units;
  double num_local_scheduled_units;
  double M1,M2;
  double total_exec_time;
  double total_local_exec_time;
};


// ****************************************************************************************************************************************************
struct intermediate_stats{
  intermediate_stats(const kernel_key_id& index, const std::vector<kernel_batch>& active_batches);
  void generate(const kernel& p, const std::vector<kernel_batch>& active_batches);

  int num_schedules;
  int num_local_schedules;
  double M1,M2;
  double num_scheduled_units;
  double num_local_scheduled_units;
  double total_exec_time;
  double total_local_exec_time;
};

// ****************************************************************************************************************************************************
extern std::ofstream stream;
extern int tuning_delta;
extern int comm_sample_aggregation_mode;
extern int comm_state_aggregation_mode;
extern int comp_sample_aggregation_mode;
extern int comp_state_aggregation_mode;
extern int sample_constraint_mode;
//extern int sample_reset_mode;
extern int schedule_kernels;
extern int update_analysis;
extern MPI_Datatype kernel_type;
extern MPI_Datatype batch_type;
extern size_t kernel_count_limit;
extern double kernel_time_limit;
extern double kernel_error_limit;
extern int comp_kernel_transfer_id;
extern int comm_kernel_transfer_id;
extern int comp_kernel_buffer_id;
extern std::map<comm_kernel_key,kernel_key_id> comm_kernel_map;
extern std::map<comp_kernel_key,kernel_key_id> comp_kernel_map;
extern std::vector<comm_kernel_key> active_comm_kernel_keys;
extern std::vector<comp_kernel_key> active_comp_kernel_keys;
extern std::vector<kernel> active_kernels;
extern sample_propagation_forest spf;
extern std::map<comm_kernel_key,std::vector<kernel_batch>> comm_batch_map;
extern std::map<comp_kernel_key,std::vector<kernel_batch>> comp_batch_map;
extern std::set<intptr_t> skip_ptr_set;

extern std::ofstream stream_tune,stream_reconstruct;
extern std::vector<double> intercept_overhead;
extern std::vector<double> global_intercept_overhead;
extern std::vector<double> global_comp_kernel_stats;
extern std::vector<double> global_comm_kernel_stats;
extern std::vector<double> local_comp_kernel_stats;
extern std::vector<double> local_comm_kernel_stats;
extern std::vector<double> save_comp_kernel_stats;
extern std::vector<double> save_comm_kernel_stats;
extern size_t num_critical_path_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime, CompTime, RunTime
extern size_t num_per_process_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, RunTime
extern size_t num_volume_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, RunTime
extern size_t num_tracker_critical_path_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t num_tracker_per_process_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t num_tracker_volume_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t critical_path_costs_size;
extern size_t per_process_costs_size;
extern size_t volume_costs_size;
extern std::vector<double> critical_path_costs;
extern std::vector<double> max_per_process_costs;
extern std::vector<double> volume_costs;
extern std::vector<double_int> info_sender;
extern std::vector<double_int> info_receiver;
extern std::vector<char> eager_pad;
extern double comp_start_time;
extern std::map<MPI_Request,std::pair<bool,int>> internal_comm_info;
extern std::map<MPI_Request,std::pair<MPI_Comm,int>> internal_comm_comm;
extern std::map<MPI_Request,std::pair<double,double>> internal_comm_data;
extern std::vector<std::pair<double*,int>> internal_comm_prop;
extern std::vector<MPI_Request> internal_comm_prop_req;
extern std::vector<int*> internal_timer_prop_int;
extern std::vector<double*> internal_timer_prop_double;
extern std::vector<double_int*> internal_timer_prop_double_int;
extern std::vector<char*> internal_timer_prop_char;
extern std::vector<MPI_Request> internal_timer_prop_req;
extern std::vector<bool> decisions;
extern std::map<std::string,std::vector<double>> save_info;
extern std::vector<double> new_cs;
extern size_t mode_1_width;
extern size_t mode_2_width;
extern int internal_tag;
extern int internal_tag1;
extern int internal_tag2;
extern int internal_tag3;
extern int internal_tag4;
extern int internal_tag5;
extern bool is_first_iter;

// ****************************************************************************************************************************************************
bool is_key_skipable(const comm_kernel_key& key);
bool is_key_skipable(const comp_kernel_key& key);

int get_skel_count(const comm_kernel_key& key);
int get_skel_count(const comp_kernel_key& key);

double get_estimate(const kernel& p, int analysis_param, double unit_count=1.);
double get_estimate(const kernel_propagate& p, int analysis_param, double unit_count=1.);
double get_estimate(const kernel_key_id& index, int analysis_param, double unit_count=1.);
double get_estimate(const kernel_key_id& index, const std::vector<kernel_batch>& active_batches, int analysis_param, double unit_count);
double get_estimate(const intermediate_stats& p, int analysis_param, double unit_count=1.);

double get_arithmetic_mean(const kernel& p);
double get_arithmetic_mean(const kernel_propagate& p);
double get_arithmetic_mean(const kernel_key_id& index);
double get_arithmetic_mean(const intermediate_stats& p);

double get_harmonic_mean(const kernel& p);
double get_harmonic_mean(const kernel_propagate& p);
double get_harmonic_mean(const kernel_key_id& index);
double get_harmonic_mean(const intermediate_stats& p);

double get_variance(const kernel& p, int analysis_param);
double get_variance(const kernel_propagate& p, int analysis_param);
double get_variance(const kernel_key_id& index, int analysis_param);
double get_variance(const intermediate_stats& p, int analysis_param);

double get_std_dev(const kernel& p, int analysis_param);
double get_std_dev(const kernel_propagate& p, int analysis_param);
double get_std_dev(const kernel_key_id& index, int analysis_param);
double get_std_dev(const intermediate_stats& p, int analysis_param);

int get_std_error_count(const kernel& p);
int get_std_error_count(const kernel_propagate& p);
int get_std_error_count(const kernel_key_id& index);
int get_std_error_count(const intermediate_stats& p);

double get_std_error(const comm_kernel_key& key, const kernel& p, int analysis_param);
double get_std_error(const comm_kernel_key& key, const kernel_propagate& p, int analysis_param);
double get_std_error(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param);
double get_std_error(const comm_kernel_key& key, const intermediate_stats& p, int analysis_param);
double get_std_error(const comp_kernel_key& key, const kernel& p, int analysis_param);
double get_std_error(const comp_kernel_key& key, const kernel_propagate& p, int analysis_param);
double get_std_error(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param);
double get_std_error(const comp_kernel_key& key, const intermediate_stats& p, int analysis_param);

double get_confidence_interval(const comm_kernel_key& key, const kernel& p, int analysis_param, double level = .95);
double get_confidence_interval(const comm_kernel_key& key, const kernel_propagate& p, int analysis_param, double level = .95);
double get_confidence_interval(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param, double level = .95);
double get_confidence_interval(const comm_kernel_key& key, const intermediate_stats& p, int analysis_param, double level = .95);
double get_confidence_interval(const comp_kernel_key& key, const kernel& p, int analysis_param, double level = .95);
double get_confidence_interval(const comp_kernel_key& key, const kernel_propagate& p, int analysis_param, double level = .95);
double get_confidence_interval(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param, double level = .95);
double get_confidence_interval(const comp_kernel_key& key, const intermediate_stats& p, int analysis_param, double level = .95);

bool is_steady(const kernel& p, int analysis_param);
bool is_steady(const kernel_key_id& index, int analysis_param);
bool is_steady(const intermediate_stats& p, int analysis_param);

double get_error_estimate(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param);
double get_error_estimate(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param);
double get_error_estimate(const comm_kernel_key& key, const kernel_propagate& p, int analysis_param);
double get_error_estimate(const comp_kernel_key& key, const kernel_propagate& p, int analysis_param);

bool steady_test(const comm_kernel_key& key, const kernel& p, int analysis_param);
bool steady_test(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param);
bool steady_test(const comp_kernel_key& key, const kernel& p, int analysis_param);
bool steady_test(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param);

void update_kernel_stats(kernel& p, int analysis_param, volatile double exec_time, double unit_count);
void update_kernel_stats(const kernel_key_id& index, int analysis_param, volatile double exec_time, double unit_count);
void update_kernel_stats(kernel& dest, const kernel& src, int analysis_param);
void update_kernel_stats(kernel_batch& batch, int analysis_param, volatile double exec_time, double unit_count);
void update_kernel_stats(kernel& dest, const kernel_batch& src, int analysis_param);
void update_kernel_stats(kernel_batch& dest, const kernel_batch& src, int analysis_param);
void update_kernel_stats(const kernel_key_id& index, const intermediate_stats& stats);

int should_schedule(const kernel& p);
int should_schedule(const kernel_key_id& index);

void set_kernel_state(kernel& p, bool schedule_decision);
void set_kernel_state(const kernel_key_id& index, bool schedule_decision);

void set_kernel_state_global(kernel& p, bool schedule_decision);
void set_kernel_state_global(const kernel_key_id& index, bool schedule_decision);

void merge_batches(std::vector<kernel_batch>& batches, int analysis_param);

void allocate(MPI_Comm comm);
void open_symbol(const char* symbol, double curtime);
void close_symbol(const char* symbol, double curtime);
void final_accumulate(MPI_Comm comm, double last_time);
void reset(bool schedule_kernels_override, bool force_steady_statistical_data_overide);
void clear(int tag_count, int* distribution_tags);
void finalize();

}
}
}

#endif /*CRITTER__DISCRETIZATION__UTIL__UTIL_H_*/
