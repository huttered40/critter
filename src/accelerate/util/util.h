#ifndef CRITTER__ACCELERATE__UTIL__UTIL_H_
#define CRITTER__ACCELERATE__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace accelerate{

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
  float num_scheduled_units;
  float num_local_scheduled_units;
  float total_exec_time;
  float total_local_exec_time;
  float M1,M2;
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
  float num_scheduled_units;
  float total_exec_time;
  float M1,M2;
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
  float num_scheduled_units;
  float num_non_scheduled_units;
  float num_local_scheduled_units;
  float M1,M2;
  float total_exec_time;
  float total_local_exec_time;
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
  float num_scheduled_units;
  float num_local_scheduled_units;
  float M1,M2;
  float total_exec_time;
  float total_local_exec_time;
};


// ****************************************************************************************************************************************************
struct intermediate_stats{
  intermediate_stats(const kernel_key_id& index, const std::vector<kernel_batch>& active_batches);
  void generate(const kernel& p, const std::vector<kernel_batch>& active_batches);

  int num_schedules;
  int num_local_schedules;
  float M1,M2;
  float num_scheduled_units;
  float num_local_scheduled_units;
  float total_exec_time;
  float total_local_exec_time;
};

// ****************************************************************************************************************************************************
extern std::ofstream stream,stream_comm_kernel,stream_comp_kernel,stream_tune,stream_reconstruct;
extern bool global_schedule_decision;
extern int tuning_delta,reset_distribution_mode,reset_state_mode;
extern int comm_sample_aggregation_mode,comm_state_aggregation_mode;
extern int comp_sample_aggregation_mode,comp_state_aggregation_mode;
extern int sample_constraint_mode,schedule_kernels,update_analysis;
extern int stop_criterion_mode,debug_iter_count;
extern int comp_kernel_transfer_id,comm_kernel_transfer_id;
extern int delay_state_update,collective_state_protocol;
extern float kernel_time_limit,kernel_error_limit,kernel_percentage_limit;
extern float* save_path_data;
extern MPI_Request save_prop_req;
extern volatile double comp_start_time;
extern size_t kernel_count_limit,mode_1_width,mode_2_width;
extern size_t num_cp_measures,num_pp_measures;
extern size_t num_vol_measures,num_tracker_cp_measures;
extern size_t num_tracker_pp_measures,num_tracker_vol_measures;
extern size_t cp_costs_size,pp_costs_size,vol_costs_size;
extern int internal_tag,internal_tag1,internal_tag2,internal_tag3,internal_tag4,internal_tag5;
extern bool is_first_iter;
extern MPI_Datatype kernel_type,batch_type;
extern std::map<comm_kernel_key,kernel_key_id> comm_kernel_map;
extern std::map<comp_kernel_key,kernel_key_id> comp_kernel_map;
extern std::map<comm_kernel_key,kernel> comm_kernel_save_map;
extern std::map<comp_kernel_key,kernel> comp_kernel_save_map;
extern std::map<comm_kernel_key,kernel> comm_kernel_ref_map;
extern std::map<comp_kernel_key,kernel> comp_kernel_ref_map;
extern std::vector<comm_kernel_key> active_comm_kernel_keys;
extern std::vector<comp_kernel_key> active_comp_kernel_keys;
extern std::vector<kernel> active_kernels;
extern sample_propagation_forest spf;
extern std::map<comm_kernel_key,std::vector<kernel_batch>> comm_batch_map;
extern std::map<comp_kernel_key,std::vector<kernel_batch>> comp_batch_map;
extern std::map<comm_kernel_key,std::vector<kernel>> comm_kernel_list;
extern std::map<comp_kernel_key,std::vector<kernel>> comp_kernel_list;
extern std::vector<float> intercept_overhead;
extern std::vector<float> global_intercept_overhead;
extern std::vector<float> global_comp_kernel_stats;
extern std::vector<float> global_comm_kernel_stats;
extern std::vector<float> local_comp_kernel_stats;
extern std::vector<float> local_comm_kernel_stats;
extern std::vector<float> save_comp_kernel_stats;
extern std::vector<float> save_comm_kernel_stats;
extern std::vector<float> cp_costs;
extern std::vector<float> cp_costs_foreign;
extern std::vector<float> max_pp_costs;
extern std::vector<float> vol_costs;
extern std::vector<float> cp_costs_ref;
extern std::vector<float> max_pp_costs_ref;
extern std::vector<float> vol_costs_ref;
extern std::vector<char> eager_pad;
extern std::map<std::string,std::vector<float>> save_info;

// ****************************************************************************************************************************************************
bool is_key_skipable(const comm_kernel_key& key);
bool is_key_skipable(const comp_kernel_key& key);

int get_skel_count(const comm_kernel_key& key);
int get_skel_count(const comp_kernel_key& key);

float get_estimate(const kernel& p, int analysis_param, float unit_count=1.);
float get_estimate(const kernel_propagate& p, int analysis_param, float unit_count=1.);
float get_estimate(const kernel_key_id& index, int analysis_param, float unit_count=1.);
float get_estimate(const kernel_key_id& index, const std::vector<kernel_batch>& active_batches, int analysis_param, float unit_count);
float get_estimate(const intermediate_stats& p, int analysis_param, float unit_count=1.);

float get_arithmetic_mean(const kernel& p);
float get_arithmetic_mean(const kernel_propagate& p);
float get_arithmetic_mean(const kernel_key_id& index);
float get_arithmetic_mean(const intermediate_stats& p);

float get_harmonic_mean(const kernel& p);
float get_harmonic_mean(const kernel_propagate& p);
float get_harmonic_mean(const kernel_key_id& index);
float get_harmonic_mean(const intermediate_stats& p);

float get_variance(const kernel& p, int analysis_param);
float get_variance(const kernel_propagate& p, int analysis_param);
float get_variance(const kernel_key_id& index, int analysis_param);
float get_variance(const intermediate_stats& p, int analysis_param);

float get_std_dev(const kernel& p, int analysis_param);
float get_std_dev(const kernel_propagate& p, int analysis_param);
float get_std_dev(const kernel_key_id& index, int analysis_param);
float get_std_dev(const intermediate_stats& p, int analysis_param);

int get_std_error_count(const kernel& p);
int get_std_error_count(const kernel_propagate& p);
int get_std_error_count(const kernel_key_id& index);
int get_std_error_count(const intermediate_stats& p);

float get_std_error(const comm_kernel_key& key, const kernel& p, int analysis_param);
float get_std_error(const comm_kernel_key& key, const kernel_propagate& p, int analysis_param);
float get_std_error(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param);
float get_std_error(const comm_kernel_key& key, const intermediate_stats& p, int analysis_param);
float get_std_error(const comp_kernel_key& key, const kernel& p, int analysis_param);
float get_std_error(const comp_kernel_key& key, const kernel_propagate& p, int analysis_param);
float get_std_error(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param);
float get_std_error(const comp_kernel_key& key, const intermediate_stats& p, int analysis_param);

float get_confidence_interval(const comm_kernel_key& key, const kernel& p, int analysis_param, float level = .95);
float get_confidence_interval(const comm_kernel_key& key, const kernel_propagate& p, int analysis_param, float level = .95);
float get_confidence_interval(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param, float level = .95);
float get_confidence_interval(const comm_kernel_key& key, const intermediate_stats& p, int analysis_param, float level = .95);
float get_confidence_interval(const comp_kernel_key& key, const kernel& p, int analysis_param, float level = .95);
float get_confidence_interval(const comp_kernel_key& key, const kernel_propagate& p, int analysis_param, float level = .95);
float get_confidence_interval(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param, float level = .95);
float get_confidence_interval(const comp_kernel_key& key, const intermediate_stats& p, int analysis_param, float level = .95);

bool is_steady(const kernel& p, int analysis_param);
bool is_steady(const kernel_key_id& index, int analysis_param);
bool is_steady(const intermediate_stats& p, int analysis_param);

float get_error_estimate(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param);
float get_error_estimate(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param);
float get_error_estimate(const comm_kernel_key& key, const kernel_propagate& p, int analysis_param);
float get_error_estimate(const comp_kernel_key& key, const kernel_propagate& p, int analysis_param);

bool steady_test(const comm_kernel_key& key, const kernel& p, int analysis_param);
bool steady_test(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param);
bool steady_test(const comp_kernel_key& key, const kernel& p, int analysis_param);
bool steady_test(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param);

void update_kernel_stats(kernel& p, int analysis_param, volatile float exec_time, float unit_count);
void update_kernel_stats(const kernel_key_id& index, int analysis_param, volatile float exec_time, float unit_count);
void update_kernel_stats(kernel& dest, const kernel& src, int analysis_param);
void update_kernel_stats(kernel_batch& batch, int analysis_param, volatile float exec_time, float unit_count);
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
void init_symbol(std::vector<std::string>& symbols);
void open_symbol(const char* symbol, double curtime);
void close_symbol(const char* symbol, double curtime);
void final_accumulate(MPI_Comm comm, double last_time);
void reset(bool schedule_kernels_override, bool force_steady_statistical_data_overide);
void clear(int tag_count=0, int* distribution_tags=nullptr);
void reference_initiate();
void reference_transfer();
void finalize();

}
}
}

#endif /*CRITTER__ACCELERATE__UTIL__UTIL_H_*/
