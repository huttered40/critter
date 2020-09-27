#ifndef CRITTER__DISCRETIZATION__UTIL__UTIL_H_
#define CRITTER__DISCRETIZATION__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

struct pattern_batch;

// ****************************************************************************************************************************************************
struct channel{
  channel();
  static std::vector<std::pair<int,int>> generate_tuple(std::vector<int>& ranks, int new_comm_size);
  static void contract_tuple(std::vector<std::pair<int,int>>& tuple_list);
  static int enumerate_tuple(channel* node, std::vector<int>& process_list);
  static int duplicate_process_count(std::vector<int>& process_list);
  static int translate_rank(MPI_Comm comm, int rank);
  static bool verify_ancestor_relation(channel* comm1, channel* comm2);
  static bool verify_sibling_relation(channel* comm1, channel* comm2);
  static int span(std::pair<int,int>& id);
  static std::string generate_tuple_string(channel* comm);
 
  int offset;
  int local_hash_tag;
  int global_hash_tag;
  std::vector<std::pair<int,int>> id;
};

struct aggregate_channel : public channel{
  aggregate_channel(std::vector<std::pair<int,int>>& tuple_list, int local_hash, int global_hash, int offset, int channel_size);
  static std::string generate_hash_history(aggregate_channel* comm);

  bool is_final;
  int num_channels;
  std::set<int> channels;
};

struct solo_channel : public channel{
  solo_channel();
  static bool verify_sibling_relation(solo_channel* node, int subtree_idx, std::vector<int>& skip_indices);

  int tag;
  int frequency;
  solo_channel* parent;
  std::vector<std::vector<solo_channel*>> children;
};

struct sample_propagation_tree{
 solo_channel* root;
};

struct sample_propagation_forest{
  sample_propagation_forest();
  ~sample_propagation_forest();

  //void generate_span(channel* node, std::vector<std::pair<int,int>>& perm_tuples);
  void insert_node(solo_channel* node);
  void clear_info();
  void fill_ancestors(solo_channel* node, pattern_batch& batch);
  void fill_descendants(solo_channel* node, pattern_batch& batch);

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
struct pattern_batch{
  // If I leverage the kurtosis, I will have to utilize the arithmetic mean.
  //   Note that I'd rather utilize the geometric mean, but I'm not sure how to convert this algorithm
  //     to handle that.
  pattern_batch(channel* node = nullptr);
  pattern_batch(const pattern_batch& _copy);
  pattern_batch& operator=(const pattern_batch& _copy);

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
struct pattern_batch_propagate{
  // If I leverage the kurtosis, I will have to utilize the arithmetic mean.
  //   Note that I'd rather utilize the geometric mean, but I'm not sure how to convert this algorithm
  //     to handle that.
  pattern_batch_propagate(){};
  pattern_batch_propagate(const pattern_batch& _copy);
  pattern_batch_propagate(const pattern_batch_propagate& _copy);
  pattern_batch_propagate& operator=(const pattern_batch_propagate& _copy);

  int hash_id;
  int channel_count;
  int num_schedules;
  double num_scheduled_units;
  double total_exec_time;
  double M1,M2;
};

// ****************************************************************************************************************************************************
struct pattern{
  // If I leverage the kurtosis, I will have to utilize the arithmetic mean.
  //   Note that I'd rather utilize the geometric mean, but I'm not sure how to convert this algorithm
  //     to handle that.
  pattern();
  pattern(const pattern_batch& _copy);
  pattern(const pattern& _copy);
  pattern& operator=(const pattern& _copy);

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
};

// ****************************************************************************************************************************************************
struct pattern_key_id{

  pattern_key_id(bool _is_active=false, int _key_index=0, int _val_index=0, bool _is_updated=false);
  pattern_key_id(const pattern_key_id& _copy);
  pattern_key_id& operator=(const pattern_key_id& _copy);

  // Active just means its still being propogated. It acts as a switch betweeh steady_state arrays and active arrays
  bool is_active;
  bool is_updated;
  int key_index;
  int val_index;
};

// ****************************************************************************************************************************************************
struct intermediate_stats{
  intermediate_stats(const pattern_key_id& index, const std::vector<pattern_batch>& active_batches);
  void generate(const pattern& p, const std::vector<pattern_batch>& active_batches);

  int num_schedules;
  int num_local_schedules;
  double M1,M2;
  double num_scheduled_units;
  double num_local_scheduled_units;
  double total_exec_time;
  double total_local_exec_time;
};

// ****************************************************************************************************************************************************
extern int is_optimized;
extern int aggregation_mode;
extern int schedule_kernels;
extern int update_analysis;
extern int autotuning_test_id;
extern MPI_Datatype comm_pattern_key_type;
extern MPI_Datatype comp_pattern_key_type;
extern MPI_Datatype pattern_type;
extern MPI_Datatype batch_type;
extern size_t pattern_count_limit;
extern double pattern_time_limit;
extern double pattern_error_limit;
extern int communicator_count;
extern std::map<comm_pattern_key,pattern_key_id> comm_pattern_map;
extern std::map<comp_pattern_key,pattern_key_id> comp_pattern_map;
extern std::vector<comm_pattern_key> steady_state_comm_pattern_keys;
extern std::vector<comm_pattern_key> active_comm_pattern_keys;
extern std::vector<comp_pattern_key> steady_state_comp_pattern_keys;
extern std::vector<comp_pattern_key> active_comp_pattern_keys;
extern std::vector<pattern> steady_state_patterns;
extern std::vector<pattern> active_patterns;
extern sample_propagation_forest spf;
extern std::map<MPI_Comm,solo_channel*> comm_channel_map;
extern std::map<int,solo_channel*> p2p_channel_map;
extern std::map<comm_pattern_key,std::vector<pattern_batch>> comm_batch_map;
extern std::map<comp_pattern_key,std::vector<pattern_batch>> comp_batch_map;
extern std::map<comm_pattern_key,bool> p2p_global_state_override;
extern std::map<int,aggregate_channel*> aggregate_channel_map;
extern std::ofstream stream_overhead,stream_tune,stream_reconstruct;

// ****************************************************************************************************************************************************
bool is_key_skipable(const comm_pattern_key& key);
bool is_key_skipable(const comp_pattern_key& key);

double get_estimate(const pattern& p, int analysis_param, double unit_count=1.);
double get_estimate(const pattern_key_id& index, int analysis_param, double unit_count=1.);
double get_estimate(const pattern_key_id& index, const std::vector<pattern_batch>& active_batches, int analysis_param, double unit_count);
double get_estimate(const intermediate_stats& p, int analysis_param, double unit_count=1.);

double get_arithmetic_mean(const pattern& p);
double get_arithmetic_mean(const pattern_key_id& index);
double get_arithmetic_mean(const intermediate_stats& p);

double get_harmonic_mean(const pattern& p);
double get_harmonic_mean(const pattern_key_id& index);
double get_harmonic_mean(const intermediate_stats& p);

double get_variance(const pattern& p, int analysis_param);
double get_variance(const pattern_key_id& index, int analysis_param);
double get_variance(const intermediate_stats& p, int analysis_param);

double get_std_dev(const pattern& p, int analysis_param);
double get_std_dev(const pattern_key_id& index, int analysis_param);
double get_std_dev(const intermediate_stats& p, int analysis_param);

double get_std_error(const pattern& p, int analysis_param);
double get_std_error(const pattern_key_id& index, int analysis_param);
double get_std_error(const intermediate_stats& p, int analysis_param);

double get_confidence_interval(const pattern& p, int analysis_param, double level = .95);
double get_confidence_interval(const pattern_key_id& index, int analysis_param, double level = .95);
double get_confidence_interval(const intermediate_stats& p, int analysis_param, double level = .95);

bool is_steady(const pattern& p, int analysis_param);
bool is_steady(const pattern_key_id& index, int analysis_param);
bool is_steady(const intermediate_stats& p, int analysis_param);

double get_error_estimate(const comm_pattern_key& key, const pattern_key_id& index, int analysis_param);
double get_error_estimate(const comp_pattern_key& key, const pattern_key_id& index, int analysis_param);

bool steady_test(const comm_pattern_key& key, const pattern& p, int analysis_param);
bool steady_test(const comm_pattern_key& key, const pattern_key_id& index, int analysis_param);
bool steady_test(const comp_pattern_key& key, const pattern& p, int analysis_param);
bool steady_test(const comp_pattern_key& key, const pattern_key_id& index, int analysis_param);

void update_kernel_stats(pattern& p, int analysis_param, volatile double exec_time, double unit_count);
void update_kernel_stats(const pattern_key_id& index, int analysis_param, volatile double exec_time, double unit_count);
void update_kernel_stats(pattern& dest, const pattern& src, int analysis_param);
void update_kernel_stats(pattern_batch& batch, int analysis_param, volatile double exec_time, double unit_count);
void update_kernel_stats(pattern& dest, const pattern_batch& src, int analysis_param);
void update_kernel_stats(pattern_batch& dest, const pattern_batch& src, int analysis_param);
void update_kernel_stats(const pattern_key_id& index, const intermediate_stats& stats);

int should_schedule(const pattern& p);
int should_schedule(const pattern_key_id& index);

int should_schedule_global(const pattern& p);
int should_schedule_global(const pattern_key_id& index);

void set_kernel_state(pattern& p, bool schedule_decision);
void set_kernel_state(const pattern_key_id& index, bool schedule_decision);

void set_kernel_state_global(pattern& p, bool schedule_decision);
void set_kernel_state_global(const pattern_key_id& index, bool schedule_decision);

void merge_batches(std::vector<pattern_batch>& batches, int analysis_param);

void allocate(MPI_Comm comm);
void open_symbol(const char* symbol, double curtime);
void close_symbol(const char* symbol, double curtime);
void final_accumulate(MPI_Comm comm, double last_time);
void reset(bool schedule_kernels_override, bool force_steady_statistical_data_overide);
void reset_frequencies();
void clear();
void finalize();

}
}
}

#endif /*CRITTER__DISCRETIZATION__UTIL__UTIL_H_*/
