#ifndef CRITTER__DISCRETIZATION__UTIL__UTIL_H_
#define CRITTER__DISCRETIZATION__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

struct pattern_batch;

// ****************************************************************************************************************************************************
struct comm_channel_node{
  comm_channel_node();

  int hash_tag;
  int tag;
  int frequency;
  int offset;
  std::vector<std::pair<int,int>> id;
  comm_channel_node* parent;
  std::vector<std::vector<comm_channel_node*>> children;
};

struct sample_propagation_tree{
  comm_channel_node* root;
};

struct sample_propagation_forest{
  sample_propagation_forest();
  ~sample_propagation_forest();

  void generate_span(comm_channel_node* node, std::vector<std::pair<int,int>>& perm_tuples);
  int translate_rank(MPI_Comm comm, int rank);
  void insert_node(comm_channel_node* tree_node);
  void clear_info();
  void fill_ancestors(comm_channel_node* tree_node, pattern_batch& batch);
  void fill_descendants(comm_channel_node* tree_node, pattern_batch& batch);

  sample_propagation_tree* tree;
private:
  void delete_tree(comm_channel_node*& tree_root);
  void clear_tree_info(comm_channel_node* tree_root);
  void generate_sibling_perm(std::vector<std::pair<int,int>>& static_info, std::vector<std::pair<int,int>>& gen_info, std::vector<std::pair<int,int>>& save_info, int level, bool& valid_siblings);
  void generate_partition_perm(std::vector<std::pair<int,int>>& static_info, std::vector<std::pair<int,int>>& gen_info, int level, bool& valid_partition,
                               int parent_max_span, int parent_min_stride);
  bool is_child(comm_channel_node* tree_node, comm_channel_node* node);
  bool are_siblings(comm_channel_node* node, int subtree_idx, std::vector<int>& skip_indices);
  bool partition_test(comm_channel_node* parent, int subtree_idx);
  void find_parent(comm_channel_node* tree_root, comm_channel_node* tree_node, comm_channel_node*& parent);
  int span(std::pair<int,int>& id);
};

// ****************************************************************************************************************************************************
struct pattern{
  // If I leverage the kurtosis, I will have to utilize the arithmetic mean.
  //   Note that I'd rather utilize the geometric mean, but I'm not sure how to convert this algorithm
  //     to handle that.
  pattern();
  pattern(const pattern& _copy);
  pattern& operator=(const pattern& _copy);

  int steady_state;
  int global_steady_state;
  int num_schedules;
  int num_non_schedules;
  int num_propagations;
  int num_non_propagations;
  double num_scheduled_units;
  double num_non_scheduled_units;
  double M1,M2;
  double total_exec_time;
};

// ****************************************************************************************************************************************************
struct pattern_batch{
  // If I leverage the kurtosis, I will have to utilize the arithmetic mean.
  //   Note that I'd rather utilize the geometric mean, but I'm not sure how to convert this algorithm
  //     to handle that.
  pattern_batch(comm_channel_node* node = nullptr);
  pattern_batch(const pattern_batch& _copy);
  pattern_batch& operator=(const pattern_batch& _copy);

  int channel_count;
  int num_schedules;
  int open_channel_count;
  double num_scheduled_units;
  double total_exec_time;
  double M1,M2;
  std::set<comm_channel_node*> closed_channels;
};

// ****************************************************************************************************************************************************
struct idle_pattern{
  // If I leverage the kurtosis, I will have to utilize the arithmetic mean.
  //   Note that I'd rather utilize the geometric mean, but I'm not sure how to convert this algorithm
  //     to handle that.
  idle_pattern();
  idle_pattern(const idle_pattern& _copy);
  idle_pattern& operator=(const idle_pattern& _copy);

  int num_schedules;
  int num_non_schedules;
  double M1,M2;
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
  double M1,M2;
  double num_scheduled_units;
  double total_exec_time;
};

// ****************************************************************************************************************************************************
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
extern int comm_channel_tag_count;
extern std::map<comm_pattern_key,pattern_key_id> comm_pattern_map;
extern std::map<comp_pattern_key,pattern_key_id> comp_pattern_map;
extern std::vector<comm_pattern_key> steady_state_comm_pattern_keys;
extern std::vector<comm_pattern_key> active_comm_pattern_keys;
extern std::vector<comp_pattern_key> steady_state_comp_pattern_keys;
extern std::vector<comp_pattern_key> active_comp_pattern_keys;
extern std::vector<pattern> steady_state_patterns;
extern std::vector<pattern> active_patterns;
extern std::map<std::pair<comm_pattern_key,comm_pattern_key>,idle_pattern> comm_pattern_pair_map;
extern sample_propagation_forest spf;
extern std::map<MPI_Comm,comm_channel_node*> comm_channel_map;
extern std::map<int,comm_channel_node*> p2p_channel_map;
extern std::map<comm_pattern_key,std::vector<pattern_batch>> comm_batch_map;
extern std::map<comp_pattern_key,std::vector<pattern_batch>> comp_batch_map;
extern std::vector<comm_channel_node*> intermediate_channels;
extern std::map<comm_pattern_key,bool> p2p_global_state_override;

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

bool steady_test(const comm_pattern_key& key, const pattern& p, int analysis_param);
bool steady_test(const comm_pattern_key& key, const pattern_key_id& index, int analysis_param);
bool steady_test(const comp_pattern_key& key, const pattern& p, int analysis_param);
bool steady_test(const comp_pattern_key& key, const pattern_key_id& index, int analysis_param);

void update_kernel_stats(pattern& p, int analysis_param, volatile double exec_time, double unit_count);
void update_kernel_stats(const pattern_key_id& index, int analysis_param, volatile double exec_time, double unit_count);
void update_kernel_stats(pattern& dest, const pattern& src, int analysis_param);
void update_kernel_stats(idle_pattern& p, bool is_global_steady_state, volatile double exec_time);
void update_kernel_stats(pattern_batch& batch, int analysis_param, volatile double exec_time, double unit_count);
void update_kernel_stats(const pattern_key_id& index, const intermediate_stats& stats);

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
void reset_frequencies();
void clear();
void finalize();

}
}
}

#endif /*CRITTER__DISCRETIZATION__UTIL__UTIL_H_*/
