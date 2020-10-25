#include <limits.h>

#include "util.h"
#include "../../skeletonization/util/util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace discretization{

int tuning_delta;
int comm_sample_aggregation_mode;
int comm_state_aggregation_mode;
int comp_sample_aggregation_mode;
int comp_state_aggregation_mode;
int sample_constraint_mode;
//int sample_reset_mode;
int schedule_kernels;
int update_analysis;
MPI_Datatype kernel_type;
MPI_Datatype batch_type;
size_t kernel_count_limit;
double kernel_time_limit;
double kernel_error_limit;
int comp_kernel_transfer_id;
int comm_kernel_transfer_id;
int comp_kernel_buffer_id;
std::map<comm_kernel_key,kernel_key_id> comm_kernel_map;
std::map<comp_kernel_key,kernel_key_id> comp_kernel_map;
std::vector<comm_kernel_key> active_comm_kernel_keys;
std::vector<comp_kernel_key> active_comp_kernel_keys;
std::vector<kernel> active_kernels;
sample_propagation_forest spf;
std::map<comm_kernel_key,std::vector<kernel_batch>> comm_batch_map;
std::map<comp_kernel_key,std::vector<kernel_batch>> comp_batch_map;
std::set<intptr_t> skip_ptr_set;

std::ofstream stream,stream_tune,stream_reconstruct;
std::vector<double> intercept_overhead;
std::vector<double> global_intercept_overhead;
std::vector<double> global_comp_kernel_stats;
std::vector<double> global_comm_kernel_stats;
std::vector<double> local_comp_kernel_stats;
std::vector<double> local_comm_kernel_stats;
std::vector<double> save_comp_kernel_stats;
std::vector<double> save_comm_kernel_stats;
size_t num_critical_path_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime, CompTime, RunTime
size_t num_per_process_measures;		// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, RunTime
size_t num_volume_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, RunTime
size_t num_tracker_critical_path_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime
size_t num_tracker_per_process_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime
size_t num_tracker_volume_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime
size_t critical_path_costs_size;
size_t per_process_costs_size;
size_t volume_costs_size;
std::vector<double> critical_path_costs;
std::vector<double> max_per_process_costs;
std::vector<double> volume_costs;
std::vector<double_int> info_sender;
std::vector<double_int> info_receiver;
std::vector<char> eager_pad;
double comp_start_time;
std::map<MPI_Request,std::pair<bool,int>> internal_comm_info;
std::map<MPI_Request,std::pair<MPI_Comm,int>> internal_comm_comm;
std::map<MPI_Request,std::pair<double,double>> internal_comm_data;
std::vector<std::pair<double*,int>> internal_comm_prop;
std::vector<MPI_Request> internal_comm_prop_req;
std::vector<int*> internal_timer_prop_int;
std::vector<double*> internal_timer_prop_double;
std::vector<double_int*> internal_timer_prop_double_int;
std::vector<char*> internal_timer_prop_char;
std::vector<MPI_Request> internal_timer_prop_req;
std::vector<bool> decisions;
std::map<std::string,std::vector<double>> save_info;
std::vector<double> new_cs;
size_t mode_1_width;
size_t mode_2_width;
int internal_tag;
int internal_tag1;
int internal_tag2;
int internal_tag3;
int internal_tag4;
int internal_tag5;
bool is_first_iter;

// ****************************************************************************************************************************************************
kernel::kernel(channel* node){
  this->hash_id = 0;
  this->total_exec_time = 0;
  this->total_local_exec_time = 0;
  this->steady_state=0;
  this->global_steady_state=0;
  this->num_schedules = 0;
  this->num_non_schedules = 0;
  this->num_scheduled_units = 0;
  this->num_local_scheduled_units = 0;
  this->num_non_scheduled_units = 0;
  this->num_local_schedules=0;
  this->M1=0; this->M2=0;
  if (node != nullptr){
    this->hash_id = node->global_hash_tag;
    this->registered_channels.insert(node);
  }
}
kernel::kernel(const kernel_batch& _copy){
  this->hash_id = _copy.hash_id;
  this->total_exec_time = _copy.total_exec_time;
  this->total_local_exec_time = _copy.total_local_exec_time;
  this->steady_state=0;
  this->global_steady_state=0;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = 0;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_local_scheduled_units = _copy.num_local_scheduled_units;
  this->num_non_scheduled_units = 0;
  this->num_local_schedules=_copy.num_local_schedules;
  this->M1=_copy.M1;
  this->M2=_copy.M2;
  this->registered_channels = _copy.registered_channels;
}
kernel::kernel(const kernel& _copy){
  this->hash_id = _copy.hash_id;
  this->total_exec_time = _copy.total_exec_time;
  this->total_local_exec_time = _copy.total_local_exec_time;
  this->steady_state = _copy.steady_state;
  this->global_steady_state = _copy.global_steady_state;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = _copy.num_non_schedules;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_local_scheduled_units = _copy.num_local_scheduled_units;
  this->num_non_scheduled_units = _copy.num_non_scheduled_units;
  this->num_local_schedules=_copy.num_local_schedules;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  this->registered_channels = _copy.registered_channels;
}
kernel& kernel::operator=(const kernel& _copy){
  this->hash_id = _copy.hash_id;
  this->total_exec_time = _copy.total_exec_time;
  this->total_local_exec_time = _copy.total_local_exec_time;
  this->steady_state = _copy.steady_state;
  this->global_steady_state = _copy.global_steady_state;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = _copy.num_non_schedules;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_local_scheduled_units = _copy.num_local_scheduled_units;
  this->num_non_scheduled_units = _copy.num_non_scheduled_units;
  this->num_local_schedules=_copy.num_local_schedules;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  this->registered_channels = _copy.registered_channels;
  return *this;
}
kernel::kernel(const kernel_propagate& _copy){
  this->hash_id = _copy.hash_id;
  this->total_exec_time = _copy.total_exec_time;
  this->total_local_exec_time = _copy.total_local_exec_time;
  this->steady_state = false;
  this->global_steady_state = false;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = 0;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_local_scheduled_units = _copy.num_local_scheduled_units;
  this->num_non_scheduled_units = 0;
  this->num_local_schedules=_copy.num_local_schedules;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  if (_copy.hash_id != 0){
    for (auto& it : aggregate_channel_map[_copy.hash_id]->channels){
      this->registered_channels.insert(comm_channel_map[it]);
    }
  }
}
kernel& kernel::operator=(const kernel_propagate& _copy){
  this->hash_id = _copy.hash_id;
  this->total_exec_time = _copy.total_exec_time;
  this->total_local_exec_time = _copy.total_local_exec_time;
  this->steady_state = false;
  this->global_steady_state = false;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = 0;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_local_scheduled_units = _copy.num_local_scheduled_units;
  this->num_non_scheduled_units = 0;
  this->num_local_schedules=_copy.num_local_schedules;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  if (_copy.hash_id != 0){
    for (auto& it : aggregate_channel_map[_copy.hash_id]->channels){
      this->registered_channels.insert(comm_channel_map[it]);
    }
  }
  return *this;
}
void kernel::reset(){
  // No reason to reset distribution members ('M1','M2','steady_state','global_steady_state')
  // Do not reset the schedule_count members because those are needed to determine sample variance,
  //   which builds across schedules.
  this->num_local_schedules=0;
  this->total_local_exec_time = 0;
  this->num_local_scheduled_units = 0;
}
void kernel::clear_distribution(){
  // All other members would have already been reset (via 'reset' above)
  this->steady_state = false;
  this->global_steady_state = false;
  this->total_exec_time = 0;
  this->num_scheduled_units = 0;
  this->num_schedules = 0;
  this->M1=0; this->M2=0;
}

// ****************************************************************************************************************************************************
kernel_propagate::kernel_propagate(const kernel& _copy){
  this->hash_id = _copy.hash_id;
  this->total_exec_time = _copy.total_exec_time;
  this->total_local_exec_time = _copy.total_local_exec_time;
  this->num_schedules = _copy.num_schedules;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_local_scheduled_units = _copy.num_local_scheduled_units;
  this->num_local_schedules=_copy.num_local_schedules;
  this->M1=_copy.M1;
  this->M2=_copy.M2;
}
kernel_propagate::kernel_propagate(const kernel_propagate& _copy){
  this->hash_id = _copy.hash_id;
  this->total_exec_time = _copy.total_exec_time;
  this->total_local_exec_time = _copy.total_local_exec_time;
  this->num_schedules = _copy.num_schedules;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_local_scheduled_units = _copy.num_local_scheduled_units;
  this->num_local_schedules=_copy.num_local_schedules;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
}
kernel_propagate& kernel_propagate::operator=(const kernel_propagate& _copy){
  this->hash_id = _copy.hash_id;
  this->total_exec_time = _copy.total_exec_time;
  this->total_local_exec_time = _copy.total_local_exec_time;
  this->num_schedules = _copy.num_schedules;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_local_scheduled_units = _copy.num_local_scheduled_units;
  this->num_local_schedules=_copy.num_local_schedules;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  return *this;
}

// ****************************************************************************************************************************************************
kernel_batch::kernel_batch(channel* node){
  this->hash_id = 0;
  this->total_exec_time = 0;
  this->total_local_exec_time = 0;
  this->num_scheduled_units = 0;
  this->num_local_scheduled_units = 0;
  this->channel_count=0;
  this->num_schedules = 0;
  this->num_local_schedules = 0;
  this->M1=0; this->M2=0;
  if (node != nullptr){
    this->hash_id = node->global_hash_tag;
    this->registered_channels.insert(node);
  }
}

kernel_batch::kernel_batch(const kernel_batch& _copy){
  this->hash_id = _copy.hash_id;
  this->total_exec_time = _copy.total_exec_time;
  this->total_local_exec_time = _copy.total_local_exec_time;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_local_scheduled_units = _copy.num_local_scheduled_units;
  this->channel_count = _copy.channel_count;
  this->num_schedules = _copy.num_schedules;
  this->num_local_schedules = _copy.num_local_schedules;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  this->registered_channels = _copy.registered_channels;
}

kernel_batch& kernel_batch::operator=(const kernel_batch& _copy){
  this->hash_id = _copy.hash_id;
  this->total_exec_time = _copy.total_exec_time;
  this->total_local_exec_time = _copy.total_local_exec_time;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_local_scheduled_units = _copy.num_local_scheduled_units;
  this->channel_count = _copy.channel_count;
  this->num_schedules = _copy.num_schedules;
  this->num_local_schedules = _copy.num_local_schedules;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  this->registered_channels = _copy.registered_channels;
  return *this;
}

// ****************************************************************************************************************************************************
kernel_batch_propagate::kernel_batch_propagate(const kernel_batch& _copy){
  this->hash_id = _copy.hash_id;
  this->total_exec_time = _copy.total_exec_time;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->channel_count=_copy.channel_count;
  this->num_schedules = _copy.num_schedules;
  this->M1=_copy.M1;
  this->M2=_copy.M2;
}

kernel_batch_propagate::kernel_batch_propagate(const kernel_batch_propagate& _copy){
  this->hash_id = _copy.hash_id;
  this->total_exec_time = _copy.total_exec_time;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->channel_count = _copy.channel_count;
  this->num_schedules = _copy.num_schedules;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
}

kernel_batch_propagate& kernel_batch_propagate::operator=(const kernel_batch_propagate& _copy){
  this->hash_id = _copy.hash_id;
  this->total_exec_time = _copy.total_exec_time;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->channel_count = _copy.channel_count;
  this->num_schedules = _copy.num_schedules;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  return *this;
}

// ****************************************************************************************************************************************************

/*
void sample_propagation_forest::generate_span(channel* node, std::vector<std::pair<int,int>>& perm_tuples){
  // Assumed that perm_tuples define a permutation of channels that are communicator siblings
  //   Thus, 'perm_tuples' will be in sorted order of stride
  // Only need to fill out the 'id' member
  node->id.push_back(perm_tuples[0]);
  int index=0;
  for (int i=1; i<perm_tuples.size(); i++){
    if (node->id[index].first*node->id[index].second == perm_tuples[i].second){
      node->id[index].first *= perm_tuples[i].first;
      node->id[index].second = std::min(node->id[index].second,perm_tuples[i].second);
    }
    else if (node->id[index].first*node->id[index].second < perm_tuples[i].second){
      index++;
      node->id.push_back(perm_tuples[i]);
    }
    else { assert(0); }
  }
}
*/
/*
void sample_propagation_forest::generate_sibling_perm(std::vector<std::pair<int,int>>& static_info, std::vector<std::pair<int,int>>& gen_info,
                                                      std::vector<std::pair<int,int>>& save_info, int level, bool& valid_siblings){
  // static_info will shrink as tuples are transfered into gen_info. That means when 
  if (static_info.size()==0){
    save_info = gen_info;
    valid_siblings=true;
    return;
  }
  // At position 'level' in permutation, lets try all remaining possibilities via iterating over what remains of static_info
  for (auto i=0; i<static_info.size(); i++){
    // check if valid BEFORE recursing, except if at first level (no pruning at that level, all possibilities still valid)
    if ((level==0) || (static_info[i].second == (gen_info[level-1].first*gen_info[level-1].second))){
      gen_info.push_back(static_info[i]);
      if (i==static_info.size()-1){
        static_info.pop_back();
        this->generate_sibling_perm(static_info,gen_info,save_info,level+1,valid_siblings);
        static_info.push_back(gen_info[gen_info.size()-1]);
      } else{
        // swap with last entry and then pop
        auto temp = static_info[static_info.size()-1];
        static_info[static_info.size()-1] = static_info[i];
        static_info[i] = temp;
        static_info.pop_back();
        this->generate_sibling_perm(static_info,gen_info,save_info,level+1,valid_siblings);
        static_info.push_back(temp);
        static_info[i] = gen_info[gen_info.size()-1];
      }
      gen_info.pop_back();
    }
  }
}
*/
/*
void sample_propagation_forest::generate_partition_perm(std::vector<std::pair<int,int>>& static_info, std::vector<std::pair<int,int>>& gen_info, int level,
                                                        bool& valid_partition, int parent_max_span, int parent_min_stride){
  // static_info will shrink as tuples are transfered into gen_info. That means when 
  if (static_info.size()==0){
    if ((level>0) && (gen_info[0].second == parent_min_stride) && (gen_info[level-1].first*gen_info[level-1].second == parent_max_span)){ valid_partition=true; }
    return;
  }
  // At position 'level' in permutation, lets try all remaining possibilities via iterating over what remains of static_info
  int static_info_size=static_info.size();
  for (auto i=0; i<static_info_size; i++){
    // check if valid BEFORE recursing, except if at first level (no pruning at that level, all possibilities still valid)
    if ((level==0) || (static_info[i].second == (gen_info[level-1].first*gen_info[level-1].second))){
      //if ((level==0) && static_info[i].second!= 1) continue;// constrain initial stride in permutation
      gen_info.push_back(static_info[i]);
      if (i==static_info.size()-1){
        static_info.pop_back();
        this->generate_partition_perm(static_info,gen_info,level+1,valid_partition,parent_max_span,parent_min_stride);
        static_info.push_back(gen_info[gen_info.size()-1]);
      } else{
        // swap with last entry and then pop
        auto temp = static_info[static_info.size()-1];
        static_info[static_info.size()-1] = static_info[i];
        static_info[i] = temp;
        static_info.pop_back();
        this->generate_partition_perm(static_info,gen_info,level+1,valid_partition,parent_max_span,parent_min_stride);
        static_info.push_back(temp);
        static_info[i] = gen_info[gen_info.size()-1];
      }
      gen_info.pop_back();
    }
  }
}
*/
/*
bool sample_propagation_forest::partition_test(channel* parent, int subtree_idx){
  // Perform recursive permutation generation to identify if a permutation of tuples among siblings is valid
  // Return true if parent's children are valid siblings
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  std::vector<std::pair<int,int>> static_info;
  for (auto i=0; i<parent->children[subtree_idx].size(); i++){
    for (auto j=0; j<parent->children[subtree_idx][i]->id.size(); j++){
      if (parent->children[subtree_idx][i]->tag < ((-1)*world_size)){
        static_info.push_back(parent->children[subtree_idx][i]->id[j]);
      }
    }
  }
  std::vector<std::pair<int,int>> gen_info;
  bool valid_partition=false;
  auto parent_max_span = parent->id[parent->id.size()-1].first * parent->id[parent->id.size()-1].second;
  auto parent_min_stride = parent->id[0].second;
  generate_partition_perm(static_info,gen_info,0,valid_partition,parent_max_span,parent_min_stride);
  return valid_partition;
}
*/
void sample_propagation_forest::find_parent(solo_channel* tree_root, solo_channel* tree_node, solo_channel*& parent){
  if (tree_root==nullptr) return;
  for (auto i=0; i<tree_root->children.size(); i++){
    for (auto j=0; j<tree_root->children[i].size(); j++){
      this->find_parent(tree_root->children[i][j],tree_node,parent);// Cannot be nullptrs. Nullptr children mean the children member is empty
    }
  }
  if ((parent==nullptr) && (channel::verify_ancestor_relation(tree_node,tree_root))){
    parent = tree_root;
  }
  return;
}
void sample_propagation_forest::fill_ancestors(solo_channel* node, kernel_batch& batch){
/*
  if (node==nullptr) return;
  batch.closed_channels.insert(node);
  this->fill_ancestors(node->parent,batch);
*/
}
void sample_propagation_forest::fill_descendants(solo_channel* node, kernel_batch& batch){
/*
  if (node==nullptr) return;
  batch.closed_channels.insert(node);
  for (auto i=0; i<node->children.size(); i++){
    for (auto j=0; j<node->children[i].size(); j++){
      this->fill_descendants(node->children[i][j],batch);
    }
  }
*/
}
void sample_propagation_forest::clear_tree_info(solo_channel* tree_root){
  if (tree_root==nullptr) return;
  for (auto i=0; i<tree_root->children.size(); i++){
    for (auto j=0; j<tree_root->children[i].size(); j++){
      this->clear_tree_info(tree_root->children[i][j]);// Cannot be nullptrs. Nullptr children mean the children member is empty
    }
  }
  tree_root->frequency=0;
  return;
}
void sample_propagation_forest::delete_tree(solo_channel*& tree_root){
  if (tree_root==nullptr) return;
  for (auto i=0; i<tree_root->children.size(); i++){
    for (auto j=0; j<tree_root->children[i].size(); j++){
      this->delete_tree(tree_root->children[i][j]);// Cannot be nullptrs. Nullptr children mean the children member is empty
    }
  }
  free(tree_root);
  tree_root=nullptr;
  return;
}
sample_propagation_forest::sample_propagation_forest(){ this->tree=nullptr; }
sample_propagation_forest::~sample_propagation_forest(){
  if (this->tree == nullptr) return;
  this->delete_tree(this->tree->root);
  free(this->tree); this->tree = nullptr;
}
void sample_propagation_forest::clear_info(){
  this->clear_tree_info(this->tree->root);
}
void sample_propagation_forest::insert_node(solo_channel* node){
  // Fill in parent and children, and iterate over all trees of course.
  // Post-order traversal
  // Follow rules from paper to deduce first whether node can be a child of the current parent.
  assert(node != nullptr);
  bool is_comm = !(node->id.size()==1 && node->id[0].second==0);
  solo_channel* parent = nullptr;
  //TODO: I assume here that we care about the first SPT in the SPF. Figure out how to fix this later
  this->find_parent(this->tree->root,node,parent);
  node->parent = parent;
  assert(parent!=nullptr);
 
  // Try adding 'node' to each SPT. If none fit, append parent's children array and add it there, signifying a new tree, rooted at 'parent'
  bool valid_parent = false;
  int save_tree_idx=-1;
  for (auto i=0; i<parent->children.size(); i++){
    parent->children[i].push_back(node);
    std::vector<int> sibling_to_child_indices;
    for (auto j=0; j<parent->children[i].size()-1; j++){
      if (channel::verify_ancestor_relation(parent->children[i][j],node)){
        sibling_to_child_indices.push_back(j);
      }
    }
    bool sibling_decision = solo_channel::verify_sibling_relation(parent,i,sibling_to_child_indices);
    if (sibling_decision){
      save_tree_idx=i;
      valid_parent=true;
      for (auto j=0; j<sibling_to_child_indices.size(); j++){
        node->children[0].push_back(parent->children[i][sibling_to_child_indices[j]]);
        parent->children[i][sibling_to_child_indices[j]]->parent=node;
      }
      int skip_index=0;
      int save_index=0;
      for (auto j=0; j<parent->children[i].size()-1; j++){
        if ((skip_index<sibling_to_child_indices.size()) && (j==sibling_to_child_indices[skip_index])){
          skip_index++;
        } else{
          parent->children[i][save_index] = parent->children[i][j];
          save_index++;
        }
      }
      parent->children[i][save_index] = parent->children[i][parent->children[i].size()-1];
      save_index++;
      for (auto j=parent->children[i].size()-save_index; j>0; j--){
        parent->children[i].pop_back();
      }
      break;
    } else{
      parent->children[i].pop_back();
    }
  }
  if (!valid_parent){
    parent->children.push_back(std::vector<solo_channel*>());
    parent->children[parent->children.size()-1].push_back(node);
  }

  // sanity check for right now
  int world_comm_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_comm_rank);
  if (world_comm_rank==8){
    std::cout << "parent of node{ (" << node->local_hash_tag << "," << node->global_hash_tag << "," << node->offset;
    for (auto i=0; i<node->id.size(); i++){
      std::cout << ") (" << node->id[i].first << "," << node->id[i].second << ")";
    }
    std::cout << " } is { (" << parent->local_hash_tag << "," << parent->global_hash_tag << "," << parent->offset;
    for (auto i=0; i<parent->id.size(); i++){
      std::cout << ") (" << parent->id[i].first << "," << parent->id[i].second << ")";
    }
    std::cout << " }\n";
    for (auto i=0; i<parent->children.size(); i++){
      std::cout << "\tsubtree " << i << " contains " << parent->children[i].size() << " children\n";
      for (auto j=0; j<parent->children[i].size(); j++){
        std::cout << "\t\tchild " << j << " is { (" << parent->children[i][j]->local_hash_tag << "," << parent->children[i][j]->global_hash_tag << "," << parent->children[i][j]->offset;
        for (auto k=0; k<parent->children[i][j]->id.size(); k++){
          std::cout << ") (" << parent->children[i][j]->id[k].first << " " << parent->children[i][j]->id[k].second << ")";
        }
        std::cout << " }\n";
      }
    }
  }
}


// ****************************************************************************************************************************************************
void intermediate_stats::generate(const kernel& p, const std::vector<kernel_batch>& active_batches){
  this->M1 = p.M1;
  this->M2 = p.M2;
  this->num_schedules = p.num_schedules;
  this->num_local_schedules = p.num_local_schedules;
  this->num_scheduled_units = p.num_scheduled_units;
  this->num_local_scheduled_units = p.num_local_scheduled_units;
  this->total_exec_time = p.total_exec_time;
  this->total_local_exec_time = p.total_local_exec_time;
  for (auto i=0; i<active_batches.size(); i++){
    double M1_2 = active_batches[i].M1;
    double M2_2 = active_batches[i].M2;
    int n2 = active_batches[i].num_schedules;
    //assert(n2>0);
    double delta = this->M1 - M1_2;
    this->M1 = (this->num_schedules*this->M1 + n2*M1_2)/(this->num_schedules+n2);
    this->M2 = this->M2 + M2_2 + delta/(this->num_schedules+n2)*delta*(this->num_schedules*n2);
    this->num_schedules += n2;
    this->num_local_schedules += active_batches[i].num_local_schedules;
    this->num_scheduled_units += active_batches[i].num_scheduled_units;
    this->num_local_scheduled_units += active_batches[i].num_local_scheduled_units;
    this->total_exec_time += active_batches[i].total_exec_time;
    this->total_local_exec_time += active_batches[i].total_local_exec_time;
  }
}
intermediate_stats::intermediate_stats(const kernel_key_id& index, const std::vector<kernel_batch>& active_batches){
  this->generate(active_kernels[index.val_index],active_batches);
}

bool is_key_skipable(const comm_kernel_key& key){
  return true;
}
bool is_key_skipable(const comp_kernel_key& key){
  return true;
}

int get_skel_count(const comm_kernel_key& key){
  assert(sample_constraint_mode==3);
  int skel_val = 1;
  if (skeletonization::comm_kernel_map.find(key) != skeletonization::comm_kernel_map.end()){
    skel_val = skeletonization::active_kernels[skeletonization::comm_kernel_map[key].val_index];
    skel_val = std::max(skel_val,1);
  }
  return skel_val;
}
int get_skel_count(const comp_kernel_key& key){
  assert(sample_constraint_mode==3);
  int skel_val = 1;
  if (skeletonization::comp_kernel_map.find(key) != skeletonization::comp_kernel_map.end()){
    skel_val = skeletonization::active_kernels[skeletonization::comp_kernel_map[key].val_index];
    skel_val = std::max(skel_val,1);
  }
  return skel_val;
}


double get_estimate(const kernel& p, int analysis_param, double unit_count){
  if (analysis_param == 0){// arithmetic mean
    return get_arithmetic_mean(p);
  } else{
    return unit_count*get_harmonic_mean(p);
  }
}
double get_estimate(const kernel_propagate& p, int analysis_param, double unit_count){
  if (analysis_param == 0){// arithmetic mean
    return get_arithmetic_mean(p);
  } else{
    return unit_count*get_harmonic_mean(p);
  }
}
double get_estimate(const kernel_key_id& index, int analysis_param, double unit_count){
  if (analysis_param == 0){// arithmetic mean
    return get_arithmetic_mean(index);
  } else{
    return unit_count*get_harmonic_mean(index);
  }
}
double get_estimate(const kernel_key_id& index, const std::vector<kernel_batch>& active_batches, int analysis_param, double unit_count){
  auto stats = intermediate_stats(index,active_batches);
  if (analysis_param == 0){// arithmetic mean
    return stats.M1;
  } else{
    return unit_count*(1./stats.M1);
  }
}
double get_estimate(const intermediate_stats& p, int analysis_param, double unit_count){
  if (analysis_param == 0){// arithmetic mean
    return get_arithmetic_mean(p);
  } else{
    return unit_count*get_harmonic_mean(p);
  }
}

double get_arithmetic_mean(const kernel& p){
  // returns arithmetic mean
  return p.M1;
}
double get_arithmetic_mean(const kernel_propagate& p){
  // returns arithmetic mean
  return p.M1;
}
double get_arithmetic_mean(const kernel_key_id& index){
  // returns arithmetic mean
  return active_kernels[index.val_index].M1;
}
double get_arithmetic_mean(const intermediate_stats& p){
  // returns arithmetic mean
  return p.M1;
}

double get_harmonic_mean(const kernel& p){
  // returns arithmetic mean
  return 1./p.M1;
}
double get_harmonic_mean(const kernel_propagate& p){
  // returns arithmetic mean
  return 1./p.M1;
}
double get_harmonic_mean(const kernel_key_id& index){
  // returns arithmetic mean
  return 1./active_kernels[index.val_index].M1;
}
double get_harmonic_mean(const intermediate_stats& p){
  // returns arithmetic mean
  return 1./p.M1;
}

double get_variance(const kernel& p, int analysis_param){
  // returns variance
  size_t n = p.num_schedules;
  if (n<=1) return 1000000.;
  if (analysis_param == 0){
    return p.M2 / (n-1.);
  } else{
    return 1./p.M2 / (n-1.);
  }
}
double get_variance(const kernel_propagate& p, int analysis_param){
  // returns variance
  size_t n = p.num_schedules;
  if (n<=1) return 1000000.;
  if (analysis_param == 0){
    return p.M2 / (n-1.);
  } else{
    return 1./p.M2 / (n-1.);
  }
}
double get_variance(const kernel_key_id& index, int analysis_param){
  // returns variance
  size_t n = active_kernels[index.val_index].num_schedules;
  if (n<=1) return 1000000.;
  if (analysis_param == 0){
    return active_kernels[index.val_index].M2 / (n-1.);
  } else{
    return 1./active_kernels[index.val_index].M2 / (n-1.);
  }
}
double get_variance(const intermediate_stats& p, int analysis_param){
  // returns variance
  size_t n = p.num_schedules;
  if (n<=1) return 1000000.;
  if (analysis_param == 0){
    return p.M2 / (n-1.);
  } else{
    return 1./p.M2 / (n-1.);
  }
}

double get_std_dev(const kernel& p, int analysis_param){
  // returns variance
  return pow(get_variance(p,analysis_param),1./2.);
}
double get_std_dev(const kernel_propagate& p, int analysis_param){
  // returns variance
  return pow(get_variance(p,analysis_param),1./2.);
}
double get_std_dev(const kernel_key_id& index, int analysis_param){
  // returns variance
  return pow(get_variance(index,analysis_param),1./2.);
}
double get_std_dev(const intermediate_stats& p, int analysis_param){
  // returns variance
  return pow(get_variance(p,analysis_param),1./2.);
}

int get_std_error_count(const kernel& p){
  int n = 1;
  switch (sample_constraint_mode){
    case -1:
      n = 1; break;
    case 0:
      n=1; break;
    case 1:
      n = p.num_local_schedules; break;
    case 2:
      n = p.num_schedules; break;
  }
  n = std::max(1,n);
  return n;
}
int get_std_error_count(const kernel_propagate& p){
  int n = 1;
  switch (sample_constraint_mode){
    case -1:
      n = 1; break;
    case 0:
      n=1; break;
    case 1:
      n = p.num_local_schedules; break;
    case 2:
      n = p.num_schedules; break;
  }
  n = std::max(1,n);
  return n;
}
int get_std_error_count(const kernel_key_id& index){
  int n = 1;
  switch (sample_constraint_mode){
    case -1:
      n = 1; break;
    case 0:
      n=1; break;
    case 1:
      n = active_kernels[index.val_index].num_local_schedules; break;
    case 2:
      n = active_kernels[index.val_index].num_schedules; break;
  }
  n = std::max(1,n);
  return n;
}
int get_std_error_count(const intermediate_stats& p){
  // returns standard error
  int n = 1;
  switch (sample_constraint_mode){
    case -1:
      n = 1; break;
    case 0:
      n=1; break;
    case 1:
      n = p.num_local_schedules; break;
    case 2:
      n = p.num_schedules; break;
  }
  n = std::max(1,n);
  return n;
}

// Standard error calculation takes into account execution path analysis
double get_std_error(const comm_kernel_key& key, const kernel& p, int analysis_param){
  // returns standard error
  int n = get_std_error_count(p);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(p,analysis_param) / pow(n*1.,1./2.);
}
double get_std_error(const comp_kernel_key& key, const kernel& p, int analysis_param){
  // returns standard error
  int n = get_std_error_count(p);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(p,analysis_param) / pow(n*1.,1./2.);
}
double get_std_error(const comm_kernel_key& key, const kernel_propagate& p, int analysis_param){
  // returns standard error
  int n = get_std_error_count(p);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(p,analysis_param) / pow(n*1.,1./2.);
}
double get_std_error(const comp_kernel_key& key, const kernel_propagate& p, int analysis_param){
  // returns standard error
  int n = get_std_error_count(p);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(p,analysis_param) / pow(n*1.,1./2.);
}
double get_std_error(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param){
  // returns standard error
  int n = get_std_error_count(index);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(index,analysis_param) / pow(n*1.,1./2.);
}
double get_std_error(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param){
  // returns standard error
  int n = get_std_error_count(index);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(index,analysis_param) / pow(n*1.,1./2.);
}
double get_std_error(const comm_kernel_key& key, const intermediate_stats& p, int analysis_param){
  // returns standard error
  int n = get_std_error_count(p);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(p,analysis_param) / pow(n*1.,1./2.);
}
double get_std_error(const comp_kernel_key& key, const intermediate_stats& p, int analysis_param){
  // returns standard error
  int n = get_std_error_count(p);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(p,analysis_param) / pow(n*1.,1./2.);
}

double get_confidence_interval(const comm_kernel_key& key, const kernel& p, int analysis_param, double level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,p,analysis_param);
}
double get_confidence_interval(const comm_kernel_key& key, const kernel_propagate& p, int analysis_param, double level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,p,analysis_param);
}
double get_confidence_interval(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param, double level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,index,analysis_param);
}
double get_confidence_interval(const comm_kernel_key& key, const intermediate_stats& p, int analysis_param, double level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,p,analysis_param);
}
double get_confidence_interval(const comp_kernel_key& key, const kernel& p, int analysis_param, double level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,p,analysis_param);
}
double get_confidence_interval(const comp_kernel_key& key, const kernel_propagate& p, int analysis_param, double level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,p,analysis_param);
}
double get_confidence_interval(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param, double level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,index,analysis_param);
}
double get_confidence_interval(const comp_kernel_key& key, const intermediate_stats& p, int analysis_param, double level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,p,analysis_param);
}

bool is_steady(const comm_kernel_key& key, const kernel& p, int analysis_param){
  if (sample_constraint_mode==-1) { assert(p.num_schedules < 2); return true; }
  return ((get_confidence_interval(key,p,analysis_param) / get_estimate(p,analysis_param)) < kernel_error_limit) &&
          (p.num_schedules >= kernel_count_limit) &&
          (p.total_exec_time >= kernel_time_limit);
}
bool is_steady(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param){
  if (sample_constraint_mode==-1) { assert(active_kernels[index.val_index].num_schedules < 2); return true; }
  return ((get_confidence_interval(key,index,analysis_param) / get_estimate(index,analysis_param)) < kernel_error_limit) &&
          (active_kernels[index.val_index].num_schedules >= kernel_count_limit) &&
          (active_kernels[index.val_index].total_exec_time >= kernel_time_limit);
}
bool is_steady(const comm_kernel_key& key, const intermediate_stats& p, int analysis_param){
  if (sample_constraint_mode==-1) { assert(p.num_schedules < 2); return true; }
  return ((get_confidence_interval(key,p,analysis_param) / get_estimate(p,analysis_param)) < kernel_error_limit) &&
          (p.num_schedules >= kernel_count_limit) &&
          (p.total_exec_time >= kernel_time_limit);
}
bool is_steady(const comp_kernel_key& key, const kernel& p, int analysis_param){
  if (sample_constraint_mode==-1) { assert(p.num_schedules < 2); return true; }
  return ((get_confidence_interval(key,p,analysis_param) / get_estimate(p,analysis_param)) < kernel_error_limit) &&
          (p.num_schedules >= kernel_count_limit) &&
          (p.total_exec_time >= kernel_time_limit);
}
bool is_steady(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param){
  if (sample_constraint_mode==-1) { assert(active_kernels[index.val_index].num_schedules < 2); return true; }
  return ((get_confidence_interval(key,index,analysis_param) / get_estimate(index,analysis_param)) < kernel_error_limit) &&
          (active_kernels[index.val_index].num_schedules >= kernel_count_limit) &&
          (active_kernels[index.val_index].total_exec_time >= kernel_time_limit);
}
bool is_steady(const comp_kernel_key& key, const intermediate_stats& p, int analysis_param){
  if (sample_constraint_mode==-1) { assert(p.num_schedules < 2); return true; }
  return ((get_confidence_interval(key,p,analysis_param) / get_estimate(p,analysis_param)) < kernel_error_limit) &&
          (p.num_schedules >= kernel_count_limit) &&
          (p.total_exec_time >= kernel_time_limit);
}

double get_error_estimate(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param){
  auto& active_batches = comm_batch_map[key];
  auto stats = intermediate_stats(index,active_batches);
  return get_confidence_interval(key,stats,analysis_param) / (get_estimate(stats,analysis_param));
}
double get_error_estimate(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param){
  auto& active_batches = comp_batch_map[key];
  auto stats = intermediate_stats(index,active_batches);
  return get_confidence_interval(key,stats,analysis_param) / (get_estimate(stats,analysis_param));
}
double get_error_estimate(const comm_kernel_key& key, const kernel_propagate& p, int analysis_param){
  return get_confidence_interval(key,p,analysis_param) / (get_estimate(p,analysis_param));
}
double get_error_estimate(const comp_kernel_key& key, const kernel_propagate& p, int analysis_param){
  return get_confidence_interval(key,p,analysis_param) / (get_estimate(p,analysis_param));
}

bool steady_test(const comm_kernel_key& key, const kernel& p, int analysis_param){
  if (!is_key_skipable(key)) return false;
  return is_steady(key,p,analysis_param);
}
bool steady_test(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param){
  if (!is_key_skipable(key)) return false;
  if (comm_sample_aggregation_mode == 0 || comm_batch_map.find(key) == comm_batch_map.end()){
    return is_steady(key,index,analysis_param);
  }
  else if (comm_sample_aggregation_mode == 1){
    //TODO: We may not want to allow this
    auto& active_batches = comm_batch_map[key];
    auto stats = intermediate_stats(index,active_batches);
    return is_steady(key,stats,analysis_param);
  }
}
bool steady_test(const comp_kernel_key& key, const kernel& p, int analysis_param){
  if (!is_key_skipable(key)) return false;
  return is_steady(key,p,analysis_param);
}
bool steady_test(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param){
  if (!is_key_skipable(key)) return false;
  if (comp_sample_aggregation_mode == 0 || comp_batch_map.find(key) == comp_batch_map.end()){
    return is_steady(key,index,analysis_param);
  }
  else if (comp_sample_aggregation_mode == 1){
    //TODO: We may not want to allow this
    auto& active_batches = comp_batch_map[key];
    auto stats = intermediate_stats(index,active_batches);
    return is_steady(key,stats,analysis_param);
  }
}

void update_kernel_stats(kernel& p, int analysis_param, volatile double exec_time, double unit_count){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  if (exec_time == 0) { exec_time=1.e-9; }
  if (p.steady_state == 0){
    p.num_schedules++;
    p.num_local_schedules++;
    p.num_scheduled_units += unit_count;
    p.num_local_scheduled_units += unit_count;
    p.total_exec_time += exec_time;
    p.total_local_exec_time += exec_time;
    // Online computation of up to 4th-order central moments using compunication time samples
    size_t n1 = p.num_schedules-1;
    size_t n = p.num_schedules;
    double x;
    if (analysis_param == 0){x = exec_time; }	// prep for arithmetic mean
    else                   {x = (unit_count>0 ? unit_count : 1.)/exec_time; }	// prep for harmonic mean
    double delta = x - p.M1;
    double delta_n = delta / n;
    double delta_n2 = delta_n*delta_n;
    double term1 = delta*delta_n*n1;
    p.M1 += delta_n;
    p.M2 += term1;
  }
  else{
    p.num_non_schedules++;
    p.num_non_scheduled_units += unit_count;
  }
}
void update_kernel_stats(const kernel_key_id& index, int analysis_param, volatile double exec_time, double unit_count){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  if (exec_time == 0) { exec_time=1.e-9; }
  if (active_kernels[index.val_index].steady_state == 0){
    active_kernels[index.val_index].num_schedules++;
    active_kernels[index.val_index].num_local_schedules++;
    active_kernels[index.val_index].num_scheduled_units += unit_count;
    active_kernels[index.val_index].num_local_scheduled_units += unit_count;
    active_kernels[index.val_index].total_exec_time += exec_time;
    active_kernels[index.val_index].total_local_exec_time += exec_time;
    // Online computation of up to 4th-order central moments using compunication time samples
    size_t n1 = active_kernels[index.val_index].num_schedules-1;
    size_t n = active_kernels[index.val_index].num_schedules;
    double x;
    if (analysis_param == 0){x = exec_time; }	// prep for arithmetic mean
    else                   {x = (unit_count>0 ? unit_count : 1.)/exec_time; }	// prep for harmonic mean
    double delta = x - active_kernels[index.val_index].M1;
    double delta_n = delta / n;
    double delta_n2 = delta_n*delta_n;
    double term1 = delta*delta_n*n1;
    active_kernels[index.val_index].M1 += delta_n;
    active_kernels[index.val_index].M2 += term1;
  }
  else{
    active_kernels[index.val_index].num_non_schedules++;
    active_kernels[index.val_index].num_non_scheduled_units += unit_count;
  }
}
void update_kernel_stats(kernel& dest, const kernel& src, int analysis_param){
  // This function will implement the parallel algorithm computing the mean and variance of two partitions
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  // Online computation of up to 4th-order central moments using compunication time samples
  size_t n1 = dest.num_schedules;
  size_t n2 = src.num_schedules;
  double delta = dest.M1 - src.M1;
  dest.M1 = (n1*dest.M1 + n2*src.M1)/(n1+n2);
  dest.M2 = dest.M2 + src.M2 + delta/(n1+n2)*delta*(n1*n2);
  dest.num_schedules += src.num_schedules;
  dest.num_local_schedules += src.num_local_schedules;
  dest.num_scheduled_units += src.num_scheduled_units;
  dest.num_local_scheduled_units += src.num_local_scheduled_units;
  dest.num_non_schedules += src.num_non_schedules;
  dest.num_non_scheduled_units += src.num_non_scheduled_units;// Wouldn't these be zero?
  dest.total_exec_time += src.total_exec_time;
  dest.total_local_exec_time += src.total_local_exec_time;
}

void update_kernel_stats(kernel_batch& batch, int analysis_param, volatile double exec_time, double unit_count){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  if (exec_time == 0) { exec_time=1.e-9; }
  batch.num_schedules++;
  batch.num_local_schedules++;
  batch.num_scheduled_units += unit_count;
  batch.num_local_scheduled_units += unit_count;
  batch.total_exec_time += exec_time;
  batch.total_local_exec_time += exec_time;
  // Online computation of up to 4th-order central moments using compunication time samples
  size_t n1 = batch.num_schedules-1;
  size_t n = batch.num_schedules;
  double x;
  if (analysis_param == 0){x = exec_time; }	// prep for arithmetic mean
  else                   {x = (unit_count>0 ? unit_count : 1.)/exec_time; }	// prep for harmonic mean
  double delta = x - batch.M1;
  double delta_n = delta / n;
  double delta_n2 = delta_n*delta_n;
  double term1 = delta*delta_n*n1;
  batch.M1 += delta_n;
  batch.M2 += term1;
}
void update_kernel_stats(kernel& dest, const kernel_batch& src, int analysis_param){
  // This function will implement the parallel algorithm computing the mean and variance of two partitions
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  // Online computation of up to 4th-order central moments using compunication time samples
  size_t n1 = dest.num_schedules;
  size_t n2 = src.num_schedules;
  double delta = dest.M1 - src.M1;
  dest.M1 = (n1*dest.M1 + n2*src.M1)/(n1+n2);
  dest.M2 = dest.M2 + src.M2 + delta/(n1+n2)*delta*(n1*n2);
  dest.num_schedules += src.num_schedules;
  dest.num_local_schedules += src.num_local_schedules;
  dest.num_scheduled_units += src.num_scheduled_units;
  dest.num_local_scheduled_units += src.num_local_scheduled_units;
  dest.total_exec_time += src.total_exec_time;
  dest.total_local_exec_time += src.total_local_exec_time;
}
void update_kernel_stats(kernel_batch& dest, const kernel_batch& src, int analysis_param){
  // This function will implement the parallel algorithm computing the mean and variance of two partitions
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  // Online computation of up to 4th-order central moments using compunication time samples
  size_t n1 = dest.num_schedules;
  size_t n2 = src.num_schedules;
  double delta = dest.M1 - src.M1;
  dest.M1 = (n1*dest.M1 + n2*src.M1)/(n1+n2);
  dest.M2 = dest.M2 + src.M2 + delta/(n1+n2)*delta*(n1*n2);
  dest.num_schedules += src.num_schedules;
  dest.num_scheduled_units += src.num_scheduled_units;
  dest.total_exec_time += src.total_exec_time;
 // Don't update the three 'local' members
}
void update_kernel_stats(const kernel_key_id& index, const intermediate_stats& stats){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  active_kernels[index.val_index].num_schedules = stats.num_schedules;
  active_kernels[index.val_index].num_local_schedules = stats.num_local_schedules;
  active_kernels[index.val_index].num_scheduled_units = stats.num_scheduled_units;
  active_kernels[index.val_index].num_local_scheduled_units = stats.num_local_scheduled_units;
  active_kernels[index.val_index].total_exec_time = stats.total_exec_time;
  active_kernels[index.val_index].total_local_exec_time = stats.total_local_exec_time;
  active_kernels[index.val_index].M1 = stats.M1;
  active_kernels[index.val_index].M2 = stats.M2;
}

int should_schedule(const kernel& p){
  return (p.global_steady_state==1 ? 0 : 1);
}
int should_schedule(const kernel_key_id& index){
  return (active_kernels[index.val_index].global_steady_state==1 ? 0 : 1);
}

void set_kernel_state(kernel& p, bool schedule_decision){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  if (p.steady_state == 1) return;// Don't allow a kernel to go from steady back to unsteady
  p.steady_state = (schedule_decision==true ? 0 : 1);
}
void set_kernel_state(const kernel_key_id& index, bool schedule_decision){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  if (active_kernels[index.val_index].steady_state == 1) return;// Don't allow a kernel to go from steady back to unsteady
  active_kernels[index.val_index].steady_state = (schedule_decision==true ? 0 : 1);
}

void set_kernel_state_global(kernel& p, bool schedule_decision){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  if (p.global_steady_state == 1) return;// Don't allow a kernel to go from steady back to unsteady
  p.global_steady_state = (schedule_decision==true ? 0 : 1);
}
void set_kernel_state_global(const kernel_key_id& index, bool schedule_decision){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  if (active_kernels[index.val_index].global_steady_state == 1) return;// Don't allow a kernel to go from steady back to unsteady
  active_kernels[index.val_index].global_steady_state = (schedule_decision==true ? 0 : 1);
}

void merge_batches(std::vector<kernel_batch>& batches, int analysis_param){
  if (batches.size() == 0) return;
  // At first, I wanted to leverage the 'is_final' info from the aggregate (which can be obtained from a batch's hash_id
  //   However, I think I don't need to do this. I can just keep that state within the batch and accumulate in there.
  // Iterate over the entires to merge any two batches with the same state.
  std::sort(batches.begin(),batches.end(),[](const kernel_batch& p1, const kernel_batch& p2){return p1.hash_id < p2.hash_id;});
   // Assumption: no more than 2 batches can have the same state.
  int start_index = 0;
  for (auto i=1; i<batches.size(); i++){
    if (batches[i].hash_id == batches[start_index].hash_id){
      // Merge batches[i] into batches[start_index]
      update_kernel_stats(batches[start_index],batches[i],analysis_param);
      batches[start_index].num_local_schedules += batches[i].num_local_schedules;
      batches[start_index].num_local_scheduled_units += batches[i].num_local_scheduled_units;
      batches[start_index].total_local_exec_time += batches[i].total_local_exec_time;
    } else{
      batches[++start_index] = batches[i];
    }
  }
  for (int i=batches.size()-start_index-1; i>0; i--){
    batches.pop_back();
  }
}


void allocate(MPI_Comm comm){
  int _world_size; MPI_Comm_size(MPI_COMM_WORLD,&_world_size);
  int _world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&_world_rank);
  mode_1_width = 25;
  mode_2_width = 15;
  communicator_count=INT_MIN;// to avoid conflict with p2p, which could range from (-p,p)
  internal_tag = 31133;
  internal_tag1 = internal_tag+1;
  internal_tag2 = internal_tag+2;
  internal_tag3 = internal_tag+3;
  internal_tag4 = internal_tag+4;
  internal_tag5 = internal_tag+5;
  is_first_iter = true;
  intercept_overhead.resize(3,0);
  global_intercept_overhead.resize(3,0);
  global_comm_kernel_stats.resize(5,0);
  global_comp_kernel_stats.resize(5,0);
  local_comm_kernel_stats.resize(5,0);
  local_comp_kernel_stats.resize(5,0);
  save_comp_kernel_stats.resize(2,0);
  save_comm_kernel_stats.resize(2,0);

  solo_channel* world_node = new solo_channel();
  world_node->tag = communicator_count++;
  world_node->offset = 0;
  world_node->id.push_back(std::make_pair(_world_size,1));
  world_node->parent=nullptr;
  std::string local_channel_hash_str = ".." + std::to_string(world_node->id[0].first) + "." + std::to_string(world_node->id[0].second);
  std::string global_channel_hash_str = ".." + std::to_string(world_node->id[0].first) + "." + std::to_string(world_node->id[0].second);
  world_node->local_hash_tag = std::hash<std::string>()(local_channel_hash_str);// will avoid any local overlap.
  world_node->global_hash_tag = std::hash<std::string>()(global_channel_hash_str);// will avoid any global overlap.
  sample_propagation_tree* tree = new sample_propagation_tree;
  tree->root = world_node;
  spf.tree = tree;
  comm_channel_map[MPI_COMM_WORLD] = world_node;
  // Always treat 1-communicator channels as trivial aggregate channels.
  aggregate_channel* agg_node = new aggregate_channel(world_node->id,world_node->local_hash_tag,world_node->global_hash_tag,world_node->offset,1);
  agg_node->channels.insert(world_node->local_hash_tag);
  assert(aggregate_channel_map.find(world_node->local_hash_tag) == aggregate_channel_map.end());
  aggregate_channel_map[world_node->local_hash_tag] = agg_node;
  aggregate_channel_map[world_node->local_hash_tag]->is_final=true;
  if (_world_rank == 8){
    auto str1 = channel::generate_tuple_string(agg_node);
    auto str2 = aggregate_channel::generate_hash_history(agg_node);
    std::cout << "Process " << _world_rank << " has aggregate " << str1 << " " << str2 << " with hashes (" << agg_node->local_hash_tag << " " << agg_node->global_hash_tag << "), num_channels - " << agg_node->num_channels << std::endl;
  }

  kernel_propagate ex_3;
  MPI_Datatype kernel_internal_type[2] = { MPI_INT, MPI_DOUBLE };
  int kernel_internal_block_len[2] = { 3,6 };
  MPI_Aint kernel_internal_disp[2] = { (char*)&ex_3.hash_id-(char*)&ex_3, (char*)&ex_3.num_scheduled_units-(char*)&ex_3 };
  PMPI_Type_create_struct(2,kernel_internal_block_len,kernel_internal_disp,kernel_internal_type,&kernel_type);
  PMPI_Type_commit(&kernel_type);

  kernel_batch_propagate ex_4;
  MPI_Datatype batch_internal_type[2] = { MPI_INT, MPI_DOUBLE };
  int batch_internal_block_len[2] = { 3,4 };
  MPI_Aint batch_internal_disp[2] = { (char*)&ex_4.hash_id-(char*)&ex_4, (char*)&ex_4.num_scheduled_units-(char*)&ex_4 };
  PMPI_Type_create_struct(2,batch_internal_block_len,batch_internal_disp,batch_internal_type,&batch_type);
  PMPI_Type_commit(&batch_type);

  //TODO: Not a fan of these magic numbers '2' and '9'. Should utilize some error checking for strings that are not of proper length anyways.

  // Communication kernel time, computation kernel time, computation time, execution time
  num_critical_path_measures 		= 4;
  num_per_process_measures 		= 4;
  num_volume_measures 			= 4;

  critical_path_costs_size            	= num_critical_path_measures;
  per_process_costs_size              	= num_per_process_measures;
  volume_costs_size                   	= num_volume_measures;

  critical_path_costs.resize(critical_path_costs_size);
  max_per_process_costs.resize(per_process_costs_size);
  volume_costs.resize(volume_costs_size);
  new_cs.resize(critical_path_costs_size);
  info_sender.resize(num_critical_path_measures);
  info_receiver.resize(num_critical_path_measures);

  // Reset these global variables, as some are updated by function arguments for convenience
  if (std::getenv("CRITTER_VIZ_FILE") != NULL){
    std::string stream_name = std::getenv("CRITTER_VIZ_FILE");
    std::string stream_tune_name = std::getenv("CRITTER_VIZ_FILE");
    std::string stream_reconstruct_name = std::getenv("CRITTER_VIZ_FILE");
    stream_name += "_discretization.txt";
    stream_tune_name += "_discretization_tune.txt";
    stream_reconstruct_name += "_discretization_reconstruct.txt";
    if (is_world_root){
      stream.open(stream_name.c_str(),std::ofstream::app);
      stream_tune.open(stream_tune_name.c_str(),std::ofstream::app);
      stream_reconstruct.open(stream_reconstruct_name.c_str(),std::ofstream::app);
    }
  }
  if (std::getenv("CRITTER_AUTOTUNING_DELTA") != NULL){
    tuning_delta = atoi(std::getenv("CRITTER_AUTOTUNING_DELTA"));
    assert(tuning_delta>=0 && tuning_delta<=3);// tuning_delta==0 requires decomposition mechanism
  } else{
    tuning_delta = 0;
  }
  if (std::getenv("CRITTER_AUTOTUNING_COMM_SAMPLE_AGGREGATION_MODE") != NULL){
    comm_sample_aggregation_mode = atoi(std::getenv("CRITTER_AUTOTUNING_COMM_SAMPLE_AGGREGATION_MODE"));
    assert(comm_sample_aggregation_mode>=0 && comm_sample_aggregation_mode<=1);
  } else{
    comm_sample_aggregation_mode = 0;
  }
  if (std::getenv("CRITTER_AUTOTUNING_COMM_STATE_AGGREGATION_MODE") != NULL){
    comm_state_aggregation_mode = atoi(std::getenv("CRITTER_AUTOTUNING_COMM_STATE_AGGREGATION_MODE"));
    assert(comm_state_aggregation_mode>=0 && comm_state_aggregation_mode<=2);
  } else{
    comm_state_aggregation_mode = 0;
  }
  if (std::getenv("CRITTER_AUTOTUNING_COMP_SAMPLE_AGGREGATION_MODE") != NULL){
    comp_sample_aggregation_mode = atoi(std::getenv("CRITTER_AUTOTUNING_COMP_SAMPLE_AGGREGATION_MODE"));
    assert(comp_sample_aggregation_mode>=0 && comp_sample_aggregation_mode<=1);
  } else{
    comp_sample_aggregation_mode = 0;
  }
  if (std::getenv("CRITTER_AUTOTUNING_COMP_STATE_AGGREGATION_MODE") != NULL){
    comp_state_aggregation_mode = atoi(std::getenv("CRITTER_AUTOTUNING_COMP_STATE_AGGREGATION_MODE"));
    assert(comp_state_aggregation_mode>=0 && comp_state_aggregation_mode<=2);
  } else{
    comp_state_aggregation_mode = 0;
  }
  assert(comm_state_aggregation_mode >= comm_sample_aggregation_mode);
  assert(comp_state_aggregation_mode >= comp_sample_aggregation_mode);
  if (std::getenv("CRITTER_AUTOTUNING_SAMPLE_CONSTRAINT_MODE") != NULL){
    sample_constraint_mode = atoi(std::getenv("CRITTER_AUTOTUNING_SAMPLE_CONSTRAINT_MODE"));
    assert(sample_constraint_mode>=-1 && sample_constraint_mode<=3);
  } else{
    sample_constraint_mode = 0;
  }
/*
  if (std::getenv("CRITTER_AUTOTUNING_SAMPLE_RESET_MODE") != NULL){
    sample_reset_mode = atoi(std::getenv("CRITTER_AUTOTUNING_SAMPLE_RESET_MODE"));
    assert(sample_reset_mode>=0 && sample_reset_mode<=3);
  } else{
    sample_reset_mode = 0;
  }
*/
  if (std::getenv("CRITTER_COMM_ENVELOPE_PARAM") != NULL){
    comm_envelope_param = atoi(std::getenv("CRITTER_COMM_ENVELOPE_PARAM"));
  } else{
    comm_envelope_param = 0;
  }
  if (std::getenv("CRITTER_COMM_UNIT_PARAM") != NULL){
    comm_unit_param = atoi(std::getenv("CRITTER_COMM_UNIT_PARAM"));
  } else{
    comm_unit_param = 0;
  }
  if (std::getenv("CRITTER_COMM_ANALYSIS_PARAM") != NULL){
    comm_analysis_param = atoi(std::getenv("CRITTER_COMM_ANALYSIS_PARAM"));
  } else{
    comm_analysis_param = 0;
  }
  if (std::getenv("CRITTER_COMP_ENVELOPE_PARAM") != NULL){
    comp_envelope_param = atoi(std::getenv("CRITTER_COMP_ENVELOPE_PARAM"));
  } else{
    comp_envelope_param = 0;
  }
  if (std::getenv("CRITTER_COMP_UNIT_PARAM") != NULL){
    comp_unit_param = atoi(std::getenv("CRITTER_COMP_UNIT_PARAM"));
  } else{
    comp_unit_param = 0;
  }
  if (std::getenv("CRITTER_COMP_ANALYSIS_PARAM") != NULL){
    comp_analysis_param = atoi(std::getenv("CRITTER_COMP_ANALYSIS_PARAM"));
  } else{
    comp_analysis_param = 0;
  }
  if (std::getenv("CRITTER_PATTERN_COUNT_LIMIT") != NULL){
    kernel_count_limit = atoi(std::getenv("CRITTER_PATTERN_COUNT_LIMIT"));
  } else{
    kernel_count_limit = 2;
  }
  if (std::getenv("CRITTER_PATTERN_TIME_LIMIT") != NULL){
    kernel_time_limit = atof(std::getenv("CRITTER_PATTERN_TIME_LIMIT"));
  } else{
    kernel_time_limit = 0;
  }
  if (std::getenv("CRITTER_PATTERN_ERROR_LIMIT") != NULL){
    kernel_error_limit = atof(std::getenv("CRITTER_PATTERN_ERROR_LIMIT"));
  } else{
    kernel_error_limit = .5;
  }
  if (std::getenv("CRITTER_COMP_KERNEL_TRANSFER_ID") != NULL){
    comp_kernel_transfer_id = atof(std::getenv("CRITTER_COMP_KERNEL_TRANSFER_ID"));
  } else{
    comp_kernel_transfer_id = 0;
  }
  if (std::getenv("CRITTER_COMM_KERNEL_TRANSFER_ID") != NULL){
    comm_kernel_transfer_id = atof(std::getenv("CRITTER_COMM_KERNEL_TRANSFER_ID"));
  } else{
    comm_kernel_transfer_id = 0;
  }
  if (std::getenv("CRITTER_COMP_KERNEL_BUFFER_ID") != NULL){
    comp_kernel_buffer_id = atof(std::getenv("CRITTER_COMP_KERNEL_BUFFER_ID"));
  } else{
    comp_kernel_buffer_id = 0;
  }


  if (eager_p2p){
    int eager_msg_sizes[5];
    MPI_Pack_size(1,MPI_CHAR,comm,&eager_msg_sizes[0]);
    MPI_Pack_size(1,MPI_CHAR,comm,&eager_msg_sizes[1]);
    MPI_Pack_size(num_critical_path_measures,MPI_DOUBLE_INT,comm,&eager_msg_sizes[2]);
    MPI_Pack_size(critical_path_costs_size,MPI_DOUBLE,comm,&eager_msg_sizes[3]);
    MPI_Pack_size(1,MPI_INT,comm,&eager_msg_sizes[4]);
    int eager_pad_size = 5*MPI_BSEND_OVERHEAD;
    for (int i=0; i<8; i++) { eager_pad_size += eager_msg_sizes[i]; }
    eager_pad.resize(eager_pad_size);
  }
}

void open_symbol(const char* symbol, double curtime){}

void close_symbol(const char* symbol, double curtime){}

void final_accumulate(MPI_Comm comm, double last_time){
  assert(internal_comm_info.size() == 0);
  critical_path_costs[num_critical_path_measures-1]+=(last_time-computation_timer);	// update critical path runtime
  critical_path_costs[num_critical_path_measures-3]+=(last_time-computation_timer);	// update critical path computation time
  volume_costs[num_volume_measures-1]+=(last_time-computation_timer);			// update per-process execution time
  volume_costs[num_volume_measures-3]+=(last_time-computation_timer);			// update per-process execution time

  _wall_time = wall_timer[wall_timer.size()-1];

  double temp_costs[4+4+3+5+5+2+3+3];
  for (int i=0; i<18; i++){ temp_costs[11+i]=0; }

  discretization::_MPI_Barrier.comm = MPI_COMM_WORLD;
  discretization::_MPI_Barrier.save_comp_key.clear();
  discretization::_MPI_Barrier.save_comm_key.clear();
  for (auto& it : comp_kernel_map){
    if ((active_kernels[it.second.val_index].steady_state==1) && (should_schedule(it.second)==1)){
      discretization::_MPI_Barrier.save_comp_key.push_back(it.first);
      temp_costs[critical_path_costs.size()+max_per_process_costs.size()+13]++;
    }

    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+3] += active_kernels[it.second.val_index].num_local_schedules;
    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+4] += active_kernels[it.second.val_index].num_non_schedules;
    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+5] += active_kernels[it.second.val_index].num_local_scheduled_units;
    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+6] += active_kernels[it.second.val_index].num_non_scheduled_units;
    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+7] += active_kernels[it.second.val_index].total_local_exec_time;
    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+15] += active_kernels[it.second.val_index].num_local_schedules;
    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+16] += active_kernels[it.second.val_index].num_local_scheduled_units;
    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+17] += active_kernels[it.second.val_index].total_local_exec_time;
    // Invalid assert-> this case is possible assert(comp_batch_map.find(it.first) != comp_batch_map.end());
    // if (comp_batch_map[it.first].size() == 0) continue;
    // Complete any active and incomplete batches. Liquidate them into the pathset.
    if (comp_sample_aggregation_mode == 1){
      auto stats = intermediate_stats(it.second,comp_batch_map[it.first]);
      update_kernel_stats(comp_kernel_map[it.first],stats);
      comp_batch_map[it.first].clear();
    }
  }
  for (auto& it : comm_kernel_map){
    if ((active_kernels[it.second.val_index].steady_state==1) && (should_schedule(it.second)==1)){
      discretization::_MPI_Barrier.save_comm_key.push_back(it.first);
      temp_costs[critical_path_costs.size()+max_per_process_costs.size()+14]++;
    }

    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+8] += active_kernels[it.second.val_index].num_local_schedules;
    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+9] += active_kernels[it.second.val_index].num_non_schedules;
    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+10] += active_kernels[it.second.val_index].num_local_scheduled_units;
    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+11] += active_kernels[it.second.val_index].num_non_scheduled_units;
    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+12] += active_kernels[it.second.val_index].total_local_exec_time;
    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+18] += active_kernels[it.second.val_index].num_local_schedules;
    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+19] += active_kernels[it.second.val_index].num_local_scheduled_units;
    temp_costs[critical_path_costs.size()+max_per_process_costs.size()+20] += active_kernels[it.second.val_index].total_local_exec_time;
    // Invalid assert-> this case is possible assert(comm_batch_map.find(it.first) != comm_batch_map.end());
    //if (comm_batch_map[it.first].size() == 0) continue;
    // Complete any active and incomplete batches. Liquidate them into the pathset.
    if (comm_sample_aggregation_mode == 1){
      auto stats = intermediate_stats(it.second,comm_batch_map[it.first]);
      update_kernel_stats(comm_kernel_map[it.first],stats);
      comm_batch_map[it.first].clear();
    }
  }

  max_per_process_costs = volume_costs;// copy over the per-process measurements that exist in volume_costs
  for (auto i=0; i<critical_path_costs.size(); i++) temp_costs[i] = critical_path_costs[i];
  for (auto i=0; i<max_per_process_costs.size(); i++) temp_costs[critical_path_costs.size()+i] = max_per_process_costs[i];
  temp_costs[critical_path_costs.size()+max_per_process_costs.size()] = intercept_overhead[0];
  temp_costs[critical_path_costs.size()+max_per_process_costs.size()+1] = intercept_overhead[1];
  temp_costs[critical_path_costs.size()+max_per_process_costs.size()+2] = intercept_overhead[2];
  PMPI_Allreduce(MPI_IN_PLACE,&temp_costs[0],4+4+3+5+5+2+3+3,MPI_DOUBLE,MPI_MAX,comm);
  for (auto i=0; i<critical_path_costs.size(); i++) critical_path_costs[i] = temp_costs[i];
  for (auto i=0; i<max_per_process_costs.size(); i++) max_per_process_costs[i] = temp_costs[critical_path_costs.size()+i];
  global_intercept_overhead[0] += temp_costs[critical_path_costs.size()+max_per_process_costs.size()];
  global_intercept_overhead[1] += temp_costs[critical_path_costs.size()+max_per_process_costs.size()+1];
  global_intercept_overhead[2] += temp_costs[critical_path_costs.size()+max_per_process_costs.size()+2];
  global_comp_kernel_stats[0] += temp_costs[critical_path_costs.size()+max_per_process_costs.size()+3];
  global_comp_kernel_stats[1] = temp_costs[critical_path_costs.size()+max_per_process_costs.size()+4];
  global_comp_kernel_stats[2] += temp_costs[critical_path_costs.size()+max_per_process_costs.size()+5];
  global_comp_kernel_stats[3] = temp_costs[critical_path_costs.size()+max_per_process_costs.size()+6];
  global_comp_kernel_stats[4] += temp_costs[critical_path_costs.size()+max_per_process_costs.size()+7];
  global_comm_kernel_stats[0] += temp_costs[critical_path_costs.size()+max_per_process_costs.size()+8];
  global_comm_kernel_stats[1] = temp_costs[critical_path_costs.size()+max_per_process_costs.size()+9];
  global_comm_kernel_stats[2] += temp_costs[critical_path_costs.size()+max_per_process_costs.size()+10];
  global_comm_kernel_stats[3] = temp_costs[critical_path_costs.size()+max_per_process_costs.size()+11];
  global_comm_kernel_stats[4] += temp_costs[critical_path_costs.size()+max_per_process_costs.size()+12];
  local_comp_kernel_stats[0] = temp_costs[critical_path_costs.size()+max_per_process_costs.size()+15];
  local_comp_kernel_stats[1] = global_comp_kernel_stats[1] - save_comp_kernel_stats[0];
  local_comp_kernel_stats[2] = temp_costs[critical_path_costs.size()+max_per_process_costs.size()+16];
  local_comp_kernel_stats[3] = global_comp_kernel_stats[3] - save_comp_kernel_stats[1];
  local_comp_kernel_stats[4] = temp_costs[critical_path_costs.size()+max_per_process_costs.size()+17];
  local_comm_kernel_stats[0] = temp_costs[critical_path_costs.size()+max_per_process_costs.size()+18];
  local_comm_kernel_stats[1] = global_comm_kernel_stats[1] - save_comm_kernel_stats[0];
  local_comm_kernel_stats[2] = temp_costs[critical_path_costs.size()+max_per_process_costs.size()+19];
  local_comm_kernel_stats[3] = global_comm_kernel_stats[3] - save_comm_kernel_stats[1];
  local_comm_kernel_stats[4] = temp_costs[critical_path_costs.size()+max_per_process_costs.size()+20];
  save_comp_kernel_stats[0] = global_comp_kernel_stats[1];
  save_comp_kernel_stats[1] = global_comp_kernel_stats[3];
  save_comm_kernel_stats[0] = global_comm_kernel_stats[1];
  save_comm_kernel_stats[1] = global_comm_kernel_stats[3];
  PMPI_Allreduce(MPI_IN_PLACE,&volume_costs[0],volume_costs.size(),MPI_DOUBLE,MPI_SUM,comm);
  discretization::_MPI_Barrier.aggregate_comp_kernels = temp_costs[critical_path_costs.size()+max_per_process_costs.size()+13]>0;
  discretization::_MPI_Barrier.aggregate_comm_kernels = temp_costs[critical_path_costs.size()+max_per_process_costs.size()+24]>0;
  discretization::_MPI_Barrier.should_propagate = discretization::_MPI_Barrier.aggregate_comp_kernels>0 || discretization::_MPI_Barrier.aggregate_comm_kernels>0;
}

void reset(bool schedule_kernels_override, bool force_steady_statistical_data_overide){
  assert(internal_comm_info.size() == 0);
  memset(&critical_path_costs[0],0,sizeof(double)*critical_path_costs.size());
  memset(&max_per_process_costs[0],0,sizeof(double)*max_per_process_costs.size());
  memset(&volume_costs[0],0,sizeof(double)*volume_costs.size());
  internal::bsp_counter=0;
  memset(&intercept_overhead[0],0,sizeof(double)*intercept_overhead.size());
  memset(&local_comp_kernel_stats[0],0,sizeof(double)*local_comp_kernel_stats.size());
  memset(&local_comm_kernel_stats[0],0,sizeof(double)*local_comm_kernel_stats.size());
  skip_ptr_set.clear();

  // This reset will no longer reset the kernel state, but will reset the schedule counters
  for (auto& it : comp_kernel_map){
    active_kernels[it.second.val_index].reset();
  }
  for (auto& it : comm_kernel_map){
    active_kernels[it.second.val_index].reset();
  }

  if (std::getenv("CRITTER_MODE") != NULL){
    internal::mode = atoi(std::getenv("CRITTER_MODE"));
  } else{
    internal::mode = 1;
  }
  if (std::getenv("CRITTER_SCHEDULE_KERNELS") != NULL){
    schedule_kernels = atoi(std::getenv("CRITTER_SCHEDULE_KERNELS"));
  } else{
    schedule_kernels = 1;
  }
  if (schedule_kernels==1){ schedule_kernels = (schedule_kernels_override ? schedule_kernels : 0); }
  update_analysis = 1;
  if (force_steady_statistical_data_overide){
    // This branch is to be entered only after tuning a space of algorithmic parameterizations, in which the expectation is that all kernels,
    //   both comm and comp, have reached a sufficiently-predictable state (steady state).
    for (auto it : comm_kernel_map){
      set_kernel_state(it.second,false);
      set_kernel_state_global(it.second,false);
    }
    for (auto it : comp_kernel_map){
      set_kernel_state(it.second,false);
      set_kernel_state_global(it.second,false);
    }
    update_analysis=0;
  }
}

void clear(){
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  // I don't see any reason to clear the communicator map. In fact, doing so would be harmful
  // Actually, the batch_maps will be empty anyways, as per the loops in 'final_accumulate'.
  // The kernel keys don't really need to be updated/cleared if the active/steady buffer logic isn't being used
  //   (which it isn't), so in fact the only relevant part of this routine is the clearing of the pathsets if necessary.
  for (auto& it : comp_kernel_map){
    bool is_clear=false;
    if (schedule_tag==""){
      active_kernels[it.second.val_index].clear_distribution();
      is_clear=true;
    }
    else if (schedule_tag=="cholinv"){
      if (clear_counter % 5 == 0){
        if (it.first.tag == 101){// potrf
          active_kernels[it.second.val_index].clear_distribution();
          is_clear=true;
        }
        else if (it.first.tag == 102){// trtri
          active_kernels[it.second.val_index].clear_distribution();
          is_clear=true;
        }
        else if (it.first.tag == 200){
          active_kernels[it.second.val_index].clear_distribution();
          is_clear=true;
        }
      }
      //No-op (for now, unless necessary)
    }
    else if (schedule_tag=="cacqr"){
      //No-op (for now, unless necessary)
    }
    else if (schedule_tag=="caqr_level1"){
      //No-op (for now, unless necessary)
    }
    else if (schedule_tag=="caqr_level1pipe"){
      //No-op (for now, unless necessary)
    }
  }
  for (auto& it : comm_kernel_map){
    bool is_clear=false;
    if (schedule_tag==""){
      active_kernels[it.second.val_index].clear_distribution();
      is_clear=true;
    }
    else if (schedule_tag=="cholinv"){
      // Note that I would like to check the reset count before clearing the Allgather, but I only need to do this once anyways.
      if (clear_counter == 10){
        if (world_size==64){
          if ((it.first.tag == 5) && (it.first.dim_sizes[0] == 16) && (it.first.dim_strides[0]==4) && (it.first.partner_offset == INT_MIN)){
            active_kernels[it.second.val_index].clear_distribution();
            is_clear=true;
          }
        } else if (world_size==512){
          if ((it.first.tag == 5) && (it.first.dim_sizes[0] == 64) && (it.first.dim_strides[0]==8) && (it.first.partner_offset == INT_MIN)){
            active_kernels[it.second.val_index].clear_distribution();
            is_clear=true;
          }
        }
      }
    }
    else if (schedule_tag=="cacqr"){
      //No-op (for now, unless necessary)
    }
    else if (schedule_tag=="caqr_level1"){
      //No-op (for now, unless necessary)
    }
    else if (schedule_tag=="caqr_level1pipe"){
      //No-op (for now, unless necessary)
    }
  }
}

void finalize(){
  // 'spf' deletion should occur automatically at program shutdown
  for (auto it : aggregate_channel_map){
    free(it.second);
  }
  if (std::getenv("CRITTER_VIZ_FILE") != NULL){
    if (is_world_root){
      stream.close();
      stream_tune.close();
      stream_reconstruct.close();
    }
  }
}

}
}
}
