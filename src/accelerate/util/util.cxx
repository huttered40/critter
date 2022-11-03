#include <limits.h>

#include "util.h"
#include "../../skeletonize/util/util.h"
#include "../../profile/util/util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace accelerate{

std::ofstream stream,stream_comm_kernel,stream_comp_kernel,stream_tune,stream_reconstruct;
bool global_schedule_decision;
int tuning_delta,reset_distribution_mode,reset_state_mode;
int sample_constraint_mode,schedule_kernels,update_analysis;
int stop_criterion_mode,debug_iter_count;
int comp_kernel_transfer_id,comm_kernel_transfer_id;
int delay_state_update,collective_state_protocol;
int comm_sample_aggregation_mode,comm_state_aggregation_mode;
int comp_sample_aggregation_mode,comp_state_aggregation_mode;
float kernel_time_limit,kernel_error_limit,kernel_percentage_limit;
float* save_path_data;
MPI_Request save_prop_req;
volatile double comp_start_time;
size_t kernel_count_limit,mode_1_width,mode_2_width;
size_t num_cp_measures,num_pp_measures;
size_t num_vol_measures,num_tracker_cp_measures;
size_t num_tracker_pp_measures,num_tracker_vol_measures;
size_t cp_costs_size,pp_costs_size,vol_costs_size;
int internal_tag,internal_tag1,internal_tag2,internal_tag3,internal_tag4,internal_tag5;
bool is_first_iter;
MPI_Datatype kernel_type,batch_type;
std::map<comm_kernel_key,kernel_key_id> comm_kernel_map;
std::map<comp_kernel_key,kernel_key_id> comp_kernel_map;
std::map<comm_kernel_key,kernel> comm_kernel_save_map;
std::map<comp_kernel_key,kernel> comp_kernel_save_map;
std::map<comm_kernel_key,kernel> comm_kernel_ref_map;
std::map<comp_kernel_key,kernel> comp_kernel_ref_map;
std::vector<comm_kernel_key> active_comm_kernel_keys;
std::vector<comp_kernel_key> active_comp_kernel_keys;
std::vector<kernel> active_kernels;
sample_propagation_forest spf;
std::map<comm_kernel_key,std::vector<kernel_batch>> comm_batch_map;
std::map<comp_kernel_key,std::vector<kernel_batch>> comp_batch_map;
std::map<comm_kernel_key,std::vector<kernel>> comm_kernel_list;
std::map<comp_kernel_key,std::vector<kernel>> comp_kernel_list;
std::vector<float> intercept_overhead;
std::vector<float> global_intercept_overhead;
std::vector<float> global_comp_kernel_stats;
std::vector<float> global_comm_kernel_stats;
std::vector<float> local_comp_kernel_stats;
std::vector<float> local_comm_kernel_stats;
std::vector<float> save_comp_kernel_stats;
std::vector<float> save_comm_kernel_stats;
std::vector<float> cp_costs;
std::vector<float> cp_costs_foreign;
std::vector<float> max_pp_costs;
std::vector<float> vol_costs;
std::vector<float> cp_costs_ref;
std::vector<float> max_pp_costs_ref;
std::vector<float> vol_costs_ref;
std::vector<char> eager_pad;
std::map<std::string,std::vector<float>> save_info;

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
  // No reason to reset distribution members ('M1','M2')
  // Do not reset the schedule_count members because those are needed to determine sample variance,
  //   which builds across schedules.
  if (reset_state_mode){
    this->steady_state = false;
    this->global_steady_state = false;
  }
  this->num_local_schedules=0;
  this->total_local_exec_time = 0;
  this->num_local_scheduled_units = 0;
}
void kernel::clear_distribution(){
  // All other members would have already been reset (via 'reset' above)
  this->hash_id = 0;
  this->registered_channels.clear();
  this->steady_state = false;
  this->global_steady_state = false;
  this->total_exec_time = 0;
  this->num_scheduled_units = 0;
  this->num_non_scheduled_units = 0;
  this->num_schedules = 0;
  this->num_non_schedules = 0;
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
    float M1_2 = active_batches[i].M1;
    float M2_2 = active_batches[i].M2;
    int n2 = active_batches[i].num_schedules;
    //assert(n2>0);
    float delta = this->M1 - M1_2;
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
  if (skeletonize::comm_kernel_map.find(key) != skeletonize::comm_kernel_map.end()){
    skel_val = skeletonize::active_kernels[skeletonize::comm_kernel_map[key].val_index];
    skel_val = std::max(skel_val,1);
  }
  return skel_val;
}
int get_skel_count(const comp_kernel_key& key){
  assert(sample_constraint_mode==3);
  int skel_val = 1;
  if (skeletonize::comp_kernel_map.find(key) != skeletonize::comp_kernel_map.end()){
    skel_val = skeletonize::active_kernels[skeletonize::comp_kernel_map[key].val_index];
    skel_val = std::max(skel_val,1);
  }
  return skel_val;
}


float get_estimate(const kernel& p, int analysis_param, float unit_count){
  if (analysis_param == 0){// arithmetic mean
    return get_arithmetic_mean(p);
  } else{
    return unit_count*get_harmonic_mean(p);
  }
}
float get_estimate(const kernel_propagate& p, int analysis_param, float unit_count){
  if (analysis_param == 0){// arithmetic mean
    return get_arithmetic_mean(p);
  } else{
    return unit_count*get_harmonic_mean(p);
  }
}
float get_estimate(const kernel_key_id& index, int analysis_param, float unit_count){
  if (analysis_param == 0){// arithmetic mean
    return get_arithmetic_mean(index);
  } else{
    return unit_count*get_harmonic_mean(index);
  }
}
float get_estimate(const kernel_key_id& index, const std::vector<kernel_batch>& active_batches, int analysis_param, float unit_count){
  auto stats = intermediate_stats(index,active_batches);
  if (analysis_param == 0){// arithmetic mean
    return stats.M1;
  } else{
    return unit_count*(1./stats.M1);
  }
}
float get_estimate(const intermediate_stats& p, int analysis_param, float unit_count){
  if (analysis_param == 0){// arithmetic mean
    return get_arithmetic_mean(p);
  } else{
    return unit_count*get_harmonic_mean(p);
  }
}

float get_arithmetic_mean(const kernel& p){
  // returns arithmetic mean
  return p.M1;
}
float get_arithmetic_mean(const kernel_propagate& p){
  // returns arithmetic mean
  return p.M1;
}
float get_arithmetic_mean(const kernel_key_id& index){
  // returns arithmetic mean
  return active_kernels[index.val_index].M1;
}
float get_arithmetic_mean(const intermediate_stats& p){
  // returns arithmetic mean
  return p.M1;
}

float get_harmonic_mean(const kernel& p){
  // returns arithmetic mean
  return 1./p.M1;
}
float get_harmonic_mean(const kernel_propagate& p){
  // returns arithmetic mean
  return 1./p.M1;
}
float get_harmonic_mean(const kernel_key_id& index){
  // returns arithmetic mean
  return 1./active_kernels[index.val_index].M1;
}
float get_harmonic_mean(const intermediate_stats& p){
  // returns arithmetic mean
  return 1./p.M1;
}

float get_variance(const kernel& p, int analysis_param){
  // returns variance
  size_t n = p.num_schedules;
  if (n<=1) return 1000000.;
  if (analysis_param == 0){
    return p.M2 / (n-1.);
  } else{
    return 1./p.M2 / (n-1.);
  }
}
float get_variance(const kernel_propagate& p, int analysis_param){
  // returns variance
  size_t n = p.num_schedules;
  if (n<=1) return 1000000.;
  if (analysis_param == 0){
    return p.M2 / (n-1.);
  } else{
    return 1./p.M2 / (n-1.);
  }
}
float get_variance(const kernel_key_id& index, int analysis_param){
  // returns variance
  size_t n = active_kernels[index.val_index].num_schedules;
  if (n<=1) return 1000000.;
  if (analysis_param == 0){
    return active_kernels[index.val_index].M2 / (n-1.);
  } else{
    return 1./active_kernels[index.val_index].M2 / (n-1.);
  }
}
float get_variance(const intermediate_stats& p, int analysis_param){
  // returns variance
  size_t n = p.num_schedules;
  if (n<=1) return 1000000.;
  if (analysis_param == 0){
    return p.M2 / (n-1.);
  } else{
    return 1./p.M2 / (n-1.);
  }
}

float get_std_dev(const kernel& p, int analysis_param){
  // returns variance
  return pow(get_variance(p,analysis_param),1./2.);
}
float get_std_dev(const kernel_propagate& p, int analysis_param){
  // returns variance
  return pow(get_variance(p,analysis_param),1./2.);
}
float get_std_dev(const kernel_key_id& index, int analysis_param){
  // returns variance
  return pow(get_variance(index,analysis_param),1./2.);
}
float get_std_dev(const intermediate_stats& p, int analysis_param){
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
      n = p.num_local_schedules; break;
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
      n = p.num_local_schedules; break;
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
      n = active_kernels[index.val_index].num_local_schedules; break;
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
      n = p.num_local_schedules; break;
  }
  n = std::max(1,n);
  return n;
}

// Standard error calculation takes into account execution path analysis
float get_std_error(const comm_kernel_key& key, const kernel& p, int analysis_param){
  // returns standard error
  int n = get_std_error_count(p);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(p,analysis_param) / pow(n*1.,1./2.);
}
float get_std_error(const comp_kernel_key& key, const kernel& p, int analysis_param){
  // returns standard error
  int n = get_std_error_count(p);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(p,analysis_param) / pow(n*1.,1./2.);
}
float get_std_error(const comm_kernel_key& key, const kernel_propagate& p, int analysis_param){
  // returns standard error
  int n = get_std_error_count(p);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(p,analysis_param) / pow(n*1.,1./2.);
}
float get_std_error(const comp_kernel_key& key, const kernel_propagate& p, int analysis_param){
  // returns standard error
  int n = get_std_error_count(p);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(p,analysis_param) / pow(n*1.,1./2.);
}
float get_std_error(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param){
  // returns standard error
  int n = get_std_error_count(index);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(index,analysis_param) / pow(n*1.,1./2.);
}
float get_std_error(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param){
  // returns standard error
  int n = get_std_error_count(index);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(index,analysis_param) / pow(n*1.,1./2.);
}
float get_std_error(const comm_kernel_key& key, const intermediate_stats& p, int analysis_param){
  // returns standard error
  int n = get_std_error_count(p);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(p,analysis_param) / pow(n*1.,1./2.);
}
float get_std_error(const comp_kernel_key& key, const intermediate_stats& p, int analysis_param){
  // returns standard error
  int n = get_std_error_count(p);
  if (sample_constraint_mode == 3) n = get_skel_count(key);
  return get_std_dev(p,analysis_param) / pow(n*1.,1./2.);
}

float get_confidence_interval(const comm_kernel_key& key, const kernel& p, int analysis_param, float level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,p,analysis_param);
}
float get_confidence_interval(const comm_kernel_key& key, const kernel_propagate& p, int analysis_param, float level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,p,analysis_param);
}
float get_confidence_interval(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param, float level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,index,analysis_param);
}
float get_confidence_interval(const comm_kernel_key& key, const intermediate_stats& p, int analysis_param, float level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,p,analysis_param);
}
float get_confidence_interval(const comp_kernel_key& key, const kernel& p, int analysis_param, float level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,p,analysis_param);
}
float get_confidence_interval(const comp_kernel_key& key, const kernel_propagate& p, int analysis_param, float level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,p,analysis_param);
}
float get_confidence_interval(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param, float level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,index,analysis_param);
}
float get_confidence_interval(const comp_kernel_key& key, const intermediate_stats& p, int analysis_param, float level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(key,p,analysis_param);
}

bool is_steady(const comm_kernel_key& key, const kernel& p, int analysis_param){
  if (sample_constraint_mode==-1) { assert(p.num_schedules < 2); return true; }
  if (stop_criterion_mode == 1){
    int __skel__count = (sample_constraint_mode == 3 ? get_skel_count(key) : 0);
    return ((get_confidence_interval(key,p,analysis_param) / get_estimate(p,analysis_param)) < kernel_error_limit) &&
            (p.num_schedules >= kernel_count_limit) &&
            (p.num_schedules >= (kernel_percentage_limit * __skel__count * debug_iter_count)) &&
            (p.total_exec_time >= kernel_time_limit);
  }
  else if (stop_criterion_mode == 0){
    int __skel__count = (sample_constraint_mode == 3 ? get_skel_count(key) : 0);
    return ((p.num_schedules>0) && (p.num_schedules >= (kernel_percentage_limit * __skel__count * debug_iter_count)));
  }
  else if (stop_criterion_mode == 2){
    comm_kernel_list[key].push_back(p);
    return false;
  }
}
bool is_steady(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param){
  if (sample_constraint_mode==-1) { assert(active_kernels[index.val_index].num_schedules < 2); return true; }
  if (stop_criterion_mode == 1){
    int __skel__count = (sample_constraint_mode == 3 ? get_skel_count(key) : 0);
    return ((get_confidence_interval(key,index,analysis_param) / get_estimate(index,analysis_param)) < kernel_error_limit) &&
            (active_kernels[index.val_index].num_schedules >= kernel_count_limit) &&
            (active_kernels[index.val_index].num_schedules >= (kernel_percentage_limit * __skel__count * debug_iter_count)) &&
            (active_kernels[index.val_index].total_exec_time >= kernel_time_limit);
  }
  else if (stop_criterion_mode == 0){
    int __skel__count = (sample_constraint_mode == 3 ? get_skel_count(key) : 0);
    return ((active_kernels[index.val_index].num_schedules>0) && (active_kernels[index.val_index].num_schedules >= (kernel_percentage_limit * __skel__count * debug_iter_count)));
  }
  else if (stop_criterion_mode == 2){
    comm_kernel_list[key].push_back(active_kernels[index.val_index]);
    return false;
  }
}
bool is_steady(const comm_kernel_key& key, const intermediate_stats& p, int analysis_param){
  if (sample_constraint_mode==-1) { assert(p.num_schedules < 2); return true; }
  if (stop_criterion_mode == 1){
    int __skel__count = (sample_constraint_mode == 3 ? get_skel_count(key) : 0);
    return ((get_confidence_interval(key,p,analysis_param) / get_estimate(p,analysis_param)) < kernel_error_limit) &&
            (p.num_schedules >= kernel_count_limit) &&
            (p.num_schedules >= (kernel_percentage_limit * __skel__count * debug_iter_count)) &&
            (p.total_exec_time >= kernel_time_limit);
  }
  else if (stop_criterion_mode == 0){
    int __skel__count = (sample_constraint_mode == 3 ? get_skel_count(key) : 0);
    return ((p.num_schedules>0) && (p.num_schedules >= (kernel_percentage_limit * __skel__count * debug_iter_count)));
  }
  else if (stop_criterion_mode == 2){
    return false;
  }
}
bool is_steady(const comp_kernel_key& key, const kernel& p, int analysis_param){
  if (sample_constraint_mode==-1) { assert(p.num_schedules < 2); return true; }
  if (stop_criterion_mode == 1){
    int __skel__count = (sample_constraint_mode == 3 ? get_skel_count(key) : 0);
    return ((get_confidence_interval(key,p,analysis_param) / get_estimate(p,analysis_param)) < kernel_error_limit) &&
            (p.num_schedules >= kernel_count_limit) &&
            (p.num_schedules >= (kernel_percentage_limit * __skel__count * debug_iter_count)) &&
            (p.total_exec_time >= kernel_time_limit);
  }
  else if (stop_criterion_mode == 0){
    int __skel__count = (sample_constraint_mode == 3 ? get_skel_count(key) : 0);
    return ((p.num_schedules>0) && (p.num_schedules >= (kernel_percentage_limit * __skel__count * debug_iter_count)));
  }
  else if (stop_criterion_mode == 2){
    comp_kernel_list[key].push_back(p);
    return false;
  }
}
bool is_steady(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param){
  if (sample_constraint_mode==-1) { assert(active_kernels[index.val_index].num_schedules < 2); return true; }
  if (stop_criterion_mode == 1){
    int __skel__count = (sample_constraint_mode == 3 ? get_skel_count(key) : 0);
    return ((get_confidence_interval(key,index,analysis_param) / get_estimate(index,analysis_param)) < kernel_error_limit) &&
            (active_kernels[index.val_index].num_schedules >= kernel_count_limit) &&
            (active_kernels[index.val_index].num_schedules >= (kernel_percentage_limit * __skel__count * debug_iter_count)) &&
            (active_kernels[index.val_index].total_exec_time >= kernel_time_limit);
  }
  else if (stop_criterion_mode == 0){
    int __skel__count = (sample_constraint_mode == 3 ? get_skel_count(key) : 0);
    return ((active_kernels[index.val_index].num_schedules>0) && (active_kernels[index.val_index].num_schedules >= (kernel_percentage_limit * __skel__count * debug_iter_count)));
  }
  else if (stop_criterion_mode == 2){
    comp_kernel_list[key].push_back(active_kernels[index.val_index]);
    return false;
  }
}
bool is_steady(const comp_kernel_key& key, const intermediate_stats& p, int analysis_param){
  if (sample_constraint_mode==-1) { assert(p.num_schedules < 2); return true; }
  if (stop_criterion_mode == 1){
    int __skel__count = (sample_constraint_mode == 3 ? get_skel_count(key) : 0);
    return ((get_confidence_interval(key,p,analysis_param) / get_estimate(p,analysis_param)) < kernel_error_limit) &&
            (p.num_schedules >= kernel_count_limit) &&
            (p.num_schedules >= (kernel_percentage_limit * __skel__count * debug_iter_count)) &&
            (p.total_exec_time >= kernel_time_limit);
  }
  else if (stop_criterion_mode == 0){
    int __skel__count = (sample_constraint_mode == 3 ? get_skel_count(key) : 0);
    return ((p.num_schedules>0) && (p.num_schedules >= (kernel_percentage_limit * __skel__count * debug_iter_count)));
  }
  else if (stop_criterion_mode == 2){
    return false;
  }
}

float get_error_estimate(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param){
  auto& active_batches = comm_batch_map[key];
  auto stats = intermediate_stats(index,active_batches);
  return get_confidence_interval(key,stats,analysis_param) / (get_estimate(stats,analysis_param));
}
float get_error_estimate(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param){
  auto& active_batches = comp_batch_map[key];
  auto stats = intermediate_stats(index,active_batches);
  return get_confidence_interval(key,stats,analysis_param) / (get_estimate(stats,analysis_param));
}
float get_error_estimate(const comm_kernel_key& key, const kernel_propagate& p, int analysis_param){
  return get_confidence_interval(key,p,analysis_param) / (get_estimate(p,analysis_param));
}
float get_error_estimate(const comp_kernel_key& key, const kernel_propagate& p, int analysis_param){
  return get_confidence_interval(key,p,analysis_param) / (get_estimate(p,analysis_param));
}

bool steady_test(const comm_kernel_key& key, const kernel& p, int analysis_param){
  if (!is_key_skipable(key)) return false;
  return is_steady(key,p,analysis_param);
}
bool steady_test(const comm_kernel_key& key, const kernel_key_id& index, int analysis_param){
  if (!is_key_skipable(key)) return false;
  return is_steady(key,index,analysis_param);
}
bool steady_test(const comp_kernel_key& key, const kernel& p, int analysis_param){
  if (!is_key_skipable(key)) return false;
  return is_steady(key,p,analysis_param);
}
bool steady_test(const comp_kernel_key& key, const kernel_key_id& index, int analysis_param){
  if (!is_key_skipable(key)) return false;
  return is_steady(key,index,analysis_param);
}

void update_kernel_stats(kernel& p, int analysis_param, volatile float exec_time, float unit_count){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  if (exec_time == 0) { exec_time=1.e-9; }
  if (p.global_steady_state == 0){
    p.num_schedules++;
    p.num_local_schedules++;
    p.num_scheduled_units += unit_count;
    p.num_local_scheduled_units += unit_count;
    p.total_exec_time += exec_time;
    p.total_local_exec_time += exec_time;
    // Online computation of up to 4th-order central moments using compunication time samples
    size_t n1 = p.num_schedules-1;
    size_t n = p.num_schedules;
    float x;
    if (analysis_param == 0){x = exec_time; }	// prep for arithmetic mean
    else                   {x = (unit_count>0 ? unit_count : 1.)/exec_time; }	// prep for harmonic mean
    float delta = x - p.M1;
    float delta_n = delta / n;
    float delta_n2 = delta_n*delta_n;
    float term1 = delta*delta_n*n1;
    p.M1 += delta_n;
    p.M2 += term1;
  }
  else{
    p.num_non_schedules++;
    p.num_non_scheduled_units += unit_count;
  }
}
void update_kernel_stats(const kernel_key_id& index, int analysis_param, volatile float exec_time, float unit_count){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  if (exec_time == 0) { exec_time=1.e-9; }
  if (active_kernels[index.val_index].global_steady_state == 0){
    active_kernels[index.val_index].num_schedules++;
    active_kernels[index.val_index].num_local_schedules++;
    active_kernels[index.val_index].num_scheduled_units += unit_count;
    active_kernels[index.val_index].num_local_scheduled_units += unit_count;
    active_kernels[index.val_index].total_exec_time += exec_time;
    active_kernels[index.val_index].total_local_exec_time += exec_time;
    // Online computation of up to 4th-order central moments using compunication time samples
    size_t n1 = active_kernels[index.val_index].num_schedules-1;
    size_t n = active_kernels[index.val_index].num_schedules;
    float x;
    if (analysis_param == 0){x = exec_time; }	// prep for arithmetic mean
    else                   {x = (unit_count>0 ? unit_count : 1.)/exec_time; }	// prep for harmonic mean
    float delta = x - active_kernels[index.val_index].M1;
    float delta_n = delta / n;
    float delta_n2 = delta_n*delta_n;
    float term1 = delta*delta_n*n1;
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
  float delta = dest.M1 - src.M1;
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

void update_kernel_stats(kernel_batch& batch, int analysis_param, volatile float exec_time, float unit_count){
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
  float x;
  if (analysis_param == 0){x = exec_time; }	// prep for arithmetic mean
  else                   {x = (unit_count>0 ? unit_count : 1.)/exec_time; }	// prep for harmonic mean
  float delta = x - batch.M1;
  float delta_n = delta / n;
  float delta_n2 = delta_n*delta_n;
  float term1 = delta*delta_n*n1;
  batch.M1 += delta_n;
  batch.M2 += term1;
}
void update_kernel_stats(kernel& dest, const kernel_batch& src, int analysis_param){
  // This function will implement the parallel algorithm computing the mean and variance of two partitions
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  // Online computation of up to 4th-order central moments using compunication time samples
  size_t n1 = dest.num_schedules;
  size_t n2 = src.num_schedules;
  float delta = dest.M1 - src.M1;
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
  float delta = dest.M1 - src.M1;
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
  mode_1_width = 25; mode_2_width = 15;
  communicator_count=INT_MIN;// to avoid conflict with p2p, which could range from (-p,p)
  internal_tag = 31133; internal_tag1 = internal_tag+1;
  internal_tag2 = internal_tag+2; internal_tag3 = internal_tag+3;
  internal_tag4 = internal_tag+4; internal_tag5 = internal_tag+5;
  is_first_iter = true;
  intercept_overhead.resize(3,0);
  global_intercept_overhead.resize(3,0);
  global_comm_kernel_stats.resize(5,0);
  global_comp_kernel_stats.resize(5,0);
  local_comm_kernel_stats.resize(5,0);
  local_comp_kernel_stats.resize(5,0);
  save_comp_kernel_stats.resize(2,0);
  save_comm_kernel_stats.resize(2,0);

  generate_initial_aggregate();

  kernel_propagate ex_3;
  MPI_Datatype kernel_internal_type[2] = { MPI_INT, MPI_FLOAT };
  int kernel_internal_block_len[2] = { 3,6 };
  MPI_Aint kernel_internal_disp[2] = { (char*)&ex_3.hash_id-(char*)&ex_3, (char*)&ex_3.num_scheduled_units-(char*)&ex_3 };
  PMPI_Type_create_struct(2,kernel_internal_block_len,kernel_internal_disp,kernel_internal_type,&kernel_type);
  PMPI_Type_commit(&kernel_type);

  // Reset these global variables, as some are updated by function arguments for convenience
  if (std::getenv("CRITTER_VIZ_FILE") != NULL){
    std::string stream_name = std::getenv("CRITTER_VIZ_FILE");
    std::string stream_tune_name = std::getenv("CRITTER_VIZ_FILE");
    std::string stream_reconstruct_name = std::getenv("CRITTER_VIZ_FILE");
    std::string stream_comm_kernel_name = std::getenv("CRITTER_VIZ_FILE");
    std::string stream_comp_kernel_name = std::getenv("CRITTER_VIZ_FILE");
    stream_name += "_accelerate.txt";
    stream_tune_name += "_accelerate_tune.txt";
    stream_reconstruct_name += "_accelerate_reconstruct.txt";
    stream_comm_kernel_name += "_commk.txt";
    stream_comp_kernel_name += "_compk.txt";
    if (is_world_root){
      stream.open(stream_name.c_str(),std::ofstream::app);
      stream_tune.open(stream_tune_name.c_str(),std::ofstream::app);
      stream_reconstruct.open(stream_reconstruct_name.c_str(),std::ofstream::app);
      stream_comm_kernel.open(stream_comm_kernel_name.c_str(),std::ofstream::app);
      stream_comp_kernel.open(stream_comp_kernel_name.c_str(),std::ofstream::app);
    }
  }
  if (std::getenv("CRITTER_DELTA") != NULL){
    tuning_delta = atoi(std::getenv("CRITTER_DELTA"));
    assert(tuning_delta>=0 && tuning_delta<=3);// tuning_delta==0 requires decomposition mechanism
  } else{
    tuning_delta = 0;
  }
  if (std::getenv("CRITTER_RESET_KERNEL_DISTRIBUTION") != NULL){
    reset_distribution_mode = atoi(std::getenv("CRITTER_RESET_KERNEL_DISTRIBUTION"));
    assert(reset_distribution_mode >=0 && reset_distribution_mode <=1);
  } else{
    reset_distribution_mode = 1;
  }
  if (std::getenv("CRITTER_RESET_KERNEL_STATE") != NULL){
    reset_state_mode = atoi(std::getenv("CRITTER_RESET_KERNEL_STATE"));
    assert(reset_state_mode >=0 && reset_state_mode <=1);
  } else{
    reset_state_mode = 1;
  }
  if (std::getenv("CRITTER_KERNEL_COUNT_LIMIT") != NULL){
    kernel_count_limit = atoi(std::getenv("CRITTER_KERNEL_COUNT_LIMIT"));
    assert(kernel_count_limit >= 0);
  } else{
    kernel_count_limit = 2;
  }
  if (std::getenv("CRITTER_KERNEL_TIME_LIMIT") != NULL){
    kernel_time_limit = atof(std::getenv("CRITTER_KERNEL_TIME_LIMIT"));
    assert(kernel_time_limit >= 0.);
  } else{
    kernel_time_limit = 0;
  }
  if (std::getenv("CRITTER_KERNEL_ERROR_LIMIT") != NULL){
    kernel_error_limit = atof(std::getenv("CRITTER_KERNEL_ERROR_LIMIT"));
    assert(kernel_error_limit >= 0.);
  } else{
    kernel_error_limit = 0;
  }
  if (std::getenv("CRITTER_KERNEL_PERCENTAGE_LIMIT") != NULL){
    kernel_percentage_limit = atof(std::getenv("CRITTER_KERNEL_PERCENTAGE_LIMIT"));
    assert(kernel_percentage_limit >= 0.);
  } else{
    kernel_percentage_limit = .001;
  }
  if (std::getenv("CRITTER_COMM_STATE_AGGREGATION_MODE") != NULL){
    comm_state_aggregation_mode = atoi(std::getenv("CRITTER_COMM_STATE_AGGREGATION_MODE"));
    assert(comm_state_aggregation_mode>=0 && comm_state_aggregation_mode<=2);
  } else{
    comm_state_aggregation_mode = 0;
  }
  if (std::getenv("CRITTER_COMP_STATE_AGGREGATION_MODE") != NULL){
    comp_state_aggregation_mode = atoi(std::getenv("CRITTER_COMP_STATE_AGGREGATION_MODE"));
    assert(comp_state_aggregation_mode>=0 && comp_state_aggregation_mode<=2);
  } else{
    comp_state_aggregation_mode = 0;
  }
  if (std::getenv("CRITTER_COMP_KERNEL_TRANSFER_MODE") != NULL){
    comp_kernel_transfer_id = atoi(std::getenv("CRITTER_COMP_KERNEL_TRANSFER_MODE"));
    assert(comp_kernel_transfer_id >=0 && comp_kernel_transfer_id <=1);
  } else{
    comp_kernel_transfer_id = 0;
  }
  if (std::getenv("CRITTER_COMM_KERNEL_TRANSFER_MODE") != NULL){
    comm_kernel_transfer_id = atoi(std::getenv("CRITTER_COMM_KERNEL_TRANSFER_MODE"));
    assert(comm_kernel_transfer_id >=0 && comm_kernel_transfer_id <=1);
  } else{
    comm_kernel_transfer_id = 0;
  }
  if (std::getenv("CRITTER_SAMPLE_CONSTRAINT_MODE") != NULL){
    sample_constraint_mode = atoi(std::getenv("CRITTER_SAMPLE_CONSTRAINT_MODE"));
    assert(sample_constraint_mode>=-1 && sample_constraint_mode<=3);
  } else{
    sample_constraint_mode = 0;
  }
  if (std::getenv("CRITTER_STOP_CRIT_MODE") != NULL){
    stop_criterion_mode = atof(std::getenv("CRITTER_STOP_CRIT_MODE"));
    assert(stop_criterion_mode >=0 && stop_criterion_mode <= 2);
  } else{
    stop_criterion_mode = 1;// confidence interval length
  }
  if (std::getenv("CRITTER_DELAY_STATE_UPDATE") != NULL){
    delay_state_update = atof(std::getenv("CRITTER_DELAY_STATE_UPDATE"));
    assert(delay_state_update >= 0 && delay_state_update <= 1);
  } else{
    delay_state_update = 0;
  }
  if (std::getenv("CRITTER_COLLECTIVE_STATE_PROTOCOL") != NULL){
    collective_state_protocol = atof(std::getenv("CRITTER_COLLECTIVE_STATE_PROTOCOL"));
    assert(collective_state_protocol >= 0 && collective_state_protocol <= 1);
  } else{
    collective_state_protocol = 1;
  }
  comm_envelope_param = 0;
  comm_stat_range = 0;
  comm_analysis_param = 0;
  comp_envelope_param = 0;
  comp_stat_range = 0;
  comp_analysis_param = 0;
  if (std::getenv("CRITTER_COMM_ANALYSIS") != NULL){
    comm_analysis_param = atof(std::getenv("CRITTER_COMM_ANALYSIS"));
    assert(comm_analysis_param >= 0 && comm_analysis_param <= 1);
  }
  if (std::getenv("CRITTER_COMP_ANALYSIS") != NULL){
    comp_analysis_param = atof(std::getenv("CRITTER_COMP_ANALYSIS"));
    assert(comp_analysis_param >= 0 && comp_analysis_param <= 1);
  }
  if (std::getenv("CRITTER_COMM_STAT_RANGE") != NULL){
    comm_stat_range = atof(std::getenv("CRITTER_COMM_STAT_RANGE"));
    assert(comm_stat_range >= 0 && comm_stat_range <= 1);
  }
  if (std::getenv("CRITTER_COMP_STAT_RANGE") != NULL){
    comp_stat_range = atof(std::getenv("CRITTER_COMP_STAT_RANGE"));
    assert(comp_stat_range >= 0 && comp_stat_range <= 1);
  }

  debug_iter_count = 1;
  // If user wants percentage-based stopping criterion, force knowledge of kernel frequency via skeletonize
  if (stop_criterion_mode==0) assert(sample_constraint_mode == 3);

  // Communication kernel time, computation kernel time, computation time, execution time
  num_cp_measures = 4;
  num_pp_measures = 4;
  num_vol_measures = 4;

  // 8 key members in 'comp_kernel' and 8 key members in 'comm_kernel'
  // +1 for each is the schedule count itself
  cp_costs_size = num_cp_measures + 14;
  if (sample_constraint_mode == 2) cp_costs_size += 9*comp_kernel_select_count + 9*comm_kernel_select_count;
  pp_costs_size = num_pp_measures;
  vol_costs_size = num_vol_measures;

  cp_costs.resize(cp_costs_size);
  cp_costs_foreign.resize(cp_costs_size);
  max_pp_costs.resize(pp_costs_size);
  vol_costs.resize(vol_costs_size);
  cp_costs_ref.resize(cp_costs_size);
  max_pp_costs_ref.resize(pp_costs_size);
  vol_costs_ref.resize(vol_costs_size);


  int eager_msg_size;
  MPI_Pack_size(cp_costs_size,MPI_FLOAT,comm,&eager_msg_size);
  int eager_pad_size = MPI_BSEND_OVERHEAD;
  eager_pad_size += eager_msg_size;
  eager_pad.resize(eager_pad_size);
}

void init_symbol(std::vector<std::string>& symbols){}

void open_symbol(const char* symbol, double curtime){}

void close_symbol(const char* symbol, double curtime){}

void final_accumulate(MPI_Comm comm, double last_time){
  assert(nonblocking_internal_info.size() == 0);
  cp_costs[num_cp_measures-1]+=(last_time-computation_timer);	// update critical path runtime
  cp_costs[num_cp_measures-3]+=(last_time-computation_timer);	// update critical path computation time
  vol_costs[num_vol_measures-1]+=(last_time-computation_timer);			// update per-process execution time
  vol_costs[num_vol_measures-3]+=(last_time-computation_timer);			// update per-process execution time

  _wall_time = wall_timer[wall_timer.size()-1];

  float temp_costs[4+4+3+5+5+2+3+3];
  for (int i=0; i<18; i++){ temp_costs[11+i]=0; }

  accelerate::_MPI_Barrier.comm = MPI_COMM_WORLD;
  accelerate::_MPI_Barrier.partner1 = -1;
  accelerate::_MPI_Barrier.partner2 = -1;
  accelerate::_MPI_Barrier.save_comp_key.clear();
  accelerate::_MPI_Barrier.save_comm_key.clear();
  for (auto& it : comp_kernel_map){
    if ((active_kernels[it.second.val_index].steady_state==1) && (should_schedule(it.second)==1)){
      accelerate::_MPI_Barrier.save_comp_key.push_back(it.first);
      temp_costs[num_cp_measures+num_vol_measures+13]++;
    }

    temp_costs[num_cp_measures+num_vol_measures+3] += active_kernels[it.second.val_index].num_local_schedules;
    temp_costs[num_cp_measures+num_vol_measures+4] += active_kernels[it.second.val_index].num_non_schedules;
    temp_costs[num_cp_measures+num_vol_measures+5] += active_kernels[it.second.val_index].num_local_scheduled_units;
    temp_costs[num_cp_measures+num_vol_measures+6] += active_kernels[it.second.val_index].num_non_scheduled_units;
    temp_costs[num_cp_measures+num_vol_measures+7] += active_kernels[it.second.val_index].total_local_exec_time;
    temp_costs[num_cp_measures+num_vol_measures+15] += active_kernels[it.second.val_index].num_local_schedules;
    temp_costs[num_cp_measures+num_vol_measures+16] += active_kernels[it.second.val_index].num_local_scheduled_units;
    temp_costs[num_cp_measures+num_vol_measures+17] += active_kernels[it.second.val_index].total_local_exec_time;
  }
  for (auto& it : comm_kernel_map){
    if ((active_kernels[it.second.val_index].steady_state==1) && (should_schedule(it.second)==1)){
      accelerate::_MPI_Barrier.save_comm_key.push_back(it.first);
      temp_costs[num_cp_measures+num_vol_measures+14]++;
    }

    temp_costs[num_cp_measures+num_vol_measures+8] += active_kernels[it.second.val_index].num_local_schedules;
    temp_costs[num_cp_measures+num_vol_measures+9] += active_kernels[it.second.val_index].num_non_schedules;
    temp_costs[num_cp_measures+num_vol_measures+10] += active_kernels[it.second.val_index].num_local_scheduled_units;
    temp_costs[num_cp_measures+num_vol_measures+11] += active_kernels[it.second.val_index].num_non_scheduled_units;
    temp_costs[num_cp_measures+num_vol_measures+12] += active_kernels[it.second.val_index].total_local_exec_time;
    temp_costs[num_cp_measures+num_vol_measures+18] += active_kernels[it.second.val_index].num_local_schedules;
    temp_costs[num_cp_measures+num_vol_measures+19] += active_kernels[it.second.val_index].num_local_scheduled_units;
    temp_costs[num_cp_measures+num_vol_measures+20] += active_kernels[it.second.val_index].total_local_exec_time;
  }

  max_pp_costs = vol_costs;// copy over the per-process measurements that exist in vol_costs
  for (auto i=0; i<num_cp_measures; i++) temp_costs[i] = cp_costs[i];
  for (auto i=0; i<max_pp_costs.size(); i++) temp_costs[num_cp_measures+i] = max_pp_costs[i];
  temp_costs[num_cp_measures+num_vol_measures] = intercept_overhead[0];
  temp_costs[num_cp_measures+num_vol_measures+1] = intercept_overhead[1];
  temp_costs[num_cp_measures+num_vol_measures+2] = intercept_overhead[2];
  PMPI_Allreduce(MPI_IN_PLACE,&temp_costs[0],4+4+3+5+5+2+3+3,MPI_FLOAT,MPI_MAX,comm);
  for (auto i=0; i<num_cp_measures; i++) cp_costs[i] = temp_costs[i];
  for (auto i=0; i<max_pp_costs.size(); i++) max_pp_costs[i] = temp_costs[num_cp_measures+i];
  if (autotuning_debug == 0){
    global_intercept_overhead[0] += temp_costs[num_cp_measures+num_vol_measures];
    global_intercept_overhead[1] += temp_costs[num_cp_measures+num_vol_measures+1];
    global_intercept_overhead[2] += temp_costs[num_cp_measures+num_vol_measures+2];
    global_comp_kernel_stats[0] += temp_costs[num_cp_measures+num_vol_measures+3];
    global_comp_kernel_stats[1] = temp_costs[num_cp_measures+num_vol_measures+4];
    global_comp_kernel_stats[2] += temp_costs[num_cp_measures+num_vol_measures+5];
    global_comp_kernel_stats[3] = temp_costs[num_cp_measures+num_vol_measures+6];
    global_comp_kernel_stats[4] += temp_costs[num_cp_measures+num_vol_measures+7];
    global_comm_kernel_stats[0] += temp_costs[num_cp_measures+num_vol_measures+8];
    global_comm_kernel_stats[1] = temp_costs[num_cp_measures+num_vol_measures+9];
    global_comm_kernel_stats[2] += temp_costs[num_cp_measures+num_vol_measures+10];
    global_comm_kernel_stats[3] = temp_costs[num_cp_measures+num_vol_measures+11];
    global_comm_kernel_stats[4] += temp_costs[num_cp_measures+num_vol_measures+12];
    local_comp_kernel_stats[0] = temp_costs[num_cp_measures+num_vol_measures+15];
    local_comp_kernel_stats[1] = global_comp_kernel_stats[1] - save_comp_kernel_stats[0];
    local_comp_kernel_stats[2] = temp_costs[num_cp_measures+num_vol_measures+16];
    local_comp_kernel_stats[3] = global_comp_kernel_stats[3] - save_comp_kernel_stats[1];
    local_comp_kernel_stats[4] = temp_costs[num_cp_measures+num_vol_measures+17];
    local_comm_kernel_stats[0] = temp_costs[num_cp_measures+num_vol_measures+18];
    local_comm_kernel_stats[1] = global_comm_kernel_stats[1] - save_comm_kernel_stats[0];
    local_comm_kernel_stats[2] = temp_costs[num_cp_measures+num_vol_measures+19];
    local_comm_kernel_stats[3] = global_comm_kernel_stats[3] - save_comm_kernel_stats[1];
    local_comm_kernel_stats[4] = temp_costs[num_cp_measures+num_vol_measures+20];
    save_comp_kernel_stats[0] = global_comp_kernel_stats[1];
    save_comp_kernel_stats[1] = global_comp_kernel_stats[3];
    save_comm_kernel_stats[0] = global_comm_kernel_stats[1];
    save_comm_kernel_stats[1] = global_comm_kernel_stats[3];
  }
  PMPI_Allreduce(MPI_IN_PLACE,&vol_costs[0],vol_costs.size(),MPI_FLOAT,MPI_SUM,comm);
  accelerate::_MPI_Barrier.aggregate_comp_kernels = temp_costs[num_cp_measures+num_vol_measures+13]>0;
  accelerate::_MPI_Barrier.aggregate_comm_kernels = temp_costs[num_cp_measures+num_vol_measures+14]>0;
  accelerate::_MPI_Barrier.should_propagate = accelerate::_MPI_Barrier.aggregate_comp_kernels>0 || accelerate::_MPI_Barrier.aggregate_comm_kernels>0;
  // Just don't propagate the final kernels.
}

void reset(bool schedule_kernels_override, bool force_steady_statistical_data_overide){
  assert(nonblocking_internal_info.size() == 0);
  memset(&cp_costs[0],0,sizeof(float)*cp_costs.size());
  memset(&cp_costs_foreign[0],0,sizeof(float)*cp_costs.size());
  memset(&max_pp_costs[0],0,sizeof(float)*max_pp_costs.size());
  memset(&vol_costs[0],0,sizeof(float)*vol_costs.size());
  bsp_counter=0;
  memset(&intercept_overhead[0],0,sizeof(float)*intercept_overhead.size());
  memset(&local_comp_kernel_stats[0],0,sizeof(float)*local_comp_kernel_stats.size());
  memset(&local_comm_kernel_stats[0],0,sizeof(float)*local_comm_kernel_stats.size());

  // This reset will no longer reset the kernel state, but will reset the schedule counters
  for (auto& it : comp_kernel_map){
    it.second.is_active = false;
    active_kernels[it.second.val_index].reset();
  }
  for (auto& it : comm_kernel_map){
    it.second.is_active = false;
    active_kernels[it.second.val_index].reset();
  }

  if (std::getenv("CRITTER_MODE") != NULL){
    mode = atoi(std::getenv("CRITTER_MODE"));
    assert(mode >=0 && mode <=1);
  } else{
    mode = 1;
  }
  if (std::getenv("CRITTER_SCHEDULE_KERNELS") != NULL){
    schedule_kernels = atoi(std::getenv("CRITTER_SCHEDULE_KERNELS"));
    assert(schedule_kernels >=0 && schedule_kernels <=1);
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

void clear(int tag_count, int* distribution_tags){

  if (reset_distribution_mode==1){
    // I don't see any reason to clear the communicator map. In fact, doing so would be harmful
    // Actually, the batch_maps will be empty anyways, as per the loops in 'final_accumulate'.
    // The kernel keys don't really need to be updated/cleared if the active/steady buffer logic isn't being used
    //   (which it isn't), so in fact the only relevant part of this routine is the clearing of the pathsets if necessary.
    for (auto& it : comp_kernel_map){
      active_kernels[it.second.val_index].clear_distribution();
    }
    for (auto& it : comm_kernel_map){
      active_kernels[it.second.val_index].clear_distribution();
    }
    for (auto& it : comp_kernel_ref_map){
      it.second.clear_distribution();
    }
    for (auto& it : comm_kernel_ref_map){
      it.second.clear_distribution();
    }
  }
  else{
    for (auto& it : comp_kernel_map){
      for (int i=0; i<tag_count; i++){
        if (it.first.tag == distribution_tags[i]){
          active_kernels[it.second.val_index].clear_distribution();
          break;
        }
      }
    }
    for (auto& it : comm_kernel_map){
      for (int i=0; i<tag_count; i++){
        if (it.first.tag == distribution_tags[i]){
          active_kernels[it.second.val_index].clear_distribution();
          break;
        }
      }
    }
    for (auto& it : comp_kernel_ref_map){
      for (int i=0; i<tag_count; i++){
        if (it.first.tag == distribution_tags[i]){
          it.second.clear_distribution();
          break;
        }
      }
    }
    for (auto& it : comm_kernel_ref_map){
      for (int i=0; i<tag_count; i++){
        if (it.first.tag == distribution_tags[i]){
          it.second.clear_distribution();
          break;
        }
      }
    }
  }
}

void reference_initiate(){
  kernel_error_limit=0;
  stop_criterion_mode = 1;
  kernel_percentage_limit=1;
  // Save the existing distributions so that the reference variant doesn't erase or infect them
  // Reset the distributions so that the reference doesn't take the autotuned distribution from previous configurations
  //   if reset_kernel_distribution == 0
  for (auto& it : comp_kernel_map){
    comp_kernel_save_map[it.first] = active_kernels[it.second.val_index];
    if (comp_kernel_ref_map.find(it.first) != comp_kernel_ref_map.end()){
      active_kernels[it.second.val_index] = comp_kernel_ref_map[it.first];
    } else{
      active_kernels[it.second.val_index].reset();
      active_kernels[it.second.val_index].clear_distribution();
      active_kernels[it.second.val_index].num_non_schedules = 0;
      active_kernels[it.second.val_index].num_non_scheduled_units = 0;
    }
  }
  for (auto& it : comm_kernel_map){
    comm_kernel_save_map[it.first] = active_kernels[it.second.val_index];
    if (comm_kernel_ref_map.find(it.first) != comm_kernel_ref_map.end()){
      active_kernels[it.second.val_index] = comm_kernel_ref_map[it.first];
    } else{
      active_kernels[it.second.val_index].reset();
      active_kernels[it.second.val_index].clear_distribution();
      active_kernels[it.second.val_index].num_non_schedules = 0;
      active_kernels[it.second.val_index].num_non_scheduled_units = 0;
    }
  }
}

void reference_transfer(){
  if (std::getenv("CRITTER_KERNEL_ERROR_LIMIT") != NULL){ kernel_error_limit = atof(std::getenv("CRITTER_KERNEL_ERROR_LIMIT")); }
  if (std::getenv("CRITTER_KERNEL_PERCENTAGE_LIMIT") != NULL){ kernel_percentage_limit = atof(std::getenv("CRITTER_KERNEL_PERCENTAGE_LIMIT")); }
  if (std::getenv("CRITTER_STOP_CRIT_MODE") != NULL){ stop_criterion_mode = atof(std::getenv("CRITTER_STOP_CRIT_MODE")); }
  std::memcpy(&cp_costs_ref[0],&cp_costs[0],num_cp_measures*sizeof(float));
  std::memcpy(&max_pp_costs_ref[0],&max_pp_costs[0],num_pp_measures*sizeof(float));
  std::memcpy(&vol_costs_ref[0],&vol_costs[0],num_vol_measures*sizeof(float));
  for (auto& it : comp_kernel_map){
    comp_kernel_ref_map[it.first] = active_kernels[it.second.val_index];
    if (comp_kernel_save_map.find(it.first) != comp_kernel_save_map.end()){
      active_kernels[it.second.val_index] = comp_kernel_save_map[it.first];
    } else{
      active_kernels[it.second.val_index].reset();
      active_kernels[it.second.val_index].clear_distribution();
      active_kernels[it.second.val_index].num_non_schedules = 0;
      active_kernels[it.second.val_index].num_non_scheduled_units = 0;
    }
  }
  for (auto& it : comm_kernel_map){
    comm_kernel_ref_map[it.first] = active_kernels[it.second.val_index];
    if (comm_kernel_save_map.find(it.first) != comm_kernel_save_map.end()){
      active_kernels[it.second.val_index] = comm_kernel_save_map[it.first];
    } else{
      active_kernels[it.second.val_index].reset();
      active_kernels[it.second.val_index].clear_distribution();
      active_kernels[it.second.val_index].num_non_schedules = 0;
      active_kernels[it.second.val_index].num_non_scheduled_units = 0;
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
      stream_comm_kernel.close();
      stream_comp_kernel.close();
    }
  }
}

}
}
}
