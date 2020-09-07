#include <limits.h>

#include "util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace discretization{

int analysis_mode;
int is_optimized;
int autotuning_propagate;
int schedule_kernels;
int update_analysis;
int comm_sample_include_idle;
MPI_Datatype comm_pattern_key_type;
MPI_Datatype comp_pattern_key_type;
MPI_Datatype pattern_type;
size_t pattern_count_limit;
double pattern_time_limit;
double pattern_error_limit;
int comm_channel_tag_count;
std::map<comm_pattern_key,pattern_key_id> comm_pattern_map;
std::map<comp_pattern_key,pattern_key_id> comp_pattern_map;
std::vector<comm_pattern_key> steady_state_comm_pattern_keys;
std::vector<comm_pattern_key> active_comm_pattern_keys;
std::vector<comp_pattern_key> steady_state_comp_pattern_keys;
std::vector<comp_pattern_key> active_comp_pattern_keys;
std::vector<pattern> steady_state_patterns;
std::vector<pattern> active_patterns;
comm_pattern_key previous_comm_key;
std::map<std::pair<comm_pattern_key,comm_pattern_key>,idle_pattern> comm_pattern_pair_map;
sample_propagation_forest spf;
std::map<MPI_Comm,comm_channel_node*> comm_channel_map;
std::map<int,comm_channel_node*> p2p_channel_map;
std::map<comm_pattern_key,std::vector<pattern_batch>> comm_batch_map;
std::map<comp_pattern_key,std::vector<pattern_batch>> comp_batch_map;
std::vector<comm_channel_node*> intermediate_channels;

// ****************************************************************************************************************************************************
pattern::pattern(){
  this->total_exec_time = 0;
  this->steady_state=0;
  this->global_steady_state=0;
  this->num_schedules = 0;
  this->num_non_schedules = 0;
  this->num_scheduled_units = 0;
  this->num_non_scheduled_units = 0;
  this->num_propagations=0;
  this->num_non_propagations=0;
  this->M1=0; this->M2=0;
}

pattern::pattern(const pattern& _copy){
  this->total_exec_time = _copy.total_exec_time;
  this->steady_state = _copy.steady_state;
  this->global_steady_state = _copy.global_steady_state;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = _copy.num_non_schedules;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_non_scheduled_units = _copy.num_non_scheduled_units;
  this->num_propagations=_copy.num_propagations;
  this->num_non_propagations=_copy.num_non_propagations;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
}

pattern& pattern::operator=(const pattern& _copy){
  this->total_exec_time = _copy.total_exec_time;
  this->steady_state = _copy.steady_state;
  this->global_steady_state = _copy.global_steady_state;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = _copy.num_non_schedules;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_non_scheduled_units = _copy.num_non_scheduled_units;
  this->num_propagations=_copy.num_propagations;
  this->num_non_propagations=_copy.num_non_propagations;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  return *this;
}

// ****************************************************************************************************************************************************
pattern_batch::pattern_batch(comm_channel_node* node){
  this->channel_count=0;
  this->num_schedules = 0;
  this->M1=0; this->M2=0;
  this->open_channel_count=comm_channel_map.size() + p2p_channel_map.size();
  // Remove channels that are descendents and ancestors of 'node'
  if (node != nullptr){
    // TODO: Add to map those channels that are not feasible
    spf.fill_ancestors(node,*this);
    spf.fill_descendants(node,*this);
  }
  this->open_channel_count -= this->closed_channels.size();
}

pattern_batch::pattern_batch(const pattern_batch& _copy){
  this->channel_count = _copy.channel_count;
  this->num_schedules = _copy.num_schedules;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  this->open_channel_count = _copy.open_channel_count;
  this->closed_channels = _copy.closed_channels;
}

pattern_batch& pattern_batch::operator=(const pattern_batch& _copy){
  this->channel_count = _copy.channel_count;
  this->num_schedules = _copy.num_schedules;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  this->open_channel_count = _copy.open_channel_count;
  this->closed_channels = _copy.closed_channels;
  return *this;
}

// ****************************************************************************************************************************************************
idle_pattern::idle_pattern(){
  this->M1=0; this->M2=0;
  this->num_schedules = 0;
  this->num_non_schedules = 0;
}

idle_pattern::idle_pattern(const idle_pattern& _copy){
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = _copy.num_non_schedules;
}

idle_pattern& idle_pattern::operator=(const idle_pattern& _copy){
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = _copy.num_non_schedules;
  return *this;
}

// ****************************************************************************************************************************************************
pattern_key_id::pattern_key_id(bool _is_active, int _key_index, int _val_index, bool _is_updated){
  this->is_active = _is_active;
  this->key_index = _key_index;
  this->val_index = _val_index;
  this->is_updated = _is_updated;
}

pattern_key_id::pattern_key_id(const pattern_key_id& _copy){
  this->is_active = _copy.is_active;
  this->key_index = _copy.key_index;
  this->val_index = _copy.val_index;
  this->is_updated = _copy.is_updated;
}

pattern_key_id& pattern_key_id::operator=(const pattern_key_id& _copy){
  this->is_active = _copy.is_active;
  this->key_index = _copy.key_index;
  this->val_index = _copy.val_index;
  this->is_updated = _copy.is_updated;
  return *this;
}

// ****************************************************************************************************************************************************
comm_channel_node::comm_channel_node(){
  this->tag = -1;
  this->frequency=0;
  this->children.push_back(std::vector<comm_channel_node*>());
}

void sample_propagation_forest::generate_span(comm_channel_node* node, std::vector<std::pair<int,int>>& perm_tuples){
  // Assumed that perm_tuples define a permutation of comm_channel_nodes that are communicator siblings
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
bool sample_propagation_forest::is_child(comm_channel_node* tree_node, comm_channel_node* node){
  // First check that the root node 'node' is not a p2p, regardless of whether 'tree_node' is a subcomm or p2p
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  if (node->tag >= world_size) return false;
  // Returns true iff tree_node is an ancestor of node

  if (tree_node->tag < world_size){
    // Check all tuples
    for (auto i=0; i<tree_node->id.size(); i++){
      bool found_match=false;
      for (auto j=0; j<node->id.size(); j++){
        found_match = ((tree_node->id[i].second % node->id[j].second == 0) &&
                      (tree_node->id[i].first < node->id[j].first) &&
                      (tree_node->id[i].second < (node->id[j].first*node->id[j].second)) &&
                      ((tree_node->id[i].first*tree_node->id[i].second) <= (node->id[j].first*node->id[j].second)));
        if (found_match) break;
      }// As long as one match is found, tuple id[i] fits
      if (!found_match) return false;
    }
  } else{
    // p2p nodes need special machinery for detecting whether its a child of 'node'
    int rank = tree_node->offset[0] - node->offset[0];//TODO; Fix if tree_node describes a aggregated p2p (rare case, we haven't needed it yet)
    if (rank<0) return false;// corner case
    for (int i=node->id.size()-1; i>=0; i--){
      if (rank >= (node->id[i].first*node->id[i].second)){ return false; }
      if (i>0) { rank %= node->id[i].second; }
    }
    if (rank%node->id[0].second != 0) return false;
  }
  return true;
}
bool sample_propagation_forest::are_siblings(comm_channel_node* parent, int subtree_idx, std::vector<int>& skip_indices){
  // Perform recursive permutation generation to identify if a permutation of tuples among siblings is valid
  // Return true if parent's children are valid siblings
  if (parent->children[subtree_idx].size()<=1) return true;
  // Check if all are p2p, all are subcomms, or if there is a mixture.
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int p2p_count=0; int subcomm_count=0;
  for (auto i=0; i<parent->children[subtree_idx].size(); i++){
    if (parent->children[subtree_idx][i]->tag >= world_size){ p2p_count++; }
    else { subcomm_count++; }
  }
  //if ((subcomm_count>0) && (p2p_count>0)) return false;
  if (subcomm_count==0) return true;// all p2p siblings are fine
  std::vector<std::pair<int,int>> static_info;
  int skip_index=0;
  for (auto i=0; i<parent->children[subtree_idx].size(); i++){
    if ((skip_index<skip_indices.size()) && (i==skip_indices[skip_index])){
      skip_index++;
      continue;
    }
    if (parent->children[subtree_idx][i]->tag >= world_size){
      continue;
    }
    for (auto j=0; j<parent->children[subtree_idx][i]->id.size(); j++){
      static_info.push_back(parent->children[subtree_idx][i]->id[j]);
    }
  }
  std::vector<std::pair<int,int>> gen_info;
  std::vector<std::pair<int,int>> save_info;
  bool valid_siblings=false;
  generate_sibling_perm(static_info,gen_info,save_info,0,valid_siblings);
  if (valid_siblings){
    // Now we must run through the p2p, if any. I don't think we need to check if a p2p is also skipped via 'skip_indices', because if its a child of one of the communicators,
    //   it'll b a child of the generated span
    if (p2p_count>0){
      auto* span_node = new comm_channel_node;
      generate_span(span_node,save_info);
      for (auto i=0; i<parent->children[subtree_idx].size(); i++){
        if (parent->children[subtree_idx][i]->tag >= world_size){
          valid_siblings = this->is_child(parent->children[subtree_idx][i],span_node);
          if (!valid_siblings){
            delete span_node;
            return false;
          }
        }
      }
      delete span_node;
    }
  }
  return valid_siblings;
}
bool sample_propagation_forest::partition_test(comm_channel_node* parent, int subtree_idx){
  // Perform recursive permutation generation to identify if a permutation of tuples among siblings is valid
  // Return true if parent's children are valid siblings
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  std::vector<std::pair<int,int>> static_info;
  for (auto i=0; i<parent->children[subtree_idx].size(); i++){
    for (auto j=0; j<parent->children[subtree_idx][i]->id.size(); j++){
      if (parent->children[subtree_idx][i]->tag < world_size){
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
void sample_propagation_forest::find_parent(comm_channel_node* tree_root, comm_channel_node* tree_node, comm_channel_node*& parent){
  if (tree_root==nullptr) return;
  for (auto i=0; i<tree_root->children.size(); i++){
    for (auto j=0; j<tree_root->children[i].size(); j++){
      this->find_parent(tree_root->children[i][j],tree_node,parent);// Cannot be nullptrs. Nullptr children mean the children member is empty
    }
  }
  if ((parent==nullptr) && (this->is_child(tree_node,tree_root))){
    parent = tree_root;
  }
  return;
}
void sample_propagation_forest::fill_ancestors(comm_channel_node* tree_node, pattern_batch& batch){
  if (tree_node==nullptr) return;
  batch.closed_channels.insert(tree_node);
  this->fill_ancestors(tree_node->parent,batch);
}
void sample_propagation_forest::fill_descendants(comm_channel_node* tree_node, pattern_batch& batch){
  if (tree_node==nullptr) return;
  batch.closed_channels.insert(tree_node);
  for (auto i=0; i<tree_node->children.size(); i++){
    for (auto j=0; j<tree_node->children[i].size(); j++){
      this->fill_descendants(tree_node->children[i][j],batch);
    }
  }
}
void sample_propagation_forest::clear_tree_info(comm_channel_node* tree_root){
  if (tree_root==nullptr) return;
  for (auto i=0; i<tree_root->children.size(); i++){
    for (auto j=0; j<tree_root->children[i].size(); j++){
      this->clear_tree_info(tree_root->children[i][j]);// Cannot be nullptrs. Nullptr children mean the children member is empty
    }
  }
  tree_root->frequency=0;
  return;
}
void sample_propagation_forest::delete_tree(comm_channel_node*& tree_root){
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
sample_propagation_forest::sample_propagation_forest(){}
sample_propagation_forest::~sample_propagation_forest(){
  this->delete_tree(this->tree->root);
  free(this->tree);
}
void sample_propagation_forest::clear_info(){
  this->clear_tree_info(this->tree->root);
}
int sample_propagation_forest::translate_rank(MPI_Comm comm, int rank){
  // Returns new rank relative to world communicator
  auto node = comm_channel_map[comm];
  int new_rank = node->offset[0];
  for (auto i=0; i<node->id.size(); i++){
    new_rank += node->id[i].second*(rank%node->id[i].first);
    rank /= node->id[i].first;
  }
  return new_rank;
}
void sample_propagation_forest::insert_node(comm_channel_node* tree_node){
  // Fill in parent and children, and iterate over all trees of course.
  // Post-order traversal
  // Follow rules from paper to deduce first whether tree_node can be a child of the current parent.
  assert(tree_node != nullptr);
  bool is_comm = !(tree_node->id.size()==1 && tree_node->id[0].second==0);
  comm_channel_node* parent = nullptr;
  //TODO: I assume here that we care about the first SPT in the SPF. Figure out how to fix this later
  this->find_parent(this->tree->root,tree_node,parent);
  tree_node->parent = parent;
  assert(parent!=nullptr);
 
  // Try adding 'tree_node' to each SPT. If none fit, append parent's children array and add it there, signifying a new tree, rooted at 'parent'
  bool valid_parent = false;
  int save_tree_idx=-1;
  for (auto i=0; i<parent->children.size(); i++){
    parent->children[i].push_back(tree_node);
    std::vector<int> sibling_to_child_indices;
    for (auto j=0; j<parent->children[i].size()-1; j++){
      if (this->is_child(parent->children[i][j],tree_node)){
        sibling_to_child_indices.push_back(j);
      }
    }
    bool sibling_decision = this->are_siblings(parent,i,sibling_to_child_indices);
    if (sibling_decision){
      save_tree_idx=i;
      valid_parent=true;
      for (auto j=0; j<sibling_to_child_indices.size(); j++){
        tree_node->children[0].push_back(parent->children[i][sibling_to_child_indices[j]]);
        parent->children[i][sibling_to_child_indices[j]]->parent=tree_node;
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
    save_tree_idx=parent->children.size();
    parent->children.push_back(std::vector<comm_channel_node*>());
    parent->children[parent->children.size()-1].push_back(tree_node);
  }

  // sanity check for right now
  int world_comm_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_comm_rank);
  if (world_comm_rank==8){
    std::cout << "parent of tree_node{ " << tree_node->offset[0];
    for (auto i=0; i<tree_node->id.size(); i++){
      std::cout << " (" << tree_node->id[i].first << "," << tree_node->id[i].second << ")";
    }
    std::cout << " } is { " << parent->tag << " " << parent->offset[0];
    for (auto i=0; i<parent->id.size(); i++){
      std::cout << " (" << parent->id[i].first << "," << parent->id[i].second << ")";
    }
    std::cout << " }\n";
    for (auto i=0; i<parent->children.size(); i++){
      std::cout << "\tsubtree " << i << " contains " << parent->children[i].size() << " children\n";
      for (auto j=0; j<parent->children[i].size(); j++){
        std::cout << "\t\tchild " << j << " is { " << parent->children[i][j]->tag << " " << parent->children[i][j]->offset[0];
        for (auto k=0; k<parent->children[i][j]->id.size(); k++){
          std::cout << " (" << parent->children[i][j]->id[k].first << " " << parent->children[i][j]->id[k].second << ")";
        }
        std::cout << " }\n";
      }
    }
  }
}


// ****************************************************************************************************************************************************
bool is_key_skipable(const comm_pattern_key& key){
  // For now, only barriers cannot be skipped
  if (key.tag == 0){ return false; }
  return true;
}
bool is_key_skipable(const comp_pattern_key& key){
  return true;
}

double get_estimate(const pattern& p, int analysis_param, double unit_count){
  if (analysis_param == 0){// arithmetic mean
    return get_arithmetic_mean(p);
  } else{
    return unit_count*get_harmonic_mean(p);
  }
}
double get_estimate(const pattern_key_id& index, int analysis_param, double unit_count){
  if (analysis_param == 0){// arithmetic mean
    return get_arithmetic_mean(index);
  } else{
    return unit_count*get_harmonic_mean(index);
  }
}

double get_arithmetic_mean(const pattern& p){
  // returns arithmetic mean
  return p.M1;
}
double get_arithmetic_mean(const pattern_key_id& index){
  // returns arithmetic mean
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return pattern_list[index.val_index].M1;
}

double get_harmonic_mean(const pattern& p){
  // returns arithmetic mean
  return 1./p.M1;
}
double get_harmonic_mean(const pattern_key_id& index){
  // returns arithmetic mean
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return 1./pattern_list[index.val_index].M1;
}

double get_variance(const pattern& p, int analysis_param){
  // returns variance
  size_t n = p.num_schedules;
  if (n<=1) return 1000.;
  if (analysis_param == 0){
    return p.M2 / (n-1.);
  } else{
    return 1./p.M2 / (n-1.);
  }
}
double get_variance(const pattern_key_id& index, int analysis_param){
  // returns variance
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  size_t n = pattern_list[index.val_index].num_schedules;
  if (n<=1) return 1000.;
  if (analysis_param == 0){
    return pattern_list[index.val_index].M2 / (n-1.);
  } else{
    return 1./pattern_list[index.val_index].M2 / (n-1.);
  }
}

double get_std_dev(const pattern& p, int analysis_param){
  // returns variance
  return pow(get_variance(p,analysis_param),1./2.);
}
double get_std_dev(const pattern_key_id& index, int analysis_param){
  // returns variance
  return pow(get_variance(index,analysis_param),1./2.);
}

double get_std_error(const pattern& p, int analysis_param){
  // returns standard error
  size_t n = p.num_schedules;
  return get_std_dev(p,analysis_param) / pow(n*1.,1./2.);
}
double get_std_error(const pattern_key_id& index, int analysis_param){
  // returns standard error
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  size_t n = pattern_list[index.val_index].num_schedules;
  return get_std_dev(index,analysis_param) / pow(n*1.,1./2.);
}

double get_confidence_interval(const pattern& p, int analysis_param, double level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(p,analysis_param);
}
double get_confidence_interval(const pattern_key_id& index, int analysis_param, double level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(index,analysis_param);
}

bool is_steady(const pattern& p, int analysis_param){
  return ((get_confidence_interval(p,analysis_param) / (2.*get_estimate(p,analysis_param))) < pattern_error_limit) &&
          (p.num_schedules >= pattern_count_limit) &&
          (p.total_exec_time >= pattern_time_limit);
}
bool is_steady(const pattern_key_id& index, int analysis_param){
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return ((get_confidence_interval(index,analysis_param) / (2.*get_estimate(index,analysis_param))) < pattern_error_limit) &&
          (pattern_list[index.val_index].num_schedules >= pattern_count_limit) &&
          (pattern_list[index.val_index].total_exec_time >= pattern_time_limit);
}

bool steady_test(const comm_pattern_key& key, const pattern& p, int analysis_param){
  if (!is_key_skipable(key)) return false;
  return is_steady(p,analysis_param);
}
bool steady_test(const comm_pattern_key& key, const pattern_key_id& index, int analysis_param){
  if (!is_key_skipable(key)) return false;
  return is_steady(index,analysis_param);
}
bool steady_test(const comp_pattern_key& key, const pattern& p, int analysis_param){
  if (!is_key_skipable(key)) return false;
  return is_steady(p,analysis_param);
}
bool steady_test(const comp_pattern_key& key, const pattern_key_id& index, int analysis_param){
  if (!is_key_skipable(key)) return false;
  return is_steady(index,analysis_param);
}

void update_kernel_stats(pattern& p, int analysis_param, volatile double exec_time, double unit_count){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  if (exec_time == 0) { exec_time=1.e-9; }
  if (p.steady_state == 0){
    p.num_schedules++;
    p.num_scheduled_units += unit_count;
    p.num_propagations++;
    p.total_exec_time += exec_time;
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
    p.num_non_propagations++;
  }
}
void update_kernel_stats(const pattern_key_id& index, int analysis_param, volatile double exec_time, double unit_count){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  if (exec_time == 0) { exec_time=1.e-9; }
  if (pattern_list[index.val_index].steady_state == 0){
    pattern_list[index.val_index].num_schedules++;
    pattern_list[index.val_index].num_scheduled_units += unit_count;
    pattern_list[index.val_index].num_propagations++;
    pattern_list[index.val_index].total_exec_time += exec_time;
    // Online computation of up to 4th-order central moments using compunication time samples
    size_t n1 = pattern_list[index.val_index].num_schedules-1;
    size_t n = pattern_list[index.val_index].num_schedules;
    double x;
    if (analysis_param == 0){x = exec_time; }	// prep for arithmetic mean
    else                   {x = (unit_count>0 ? unit_count : 1.)/exec_time; }	// prep for harmonic mean
    double delta = x - pattern_list[index.val_index].M1;
    double delta_n = delta / n;
    double delta_n2 = delta_n*delta_n;
    double term1 = delta*delta_n*n1;
    pattern_list[index.val_index].M1 += delta_n;
    pattern_list[index.val_index].M2 += term1;
  }
  else{
    pattern_list[index.val_index].num_non_schedules++;
    pattern_list[index.val_index].num_non_scheduled_units += unit_count;
    pattern_list[index.val_index].num_non_propagations++;
  }
}
void update_kernel_stats(pattern& dest, const pattern& src, int analysis_param){
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
  dest.num_non_schedules += src.num_non_schedules;
  dest.num_non_scheduled_units += src.num_non_scheduled_units;
  dest.num_propagations += src.num_propagations;
  dest.num_non_propagations += src.num_non_propagations;
  dest.total_exec_time += src.total_exec_time;
}

void update_kernel_stats(idle_pattern& p, bool is_global_steady_state, volatile double exec_time){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  if (!is_global_steady_state){
    p.num_schedules++;
    // Online computation of up to 4th-order central moments using compunication time samples
    size_t n1 = p.num_schedules-1;
    size_t n = p.num_schedules;
    double x = exec_time;
    double delta = x - p.M1;
    double delta_n = delta / n;
    double delta_n2 = delta_n*delta_n;
    double term1 = delta*delta_n*n1;
    p.M1 += delta_n;
    p.M2 += term1;
  }
  else{
    p.num_non_schedules++;
  }
}

int should_schedule(const pattern& p){
  return (p.steady_state==1 ? 0 : 1);
}
int should_schedule(const pattern_key_id& index){
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return (pattern_list[index.val_index].steady_state==1 ? 0 : 1);
}

int should_schedule_global(const pattern& p){
  return (p.global_steady_state==1 ? 0 : 1);
}
int should_schedule_global(const pattern_key_id& index){
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return (pattern_list[index.val_index].global_steady_state==1 ? 0 : 1);
}

void set_kernel_state(pattern& p, bool schedule_decision){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  p.steady_state = (schedule_decision==true ? 0 : 1);
}
void set_kernel_state(const pattern_key_id& index, bool schedule_decision){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  pattern_list[index.val_index].steady_state = (schedule_decision==true ? 0 : 1);
}

void set_kernel_state_global(pattern& p, bool schedule_decision){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  p.global_steady_state = (schedule_decision==true ? 0 : 1);
}
void set_kernel_state_global(const pattern_key_id& index, bool schedule_decision){
  if (update_analysis == 0) return;// no updating of analysis -- useful when leveraging data post-autotuning phase
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  pattern_list[index.val_index].global_steady_state = (schedule_decision==true ? 0 : 1);
}


void allocate(MPI_Comm comm){
  int _world_size; MPI_Comm_size(MPI_COMM_WORLD,&_world_size);
  int _world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&_world_rank);
  mode_1_width = 25;
  mode_2_width = 15;
  comm_channel_tag_count=INT_MIN;// to avoid conflict with p2p, which could range from (-p,p)

  comm_channel_node* world_node = new comm_channel_node();
  world_node->tag = comm_channel_tag_count++;
  world_node->offset.push_back(0);
  world_node->id.push_back(std::make_pair(_world_size,1));
  world_node->parent=nullptr;
  sample_propagation_tree* tree = new sample_propagation_tree;
  tree->root = world_node;
  spf.tree = tree;
  comm_channel_map[MPI_COMM_WORLD] = world_node;

  comp_pattern_key ex_1;
  MPI_Datatype comp_pattern_key_internal_type[2] = { MPI_INT, MPI_DOUBLE };
  int comp_pattern_key_internal_type_block_len[2] = { 7,1 };
  MPI_Aint comp_pattern_key_internal_type_disp[2] = { (char*)&ex_1.tag-(char*)&ex_1, (char*)&ex_1.flops-(char*)&ex_1 };
  PMPI_Type_create_struct(2,comp_pattern_key_internal_type_block_len,comp_pattern_key_internal_type_disp,comp_pattern_key_internal_type,&comp_pattern_key_type);
  PMPI_Type_commit(&comp_pattern_key_type);

  comm_pattern_key ex_2;
  MPI_Datatype comm_pattern_key_internal_type[2] = { MPI_INT, MPI_DOUBLE };
  int comm_pattern_key_internal_type_block_len[2] = { 9,1 };
  MPI_Aint comm_pattern_key_internal_type_disp[2] = { (char*)&ex_2.tag-(char*)&ex_2, (char*)&ex_2.msg_size-(char*)&ex_2 };
  PMPI_Type_create_struct(2,comm_pattern_key_internal_type_block_len,comm_pattern_key_internal_type_disp,comm_pattern_key_internal_type,&comm_pattern_key_type);
  PMPI_Type_commit(&comm_pattern_key_type);

  pattern ex_3;
  MPI_Datatype pattern_internal_type[2] = { MPI_INT, MPI_DOUBLE };
  int pattern_internal_block_len[2] = { 6,5 };
  MPI_Aint pattern_internal_disp[2] = { (char*)&ex_3.steady_state-(char*)&ex_3, (char*)&ex_3.num_scheduled_units-(char*)&ex_3 };
  PMPI_Type_create_struct(2,pattern_internal_block_len,pattern_internal_disp,pattern_internal_type,&pattern_type);
  PMPI_Type_commit(&pattern_type);

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
  critical_path_costs[num_critical_path_measures-1]+=(last_time-computation_timer);	// update critical path runtime
  volume_costs[num_volume_measures-1]+=(last_time-computation_timer);			// update runtime volume

  PMPI_Allreduce(MPI_IN_PLACE,&critical_path_costs[0],critical_path_costs.size(),MPI_DOUBLE,MPI_MAX,comm);
  PMPI_Allreduce(MPI_IN_PLACE,&volume_costs[0],volume_costs.size(),MPI_DOUBLE,MPI_MAX,comm);
  // Find the max per-process overhead
  double intercept_overhead[3] = {comm_intercept_overhead_stage1,comm_intercept_overhead_stage2,
                                  comp_intercept_overhead};
  PMPI_Allreduce(MPI_IN_PLACE,&intercept_overhead,3,MPI_DOUBLE,MPI_MAX,comm);
  comm_intercept_overhead_stage1 = intercept_overhead[0];
  comm_intercept_overhead_stage2 = intercept_overhead[1];
  comp_intercept_overhead = intercept_overhead[2];
}

void reset(bool track_statistical_data_override, bool schedule_kernels_override, bool force_steady_statistical_data_overide, bool update_statistical_data_overide){
  for (auto i=0; i<list_size; i++){ list[i]->init(); }
  memset(&critical_path_costs[0],0,sizeof(double)*critical_path_costs.size());
  memset(&max_per_process_costs[0],0,sizeof(double)*max_per_process_costs.size());
  memset(&volume_costs[0],0,sizeof(double)*volume_costs.size());

  // Reset these global variables, as some are updated by function arguments for convenience
  if (std::getenv("CRITTER_AUTOTUNING_MODE") != NULL){
    analysis_mode = atoi(std::getenv("CRITTER_AUTOTUNING_MODE"));
    assert(analysis_mode>0 && analysis_mode<=3);
  } else{
    analysis_mode = 0;
  }
  if (std::getenv("CRITTER_AUTOTUNING_OPTIMIZE") != NULL){
    is_optimized = atoi(std::getenv("CRITTER_AUTOTUNING_OPTIMIZE"));
    assert(is_optimized>=0 && is_optimized<=1);
  } else{
    is_optimized = 0;
  }
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
    pattern_count_limit = atoi(std::getenv("CRITTER_PATTERN_COUNT_LIMIT"));
  } else{
    pattern_count_limit = 1;
  }
  if (std::getenv("CRITTER_PATTERN_TIME_LIMIT") != NULL){
    pattern_time_limit = atof(std::getenv("CRITTER_PATTERN_TIME_LIMIT"));
  } else{
    pattern_time_limit = .00001;
  }
  if (std::getenv("CRITTER_PATTERN_ERROR_LIMIT") != NULL){
    pattern_error_limit = atof(std::getenv("CRITTER_PATTERN_ERROR_LIMIT"));
  } else{
    pattern_error_limit = .5;
  }
  if (std::getenv("CRITTER_SCHEDULE_KERNELS") != NULL){
    schedule_kernels = atoi(std::getenv("CRITTER_SCHEDULE_KERNELS"));
  } else{
    schedule_kernels = 1;
  }
  if (std::getenv("CRITTER_COMM_SAMPLE_INCLUDE_IDLE") != NULL){
    comm_sample_include_idle = atoi(std::getenv("CRITTER_COMM_SAMPLE_INCLUDE_IDLE"));
  } else{
    comm_sample_include_idle = 0;
  }
  if (analysis_mode>0){ analysis_mode = (track_statistical_data_override ? analysis_mode : 0); }
  if (schedule_kernels==1){ schedule_kernels = (schedule_kernels_override ? schedule_kernels : 0); }
  update_analysis = (update_statistical_data_overide ? 1 : 0);
  autotuning_propagate=1;// means nothing if analysis_mode != 3
  if (force_steady_statistical_data_overide && analysis_mode<3){// DO NOT SET TO STEADY_STATE WHEN PERFORMING CRITICAL PATH ANALYSIS
    // This branch is to be entered only after tuning a space of algorithmic parameterizations, in which the expectation is that all kernels,
    //   both comm and comp, have reached a sufficiently-predictable state (steady state). It is also only valid for per-process or volumetric analysis.
/*
    .. my problem with this now is that it breaks the invariant that kernels in global steady state must be both inactive and in the steady state buffers.
    ..   I have recently modified the routines in path.cxx to expect this, which I did not do before (and thus is the reason the code below worked for pp/vol analysis).
    .. I cannot simply flush each pattern that is not aleady in global steady state, as 
    // set all kernels into global steady state -- note this is reasonable for now,
    //   particularly as a debugging technique, but not sure if this makes sense permanently
    for (auto it : comm_pattern_map){
      set_kernel_state(it.second,false);
      set_kernel_state_global(it.second,false);
    }
    for (auto it : comp_pattern_map){
      set_kernel_state(it.second,false);
      set_kernel_schedule_global(it.second,false);
    }
*/
  }
  else{
    update_analysis=1;// Its too risky to allow no updating to pattern members when performing critical path analysis.
/*
    .. yes but how do we protect the statistics generated when tuning the parameterization space from being corrupted, even if most kernels are already steady?
    ..   because I think they would still update the "non" pattern members, right?
*/
  }
  previous_comm_key = comm_pattern_key();
}

void reset_frequencies(){
  spf.clear_info();
}

void clear(){
  // I don't see any reason to clear the communicator map. In fact, doing so would be harmful
  comm_pattern_map.clear();
  comp_pattern_map.clear();
  steady_state_comm_pattern_keys.clear();
  active_comm_pattern_keys.clear();
  steady_state_comp_pattern_keys.clear();
  active_comp_pattern_keys.clear();
  steady_state_patterns.clear();
  active_patterns.clear();
  comm_pattern_pair_map.clear();
  comm_intercept_overhead_stage1=0;
  comm_intercept_overhead_stage2=0;
  comp_intercept_overhead=0;
  previous_comm_key = comm_pattern_key();

  spf.clear_info();
}

void finalize(){
  // 'spf' deletion should occur automatically at program shutdown
  //TODO: Fix memory leak of allocated comm_channel_nodes in comm_channel_map and p2p_channel_map.
  for (auto it : intermediate_channels){
    delete it;
  }
}

}
}
}
