#include "kernel_tracker.h"
#include "../util.h"

namespace internal{

std::unordered_map<std::string,kernel_tracker> timers;

kernel_tracker::kernel_tracker(const std::string& name_){
  if (path_count==0) return;
  this->name = name_;
  this->cp_exclusive_contributions.resize(path_count,nullptr);
  this->cp_exclusive_measure.resize(path_count,nullptr);
  this->cp_numcalls.resize(path_count,nullptr);
  this->cp_incl_measure.resize(path_count,nullptr);
  this->cp_excl_measure.resize(path_count,nullptr);

  size_t cp_path_select_offset = num_kernel_ds*num_decomp_cp_measures+1;
  size_t pp_path_select_offset = num_kernel_ds*num_decomp_pp_measures+1;
  size_t vol_path_select_offset = num_kernel_ds*num_decomp_vol_measures+1;
  size_t cp_offset = num_cp_measures;
  cp_offset += (timers.size()-1)*path_count*cp_path_select_offset;
  size_t pp_offset = num_pp_measures;
  pp_offset += (timers.size()-1)*path_count*pp_path_select_offset;
  size_t vol_offset = num_vol_measures;
  vol_offset += (timers.size()-1)*vol_path_select_offset;
  this->pp_numcalls = &max_pp_costs[pp_offset];
  this->pp_excl_measure = &max_pp_costs[pp_offset+1];
  this->pp_incl_measure = exclusive_only==1 ? &scratch_pad[0] : &max_pp_costs[pp_offset+num_decomp_pp_measures+1];
  this->pp_exclusive_contributions = exclusive_only==1 ? &scratch_pad[0] : &max_pp_costs[pp_offset+2*num_decomp_pp_measures+1];
  this->pp_exclusive_measure = exclusive_only==1 ? &scratch_pad[0] : &max_pp_costs[pp_offset+3*num_decomp_pp_measures+1];

  this->vol_numcalls = &vol_costs[vol_offset];
  this->vol_excl_measure = &vol_costs[vol_offset+1];
  this->vol_incl_measure = exclusive_only==1 ? &scratch_pad[0] : &vol_costs[vol_offset+num_decomp_vol_measures+1];

  for (auto i=0; i<path_count; i++){
    this->cp_numcalls[i] = &cp_costs[cp_offset+cp_path_select_offset*i];
    this->cp_excl_measure[i] = &cp_costs[cp_offset+1+cp_path_select_offset*i];
    this->cp_incl_measure[i] = exclusive_only==1 ? &scratch_pad[0] : &cp_costs[cp_offset+num_decomp_cp_measures+1+cp_path_select_offset*i];
    this->cp_exclusive_contributions[i] = exclusive_only==1 ? &scratch_pad[0] : &cp_costs[cp_offset+2*num_decomp_cp_measures+1+cp_path_select_offset*i];
    this->cp_exclusive_measure[i] = exclusive_only==1 ? &scratch_pad[0] : &cp_costs[cp_offset+3*num_decomp_cp_measures+1+cp_path_select_offset*i];
  }
  // memset not necessary because these members point into global arrays
}

void kernel_tracker::start(double save_time){
  if (path_count==0) return;
  if (timer_stack.size()>0){
    // Add up exclusive time for most-recent symbol (the one about to be supplanted for this new one).
    auto last_timer_time = save_time-timers[timer_stack.top()].start_timer.top();
    for (auto i=0; i<path_count; i++){
      for (auto j=0; j<path_measure_index.size(); j++){
        if (path_measure_index[j] == path_measure_select.size()-1){
          timers[timer_stack.top()].cp_exclusive_measure[i][j] += last_timer_time;
          timers[timer_stack.top()].cp_excl_measure[i][j] += last_timer_time;
        }
      }
    }
    timers[timer_stack.top()].pp_exclusive_measure[num_decomp_pp_measures-1] += last_timer_time;
    timers[timer_stack.top()].pp_excl_measure[num_decomp_pp_measures-1] += last_timer_time;
  }
  auto last_comp_time = save_time - computation_timer;
  cp_costs[num_cp_measures-1] += last_comp_time;
  vol_costs[num_vol_measures-1] += last_comp_time;
  timer_stack.push(this->name);
  computation_timer = MPI_Wtime();
  this->start_timer.push((double)computation_timer);
}

void kernel_tracker::stop(double save_time){
  if (path_count==0) return;
  assert(this->start_timer.size()>0);
  auto last_timer_time = save_time-this->start_timer.top();
  // Add up exclusive time for most-recent symbol (the one about to be removed from top-of-stack).
  this->pp_exclusive_measure[num_decomp_pp_measures-1] += last_timer_time;
  this->pp_excl_measure[num_decomp_pp_measures-1] += last_timer_time;
  *this->pp_numcalls += 1.; *this->vol_numcalls += 1.;
  for (auto i=0; i<num_decomp_pp_measures; i++){
    this->pp_exclusive_contributions[i] += this->pp_exclusive_measure[i];
  }
  for (auto i=0; i<path_count; i++){
    *this->cp_numcalls[i] += 1.;
    for (auto j=0; j<path_measure_index.size(); j++){
      if (path_measure_index[j] == path_measure_select.size()-1){
        this->cp_exclusive_measure[i][j] += last_timer_time;
        this->cp_excl_measure[i][j] += last_timer_time;
      }
      this->cp_exclusive_contributions[i][j] += this->cp_exclusive_measure[i][j];
    }
    memset(this->cp_exclusive_measure[i],0,sizeof(float)*num_decomp_cp_measures);
  }
  memset(this->pp_exclusive_measure,0,sizeof(float)*num_decomp_pp_measures);
  auto save_symbol = timer_stack.top();
  this->start_timer.pop(); timer_stack.pop();

  if (exclusive_only == 0){
    if (timer_stack.size() > 0 && (save_symbol != timer_stack.top())){
      for (auto j=0; j<path_count; j++){
	for (auto i=0; i<num_decomp_cp_measures; i++){
	  timers[timer_stack.top()].cp_exclusive_contributions[j][i] += this->cp_exclusive_contributions[j][i];
	  this->cp_incl_measure[j][i] += this->cp_exclusive_contributions[j][i];
	}
	memset(this->cp_exclusive_contributions[j],0,sizeof(float)*num_decomp_cp_measures);
      }
      for (auto i=0; i<num_decomp_pp_measures; i++){
	timers[timer_stack.top()].pp_exclusive_contributions[i] += this->pp_exclusive_contributions[i];
	this->pp_incl_measure[i] += this->pp_exclusive_contributions[i];
      }
      memset(this->pp_exclusive_contributions,0,sizeof(float)*num_decomp_pp_measures);
    } else if (timer_stack.size() == 0){
      for (auto j=0; j<path_count; j++){
	for (auto i=0; i<num_decomp_cp_measures; i++){
	  this->cp_incl_measure[j][i] += this->cp_exclusive_contributions[j][i];
	}
	memset(this->cp_exclusive_contributions[j],0,sizeof(float)*num_decomp_cp_measures);
      }
      for (auto i=0; i<num_decomp_pp_measures; i++){
	this->pp_incl_measure[i] += this->pp_exclusive_contributions[i];
      }
      memset(this->pp_exclusive_contributions,0,sizeof(float)*num_decomp_pp_measures);
    }
  }
  auto last_comp_time = save_time - computation_timer;
  cp_costs[num_cp_measures-1] += last_comp_time;
  vol_costs[num_vol_measures-1] += last_comp_time;
  computation_timer = MPI_Wtime();
  if (timer_stack.size()>0){ timers[timer_stack.top()].start_timer.top() = computation_timer; }
}

}
