#ifndef CRITTER__UTIL__UTIL_H_
#define CRITTER__UTIL__UTIL_H_

#include <mpi.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <vector>
#include <stack>
#include <stdint.h>
#include <functional>
#include <map>
#include <set>
#include <unordered_map>
#include <cmath>
#include <assert.h>

namespace critter{
namespace internal{

// ****************************************************************************************************************************************************
struct double_int{
  double_int(){first=0; second=0;}
  double_int(double one, int two){first=one; second=two;}
  double first; int second;
};
struct int_int_double{
  int_int_double(){first=0; second=0; third=0;}
  int_int_double(int one, int two, double three){first=one; second=two; third=three;}
  int first; int second; double third;
};

// ****************************************************************************************************************************************************
struct event{
  event(std::string _kernel, std::vector<double> _measurements){
    tag = -1; measurements = _measurements;
    kernel = _kernel;
  }
  event(std::string _kernel, std::vector<double> _measurements, int _tag, MPI_Comm _comm, int _partner1, int _partner2, bool _is_sender, bool _is_eager){
    measurements = _measurements;
    tag = _tag; comm = _comm; partner1 = _partner1; partner2 = _partner2; is_sender = _is_sender; is_eager = _is_eager;
    kernel = _kernel;
  }
  event(std::string _kernel, std::vector<double> _measurements, int _tag, MPI_Comm _comm, int _partner1, bool _is_sender, bool _is_eager, int _match_id, bool blah){
    measurements = _measurements;
    tag = _tag; comm = _comm; partner1 = _partner1; is_sender = _is_sender; is_eager = _is_eager; match_id = _match_id; is_close = false;;
    kernel = _kernel;
  }
  event(std::string _kernel,std::vector<double> _measurements, std::vector<int> _match_vec){
    measurements = _measurements;
    tag=18; match_vec = _match_vec; match_size = _match_vec.size(); is_close = true;
    kernel = _kernel;
  }
  MPI_Comm comm;
  int tag,partner1,partner2,match_size,match_id;
  bool is_sender,is_eager,is_close;
  std::vector<int> match_vec;
  std::vector<double> measurements;
  std::string kernel;
};

// ****************************************************************************************************************************************************
struct comm_pattern_param1_key{
  friend bool operator==(const comm_pattern_param1_key& ref1, const comm_pattern_param1_key& ref2){
    if ((ref1.tag==ref2.tag) && (ref1.comm_size == ref2.comm_size) && (ref1.comm_color == ref2.comm_color) && (ref1.msg_size == ref2.msg_size) && (ref1.partner_offset == ref2.partner_offset)) return true;
    else return false;
  }
  friend bool operator<(const comm_pattern_param1_key& ref1, const comm_pattern_param1_key& ref2){
    if (ref1.tag < ref2.tag) return true;
    else if (ref1.tag > ref2.tag) return false;
    if (ref1.comm_size < ref2.comm_size) return true;
    else if (ref1.comm_size > ref2.comm_size) return false;
    if (ref1.comm_color < ref2.comm_color) return true;
    else if (ref1.comm_color > ref2.comm_color) return false;
    if (ref1.msg_size < ref2.msg_size) return true;
    else if (ref1.msg_size > ref2.msg_size) return false;
    if (ref1.partner_offset < ref2.partner_offset) return true;
    else if (ref1.partner_offset > ref2.partner_offset) return false;
    return false;
  }
  int tag;
  int comm_color;
  int comm_size;
  double msg_size;
  int partner_offset;
};

// ****************************************************************************************************************************************************
struct comm_pattern_param1_val{
  // If I leverage the kurtosis, I will have to utilize the arithmetic mean.
  //   Note that I'd rather utilize the geometric mean, but I'm not sure how to convert this algorithm
  //     to handle that.
  comm_pattern_param1_val(size_t count_limit=0, double error_limit=0., double time_limit=0.){
    this->count_limit = count_limit;
    this->error_limit = error_limit;
    this->time_limit = time_limit;
    this->steady_state=false;
    this->num_schedules = 0;
    this->num_non_schedules = 0;
    this->num_scheduled_bytes = 0;
    this->num_non_scheduled_bytes = 0;
    this->M1=0; this->M2=0; this->M3=0; this->M4=0;
  }
  double get_confidence_interval(double level = .95){
    // returns confidence interval length with 95% confidence level
    return 1.96*this->get_std_error();
  }
  double get_std_error(){
    // returns standard error
    size_t n = this->num_schedules;
    return this->get_std_dev() / pow(n*1.,1./2.);
  }
  double get_arithmetic_mean(){
    // returns arithmetic mean
    return this->M1;
  }
  double get_std_dev(){
    // returns variance
    return pow(this->get_variance(),1./2.);
  }
  double get_variance(){
    // returns variance
    size_t n = this->num_schedules;
    if (n<=1) return 1000.;
    return this->M2 / (n-1.);
  }
  double get_skewness(){
    // returns skewness
    size_t n = this->num_schedules;
    if (n<=2) return 1000.;
    return std::pow(n*1.,1./2.) * this->M3 / std::pow(this->M2,1.5);
  }
  double get_kurtosis(){
    // returns the excess kurtosis
    size_t n = this->num_schedules;
    if (n<=3) return 1000.;
    return (n*this->M4) / (this->M2*this->M2) - 3.;
  }
  double get_jacque_barra(){
    // Note this test may be too sensitive for our purposes. But we can try this first.
    size_t n = this->num_schedules;
    double k = this->get_kurtosis();
    double s = this->get_skewness();
    double jb_test_stat = (n/6.)*(s*s + (1./4.)*(k*k));
    return jb_test_stat;
  }
  void error_test(){
    bool decision = ((get_confidence_interval() / (2.*get_arithmetic_mean())) < this->error_limit) || (this->num_schedules > this->count_limit) || (this->total_comm_time > this->time_limit);
    if (decision){
      this->steady_state = true;
    } else{
      this->steady_state = false;
    }
  }
  bool should_schedule(){
    return !this->steady_state;
  }
  void update(volatile double comm_time, double byte_count){
    if (!this->steady_state){
      this->num_schedules++;
      this->num_scheduled_bytes += byte_count;
      this->total_comm_time += comm_time;
      save_comm_times.push_back((double)comm_time);	// only temporary
      // Online computation of up to 4th-order central moments using communication time samples
      size_t n1 = this->num_schedules-1;
      size_t n = this->num_schedules;
      double x = comm_time;
      double delta = x - M1;
      double delta_n = delta / n;
      double delta_n2 = delta_n*delta_n;
      double term1 = delta*delta_n*n1;
      this->M1 += delta_n;
      this->M4 = this->M4 + term1*delta_n2*(n*n - 3*n + 3) + 6*delta_n2*this->M2-4*delta_n*this->M3;
      this->M3 = this->M3 + term1*delta_n*(n-2)-3*delta_n*this->M2;
      this->M2 += term1;
      error_test();
      this->save_arithmetic_means.push_back(get_arithmetic_mean());
      this->save_std_dev.push_back(get_std_dev());
      this->save_std_error.push_back(get_std_error());
      this->save_confidence_interval.push_back(get_confidence_interval());
      //this->save_skewness.push_back(get_skewness());
      //this->save_kurtosis.push_back(get_kurtosis());
      //this->save_jb.push_back(get_jacque_barra());
    }
    else{
      this->num_non_schedules++;
      this->num_non_scheduled_bytes += byte_count;
    }
  }

  size_t count_limit;
  double error_limit;
  double time_limit;
  size_t num_schedules;
  size_t num_non_schedules;
  double num_scheduled_bytes;
  double num_non_scheduled_bytes;
  double M1,M2,M3,M4;
  bool steady_state;
  double total_comm_time;
  std::vector<double> save_comm_times;
  std::vector<double> save_arithmetic_means;
  std::vector<double> save_std_dev;
  std::vector<double> save_std_error;
  std::vector<double> save_confidence_interval;
  //std::vector<double> save_skewness;
  //std::vector<double> save_kurtosis;
  //std::vector<double> save_jb;
};

// ****************************************************************************************************************************************************
struct comp_pattern_param1_key{
  friend bool operator==(const comp_pattern_param1_key& ref1, const comp_pattern_param1_key& ref2){
    if ((ref1.tag==ref2.tag) && (ref1.param1 == ref2.param1) && (ref1.param2 == ref2.param2) && (ref1.param3 == ref2.param3) && (ref1.param4 == ref2.param4) && (ref1.param5 == ref2.param5)) return true;
    else return false;
  }
  friend bool operator<(const comp_pattern_param1_key& ref1, const comp_pattern_param1_key& ref2){
    if (ref1.tag < ref2.tag) return true;
    else if (ref1.tag > ref2.tag) return false;
    if (ref1.param1 < ref2.param1) return true;
    else if (ref1.param1 > ref2.param1) return false;
    if (ref1.param2 < ref2.param2) return true;
    else if (ref1.param2 > ref2.param2) return false;
    if (ref1.param3 < ref2.param3) return true;
    else if (ref1.param3 > ref2.param3) return false;
    if (ref1.param4 < ref2.param4) return true;
    else if (ref1.param4 > ref2.param4) return false;
    if (ref1.param5 < ref2.param5) return true;
    else if (ref1.param5 > ref2.param5) return false;
    return false;
  }
  int tag;
  double flops;
  int param1,param2,param3,param4,param5;
};

// ****************************************************************************************************************************************************
struct comp_pattern_param1_val{
  // If I leverage the kurtosis, I will have to utilize the arithmetic mean.
  //   Note that I'd rather utilize the geometric mean, but I'm not sure how to convert this algorithm
  //     to handle that.
  comp_pattern_param1_val(size_t count_limit=0, double error_limit=0., double time_limit=0.){
    this->count_limit = count_limit;
    this->error_limit = error_limit;
    this->time_limit = time_limit;
    this->steady_state=false;
    this->num_schedules = 0;
    this->num_non_schedules = 0;
    this->num_scheduled_flops = 0;
    this->num_non_scheduled_flops = 0;
    this->M1=0; this->M2=0; this->M3=0; this->M4=0;
  }
  double get_confidence_interval(double level = .95){
    // returns confidence interval length with 95% confidence level
    return 1.96*this->get_std_error();
  }
  double get_std_error(){
    // returns standard error
    size_t n = this->num_schedules;
    return this->get_std_dev() / pow(n*1.,1./2.);
  }
  double get_arithmetic_mean(){
    // returns arithmetic mean
    return this->M1;
  }
  double get_std_dev(){
    // returns variance
    return pow(this->get_variance(),1./2.);
  }
  double get_variance(){
    // returns variance
    size_t n = this->num_schedules;
    if (n<=1) return 1000.;
    return this->M2 / (n-1.);
  }
  double get_skewness(){
    // returns skewness
    size_t n = this->num_schedules;
    if (n<=2) return 1000.;
    return std::pow(n*1.,1./2.) * this->M3 / std::pow(this->M2,1.5);
  }
  double get_kurtosis(){
    // returns the excess kurtosis
    size_t n = this->num_schedules;
    if (n<=3) return 1000.;
    return (n*this->M4) / (this->M2*this->M2) - 3.;
  }
  double get_jacque_barra(){
    // Note this test may be too sensitive for our purposes. But we can try this first.
    size_t n = this->num_schedules;
    double k = this->get_kurtosis();
    double s = this->get_skewness();
    double jb_test_stat = (n/6.)*(s*s + (1./4.)*(k*k));
    return jb_test_stat;
  }
  void error_test(){
    bool decision = ((get_confidence_interval() / (2.*get_arithmetic_mean())) < this->error_limit) || (this->num_schedules > this->count_limit) || (this->total_comp_time > this->time_limit);
    if (decision){
      this->steady_state = true;
    } else{
      this->steady_state = false;
    }
  }
  bool should_schedule(){
    return !this->steady_state;
  }
  void update(volatile double comp_time, double flop_count){
    if (!this->steady_state){
      this->num_schedules++;
      this->num_scheduled_flops += flop_count;
      this->total_comp_time += comp_time;
      save_comp_times.push_back((double)comp_time);	// only temporary
      // Online computation of up to 4th-order central moments using compunication time samples
      size_t n1 = this->num_schedules-1;
      size_t n = this->num_schedules;
      double x = comp_time;
      double delta = x - M1;
      double delta_n = delta / n;
      double delta_n2 = delta_n*delta_n;
      double term1 = delta*delta_n*n1;
      this->M1 += delta_n;
      this->M4 = this->M4 + term1*delta_n2*(n*n - 3*n + 3) + 6*delta_n2*this->M2-4*delta_n*this->M3;
      this->M3 = this->M3 + term1*delta_n*(n-2)-3*delta_n*this->M2;
      this->M2 += term1;
      error_test();
      this->save_arithmetic_means.push_back(get_arithmetic_mean());
      this->save_std_dev.push_back(get_std_dev());
      this->save_std_error.push_back(get_std_error());
      this->save_confidence_interval.push_back(get_confidence_interval());
      //this->save_skewness.push_back(get_skewness());
      //this->save_kurtosis.push_back(get_kurtosis());
      //this->save_jb.push_back(get_jacque_barra());
    }
    else{
      this->num_non_schedules++;
      this->num_non_scheduled_flops += flop_count;
    }
  }

  size_t count_limit;
  double error_limit;
  double time_limit;
  size_t num_schedules;
  size_t num_non_schedules;
  double num_scheduled_flops;
  double num_non_scheduled_flops;
  double M1,M2,M3,M4;
  bool steady_state;
  double total_comp_time;
  std::vector<double> save_comp_times;
  std::vector<double> save_arithmetic_means;
  std::vector<double> save_std_dev;
  std::vector<double> save_std_error;
  std::vector<double> save_confidence_interval;
  //std::vector<double> save_skewness;
  //std::vector<double> save_kurtosis;
  //std::vector<double> save_jb;
};

// ****************************************************************************************************************************************************
extern size_t path_pattern_param;
extern size_t path_pattern_comm_count_limit;
extern size_t path_pattern_comp_count_limit;
extern double path_pattern_comm_time_limit;
extern double path_pattern_comp_time_limit;
extern double path_pattern_comm_error_limit;
extern double path_pattern_comp_error_limit;
extern std::map<MPI_Comm,std::pair<int,int>> communicator_map;
extern std::vector<int> local_communicator_list;
extern std::map<comm_pattern_param1_key,comm_pattern_param1_val> comm_pattern_cache_param1;
extern std::map<comp_pattern_param1_key,comp_pattern_param1_val> comp_pattern_cache_param1;
extern double comp_start_time;
extern size_t cp_symbol_class_count;
extern size_t pp_symbol_class_count;
extern size_t vol_symbol_class_count;
extern size_t mode_1_width;
extern size_t mode_2_width;
extern size_t max_num_symbols;
extern size_t max_timer_name_length;
extern std::string _cost_models_,_symbol_path_select_,_comm_path_select_;
extern size_t cost_model_size;
extern size_t symbol_path_select_size;
extern size_t comm_path_select_size;
extern size_t auto_capture;
extern std::vector<char> cost_models;
extern std::vector<char> symbol_path_select;
extern std::vector<char> comm_path_select;
extern size_t num_critical_path_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime, CompTime, RunTime
extern size_t num_per_process_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, RunTime
extern size_t num_volume_measures;			// CommCost*, SynchCost*, IdleTime, CommTime, SynchTime, CompTime, RunTime
extern size_t num_tracker_critical_path_measures;	// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t num_tracker_per_process_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t num_tracker_volume_measures;		// CommCost*, SynchCost*,           CommTime, SynchTime
extern size_t critical_path_costs_size;
extern size_t per_process_costs_size;
extern size_t volume_costs_size;
extern std::string stream_name,file_name;
extern bool flag,is_first_iter,is_world_root,need_new_line,opt;
extern size_t mechanism,mode,stack_id;
extern std::ofstream stream;
extern volatile double computation_timer;
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
extern std::vector<double> critical_path_costs;
extern std::vector<double> max_per_process_costs;
extern std::vector<double> volume_costs;
extern std::map<std::string,std::vector<double>> save_info;
extern std::vector<double> new_cs;
extern double scratch_pad;
extern std::vector<char> synch_pad_send;
extern std::vector<char> synch_pad_recv;
extern std::vector<char> barrier_pad_send;
extern std::vector<char> barrier_pad_recv;
extern std::vector<char> symbol_pad_cp;
extern std::vector<char> symbol_pad_ncp1;
extern std::vector<char> symbol_pad_ncp2;
extern std::vector<int> symbol_len_pad_cp;
extern std::vector<int> symbol_len_pad_ncp1;
extern std::vector<int> symbol_len_pad_ncp2;
extern std::vector<double> symbol_timer_pad_local_cp;
extern std::vector<double> symbol_timer_pad_global_cp;
extern std::vector<double> symbol_timer_pad_global_cp2;
extern std::vector<double> symbol_timer_pad_local_pp;
extern std::vector<double> symbol_timer_pad_global_pp;
extern std::vector<double> symbol_timer_pad_local_vol;
extern std::vector<double> symbol_timer_pad_global_vol;
extern std::stack<std::string> symbol_stack;
extern std::vector<std::string> symbol_order;
extern std::vector<double_int> info_sender;
extern std::vector<double_int> info_receiver;
extern std::vector<int> symbol_path_select_index;
extern bool wait_id;
extern int internal_tag;
extern int internal_tag1;
extern int internal_tag2;
extern int internal_tag3;
extern int internal_tag4;
extern int internal_tag5;
extern size_t track_blas;
extern size_t track_lapack;
extern size_t track_collective;
extern size_t track_p2p;
extern size_t track_p2p_idle;
extern size_t eager_p2p;
extern size_t delete_comm;
extern std::vector<char> eager_pad;
extern std::vector<event> event_list;
extern std::vector<int> opt_req_match;
extern std::vector<double> opt_measure_match;
extern size_t event_list_size;
extern size_t opt_max_iter;
extern size_t gradient_jump_size;
extern size_t num_gradient_points;
extern size_t
         _MPI_Send__id,
         _MPI_Ssend__id,
         _MPI_Bsend__id,
         _MPI_Recv__id,
         _MPI_Sendrecv__id,
         _MPI_Sendrecv_replace__id,
         _MPI_Barrier__id,
         _MPI_Bcast__id,
         _MPI_Reduce__id,
         _MPI_Allreduce__id,
         _MPI_Gather__id,
         _MPI_Allgather__id,
         _MPI_Scatter__id,
         _MPI_Reduce_scatter__id,
         _MPI_Alltoall__id,
         _MPI_Gatherv__id,
         _MPI_Allgatherv__id,
         _MPI_Scatterv__id,
         _MPI_Alltoallv__id,
         _MPI_Isend__id,
         _MPI_Irecv__id,
         _MPI_Ibcast__id,
         _MPI_Iallreduce__id,
         _MPI_Ireduce__id,
         _MPI_Igather__id,
         _MPI_Igatherv__id,
         _MPI_Iallgather__id,
         _MPI_Iallgatherv__id,
         _MPI_Iscatter__id,
         _MPI_Iscatterv__id,
         _MPI_Ireduce_scatter__id,
         _MPI_Ialltoall__id,
         _MPI_Ialltoallv__id;
extern size_t
	_BLAS_axpy__id,
	_BLAS_scal__id,
	_BLAS_ger__id,
	_BLAS_gemm__id,
	_BLAS_trmm__id,
	_BLAS_trsm__id,
	_BLAS_syrk__id;
extern size_t
	_LAPACK_getrf__id,
	_LAPACK_potrf__id,
	_LAPACK_trtri__id,
	_LAPACK_geqrf__id,
	_LAPACK_orgqr__id,
	_LAPACK_ormqr__id,
	_LAPACK_getri__id,
	_LAPACK_tpqrt__id,
	_LAPACK_tpmqrt__id;
extern std::map<std::pair<std::string,size_t>,bool> schedule_map;

}
}

#endif /*CRITTER__UTIL__UTIL_H_*/
