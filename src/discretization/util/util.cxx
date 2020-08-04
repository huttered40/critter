#include "util.h"

namespace critter{
namespace internal{
namespace discretization{

double get_arithmetic_mean(pattern_key_id index){
  // returns arithmetic mean
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return pattern_list[index.val_index].M1;
}
double get_variance(pattern_key_id index){
  // returns variance
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  size_t n = pattern_list[index.val_index].num_schedules;
  if (n<=1) return 1000.;
  return pattern_list[index.val_index].M2 / (n-1.);
}
double get_std_dev(pattern_key_id index){
  // returns variance
  return pow(get_variance(index),1./2.);
}
double get_std_error(pattern_key_id index){
  // returns standard error
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  size_t n = pattern_list[index.val_index].num_schedules;
  return get_std_dev(index) / pow(n*1.,1./2.);
}
double get_confidence_interval(pattern_key_id index, double level){
  // returns confidence interval length with 95% confidence level
  return 1.96*get_std_error(index);
}
/*
double get_skewness(pattern_key_id index){
  // returns skewness
  size_t n = pattern_list[index].num_schedules;
  if (n<=2) return 1000.;
  return std::pow(n*1.,1./2.) * pattern_list[index].M3 / std::pow(pattern_list[index].M2,1.5);
}
double get_kurtosis(pattern_key_id index){
  // returns the excess kurtosis
  size_t n = pattern_list[index].num_schedules;
  if (n<=3) return 1000.;
  return (n*pattern_list[index].M4) / (pattern_list[index].M2*pattern_list[index].M2) - 3.;
}
double get_jacque_barra(pattern_key_id index){
  // Note this test may be too sensitive for our purposes. But we can try this first.
  size_t n = pattern_list[index].num_schedules;
  double k = get_kurtosis();
  double s = get_skewness();
  double jb_test_stat = (n/6.)*(s*s + (1./4.)*(k*k));
  return jb_test_stat;
}
*/
void error_test(pattern_key_id index){
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  bool decision = ((get_confidence_interval(index) / (2.*get_arithmetic_mean(index))) < pattern_error_limit) &&
                  (pattern_list[index.val_index].num_schedules >= pattern_count_limit) &&
                  (pattern_list[index.val_index].total_exec_time >= pattern_time_limit);
  if (decision){
    pattern_list[index.val_index].steady_state = true;
  } else{
    pattern_list[index.val_index].steady_state = false;
  }
}
bool should_schedule(pattern_key_id index){
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  return !pattern_list[index.val_index].steady_state;
}
void update(pattern_key_id index, volatile double exec_time, double unit_count){
  auto& pattern_list = index.is_active == true ? active_patterns : steady_state_patterns;
  if (!pattern_list[index.val_index].steady_state){
    pattern_list[index.val_index].num_schedules++;
    pattern_list[index.val_index].num_scheduled_units += unit_count;
    pattern_list[index.val_index].total_exec_time += exec_time;
//    pattern_list[index].save_exec_times.push_back((double)exec_time);	// only temporary
    // Online computation of up to 4th-order central moments using compunication time samples
    size_t n1 = pattern_list[index.val_index].num_schedules-1;
    size_t n = pattern_list[index.val_index].num_schedules;
    double x = exec_time;
    double delta = x - pattern_list[index.val_index].M1;
    double delta_n = delta / n;
    double delta_n2 = delta_n*delta_n;
    double term1 = delta*delta_n*n1;
    pattern_list[index.val_index].M1 += delta_n;
//    this->M4 = this->M4 + term1*delta_n2*(n*n - 3*n + 3) + 6*delta_n2*this->M2-4*delta_n*this->M3;
//    this->M3 = this->M3 + term1*delta_n*(n-2)-3*delta_n*this->M2;
    pattern_list[index.val_index].M2 += term1;
    error_test(index);
//    pattern_list[index].save_arithmetic_means.push_back(get_arithmetic_mean());
//    pattern_list[index].save_std_dev.push_back(get_std_dev());
//    pattern_list[index].save_std_error.push_back(get_std_error());
//    pattern_list[index].save_confidence_interval.push_back(get_confidence_interval());
//    pattern_list[index].save_skewness.push_back(get_skewness());
//    pattern_list[index].save_kurtosis.push_back(get_kurtosis());
//    pattern_list[index].save_jb.push_back(get_jacque_barra());
  }
  else{
    pattern_list[index.val_index].num_non_schedules++;
    pattern_list[index.val_index].num_non_scheduled_units += unit_count;
  }
}

}
}
}
