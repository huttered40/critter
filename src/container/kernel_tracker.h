#ifndef CRITTER__CONTAINER__KERNEL_TRACKER_H_
#define CRITTER__CONTAINER__KERNEL_TRACKER_H_

#include <string>
#include <vector>
#include <unordered_map>
#include <stack>

namespace internal{

// One instance for each unique kernel
class kernel_tracker{
  public:
    kernel_tracker() {}
    kernel_tracker(const std::string& name_);
    void stop(double save_time);
    void start(double save_time);
    bool operator<(const kernel_tracker& w) const ;

    std::string name;
    std::stack<double> start_timer;
    std::vector<float*> cp_exclusive_contributions;
    std::vector<float*> cp_exclusive_measure;
    std::vector<float*> cp_numcalls;
    std::vector<float*> cp_incl_measure;
    std::vector<float*> cp_excl_measure;
    float* pp_exclusive_contributions;
    float* pp_exclusive_measure;
    float* pp_numcalls,*vol_numcalls;
    float* pp_incl_measure,*pp_excl_measure;
    float* vol_incl_measure,*vol_excl_measure;
};

extern std::unordered_map<std::string,kernel_tracker> timers;

}

#endif /*CRITTER__CONTAINER__KERNEL_TRACKER_H_*/
