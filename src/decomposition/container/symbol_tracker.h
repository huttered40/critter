#ifndef CRITTER__DECOMPOSITION__CONTAINER__SYMBOL_TRACKER_H_
#define CRITTER__DECOMPOSITION__CONTAINER__SYMBOL_TRACKER_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace decomposition{

// One instance for each unique symbol
class symbol_tracker{
  public:
    symbol_tracker() {}
    symbol_tracker(std::string name_);
    void stop(double save_time);
    void start(double save_time);
    bool operator<(const symbol_tracker& w) const ;

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

extern std::unordered_map<std::string,symbol_tracker> symbol_timers;

}
}
}

#endif /*CRITTER__DECOMPOSITION__CONTAINER__SYMBOL_TRACKER_H_*/
