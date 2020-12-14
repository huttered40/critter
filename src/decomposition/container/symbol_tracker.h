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
    void stop(float save_time);
    void start(float save_time);
    bool operator<(const symbol_tracker& w) const ;

    std::string name;
    std::stack<float> start_timer;
    std::vector<float*> cp_exclusive_contributions;
    float* pp_exclusive_contributions;
    std::vector<float*> cp_exclusive_measure;
    float* pp_exclusive_measure;
    std::vector<float*> cp_numcalls;
    float* pp_numcalls;
    float* vol_numcalls;
    std::vector<float*> cp_incl_measure;
    std::vector<float*> cp_excl_measure;
    float* pp_incl_measure;
    float* pp_excl_measure;
    float* vol_incl_measure;
    float* vol_excl_measure;
    bool has_been_processed;
};

extern std::unordered_map<std::string,symbol_tracker> symbol_timers;

}
}
}

#endif /*CRITTER__DECOMPOSITION__CONTAINER__SYMBOL_TRACKER_H_*/
