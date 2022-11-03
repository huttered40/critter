#ifndef CRITTER__ACCELERATE__CONTAINER__SYMBOL_TRACKER_H_
#define CRITTER__ACCELERATE__CONTAINER__SYMBOL_TRACKER_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace accelerate{

// One instance for each unique symbol
class symbol_tracker{
  public:
    symbol_tracker() {}
    symbol_tracker(std::string name_);
    bool operator<(const symbol_tracker& w) const ;
    void start(float save_time);
    void stop(float save_time);
};

extern std::unordered_map<std::string,symbol_tracker> symbol_timers;

}
}
}

#endif /*CRITTER__ACCELERATE__CONTAINER__SYMBOL_TRACKER_H_*/
