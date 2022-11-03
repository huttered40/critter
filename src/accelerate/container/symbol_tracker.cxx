#include "symbol_tracker.h"

namespace critter{
namespace internal{
namespace accelerate{

// Global namespace variable 'symbol_timers' must be defined here, rather than in src/util.cxx with the rest, to avoid circular dependence between this file and src/util.h
std::unordered_map<std::string,symbol_tracker> symbol_timers;

symbol_tracker::symbol_tracker(std::string name_){}

void symbol_tracker::start(float save_time){}

void symbol_tracker::stop(float save_time){}

}
}
}
