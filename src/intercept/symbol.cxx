#include "symbol.h"
#include "../util.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{

void symbol_start(const char* symbol){
  if (mode && symbol_path_select_size>0){
    volatile double save_time = MPI_Wtime();
    if (symbol_timers.find(symbol) == symbol_timers.end()){
      symbol_timers[symbol] = symbol_tracker(symbol);
      symbol_order[symbol_timers.size()-1] = symbol;
      symbol_timers[symbol].start(save_time);
    }
    else{
      symbol_timers[symbol].start(save_time);
    }
  }
}

void symbol_stop(const char* symbol){
  if (mode && symbol_path_select_size>0){
    volatile double save_time = MPI_Wtime();
    if (symbol_timers.find(symbol) == symbol_timers.end()){
      assert(0);
    }
    else{
      symbol_timers[symbol].stop(save_time);
    }
  }
}

}
}
