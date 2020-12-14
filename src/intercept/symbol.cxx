#include "symbol.h"
#include "../util/util.h"
#include "../dispatch/dispatch.h"

namespace critter{

void symbol_invoke(const char* symbol, float flops, int param1, int param2, int param3, int param4, int param5){
  if (internal::mode){
    volatile float curtime = MPI_Wtime();
    std::string _symbol_ = symbol;
    assert(internal::symbol_id_map.find(_symbol_) != internal::symbol_id_map.end());
    bool schedule_decision = internal::initiate_comp(internal::symbol_id_map[_symbol_],curtime,flops,param1,param2,param3,param4,param5);
    if (schedule_decision){
      internal::symbol_function();
    }
    internal::complete_comp(0,internal::symbol_id_map[_symbol_],flops,param1,param2,param3,param4,param5);
  }
  else{
    internal::symbol_function();
  }
}

namespace internal{

void symbol_start(const char* symbol){
  if (mode){
    volatile float save_time = MPI_Wtime();
    open_symbol(symbol,save_time);
  }
}

void symbol_stop(const char* symbol){
  if (mode){
    volatile float save_time = MPI_Wtime();
    close_symbol(symbol,save_time);
  }
}

void symbol_conditional_start(const char* symbol, std::function<void(void)> f){
  std::string _symbol_ = symbol;
  if (internal::symbol_id_map.find(_symbol_) == internal::symbol_id_map.end()){
    internal::symbol_id_map[_symbol_] = internal::symbol_id_count++;
  }
  internal::symbol_function = f;
}

}
}
