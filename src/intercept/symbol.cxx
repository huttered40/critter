#include "symbol.h"
#include "../util/util.h"
#include "../dispatch/dispatch.h"

namespace critter{
namespace internal{

void symbol_start(const char* symbol){
  if (mode && symbol_path_select_size>0){
    volatile double save_time = MPI_Wtime();
    open_symbol(symbol,save_time);
  }
}

void symbol_stop(const char* symbol){
  if (mode && symbol_path_select_size>0){
    volatile double save_time = MPI_Wtime();
    close_symbol(symbol,save_time);
  }
}

}
}
