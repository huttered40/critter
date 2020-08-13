#ifndef CRITTER__DECOMPOSITION__UTIL__UTIL_H_
#define CRITTER__DECOMPOSITION__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace decomposition{

extern size_t cp_symbol_class_count;
extern size_t pp_symbol_class_count;
extern size_t vol_symbol_class_count;
extern size_t max_num_symbols;
extern size_t max_timer_name_length;
extern std::string _cost_models_,_symbol_path_select_,_comm_path_select_;
extern size_t cost_model_size;
extern size_t symbol_path_select_size;
extern size_t comm_path_select_size;
extern std::vector<char> cost_models;
extern std::vector<char> symbol_path_select;
extern std::vector<char> comm_path_select;
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
extern std::vector<int> symbol_path_select_index;
extern std::vector<event> event_list;
extern std::vector<int> opt_req_match;
extern std::vector<double> opt_measure_match;
extern size_t event_list_size;
extern size_t opt_max_iter;
extern size_t gradient_jump_size;
extern size_t num_gradient_points;

void allocate(MPI_Comm comm);
void reset();
void open_symbol(const char* symbol, double curtime);
void close_symbol(const char* symbol, double curtime);
void final_accumulate(MPI_Comm comm, double last_time);
void clear();

}
}
}

#endif /*CRITTER__DECOMPOSITION__UTIL__UTIL_H_*/
