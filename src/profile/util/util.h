#ifndef CRITTER__PROFILE__UTIL__UTIL_H_
#define CRITTER__PROFILE__UTIL__UTIL_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace profile{

extern std::ofstream stream;
extern int cost_model,path_decomposition;
extern int internal_tag,internal_tag1,internal_tag2;
extern int internal_tag3,internal_tag4,internal_tag5;
extern bool is_first_request,is_first_iter;
extern size_t decomp_text_width,text_width,path_count,path_measure_count;
extern double comp_start_time;
// CommCost, SynchCost, CompCost,
// CommTime, ExecTime
extern size_t num_cp_measures;
// CommCost, SynchCost, CompCost,
// CommTime, ExecTime
extern size_t num_pp_measures;
// CommCost, SynchCost, CompCost,
// CommTime, ExecTime
extern size_t num_vol_measures;
extern size_t num_decomp_cp_measures;
extern size_t num_decomp_pp_measures;
extern size_t num_decomp_vol_measures;
extern size_t cp_costs_size,pp_costs_size,vol_costs_size;
extern std::string path_select,path_measure_select;
extern std::vector<bool> path_decisions;
extern std::vector<float> cp_costs,max_pp_costs,vol_costs;
extern std::vector<float> cp_costs_foreign;
extern std::vector<char> eager_pad;
extern std::vector<int> path_index,path_measure_index;
extern std::stack<std::string> symbol_stack;
extern std::vector<std::string> symbol_order;
//extern std::vector<float> intercept_overhead,global_intercept_overhead;
extern std::map<comp_kernel_key,std::pair<int,float>> comp_kernel_info;
extern std::map<comm_kernel_key,std::pair<int,float>> comm_kernel_info;

void allocate(MPI_Comm comm);
void reset();
void init_symbol(std::vector<std::string>& symbols);
void open_symbol(const char* symbol, double curtime);
void close_symbol(const char* symbol, double curtime);
void final_accumulate(MPI_Comm comm, double last_time);
void clear();
void finalize();

}
}
}

#endif /*CRITTER__PROFILE__UTIL__UTIL_H_*/
