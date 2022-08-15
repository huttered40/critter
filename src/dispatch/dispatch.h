#ifndef CRITTER__DISPATCH__DISPATCH_H_
#define CRITTER__DISPATCH__DISPATCH_H_

#include "../util/util.h"

namespace critter{
namespace internal{

void allocate(MPI_Comm comm);
void reset(bool schedule_kernels_override, bool force_steady_statistical_data_overide);
void exchange_communicators(MPI_Comm oldcomm, MPI_Comm newcomm);
bool initiate_comp(size_t id, volatile double curtime, double flop_count, const std::vector<int>& parameters);//int param1=-1, int param2=-1, int param3=-1, int param4=-1, int param5=-1);
void complete_comp(double errtime, size_t id, double flop_count, const std::vector<int>& parameters);//int param1=-1, int param2=-1, int param3=-1, int param4=-1, int param5=-1);
bool initiate_comm(size_t id, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
              bool is_sender=false, int partner1=-1, int user_tag1=-1, int partner2=-1, int user_tag2=-1);
bool inspect_comm(size_t id, volatile double curtime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
              bool is_sender=false, int partner=-1, int user_tag=-1);
void initiate_comm(size_t id, volatile double itime, int64_t nelem, MPI_Datatype t, MPI_Comm cm,
              MPI_Request* request, bool is_sender=false, int partner=-1, int user_tag=-1);
void complete_comm(size_t id);
int complete_comm(double curtime, MPI_Request* request, MPI_Status* status, int is_test=0, int* flag=nullptr);
int complete_comm(double curtime, int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status, int is_test=0, int* flag=nullptr);
int complete_comm(double curtime, int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[], int is_test=0);
int complete_comm(double curtime, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[], int is_test=0, int* flag=nullptr);
void propagate(MPI_Comm comm);
void collect(MPI_Comm comm);
void final_accumulate(MPI_Comm comm, double last_time);
void set_reference_values();
void save_reference_values();
void init_symbol(std::vector<std::string>& symbols);
void open_symbol(const char* symbol, double curtime);
void close_symbol(const char* symbol, double curtime);
void write_file(int variantID, int print_mode, double overhead_time);
void print(int variantID, int print_mode, double overhead_time);
void clear(int tag_count, int* distribution_tags);
void _finalize();
int get_critical_path_costs();
void get_critical_path_costs(float* costs);
int get_max_per_process_costs();
void get_max_per_process_costs(float* costs);
int get_volumetric_costs();
void get_volumetric_costs(float* costs);

}
}

#endif /*CRITTER__DISPATCH__DISPATCH_H_*/
