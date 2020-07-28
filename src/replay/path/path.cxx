#include "path.h"
#include "../../decomposition/container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace replay{

void invoke(){

  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  size_t n = critter::internal::decomposition::symbol_timers.size(); // number of unique symbols
  if (n==0) return;// trivial corner case
  size_t m = num_gradient_points;
  std::vector<double> gradient_scale(m);
  std::map<std::string,std::pair<int,double>> scale_map;
  std::vector<int> gradient_save_index(n);
  std::vector<double> gradient_save_val(n);
  std::vector<double> table(m*n,0.);
  std::vector<double> partner_table(m*n,0.);
  std::map<int,std::pair<std::pair<MPI_Request,MPI_Request>,double*>> req_map;
  // Fill in gradient scale
  for (auto i=0; i<gradient_scale.size(); i++){
    gradient_scale[i] = .01*(gradient_jump_size*i);
  }
  // Fill in scale_map
  for (auto i=0; i<symbol_order.size(); i++){
    scale_map[symbol_order[i]] = std::make_pair(i,1.);
  }
  // Simulation loop
  for (auto i=0; i<opt_max_iter; i++){
    req_map.clear();
    std::fill(table.begin(),table.end(),0.);
    gradient_save_index.resize(n,0);
    gradient_save_val.resize(n,0);
    // iterate over all entries in event_list
    for (auto& comm_it : event_list){

      // Note: I've been swapping around the choice of these measurements.
      double exec_time_update = comm_it.measurements[num_per_process_measures-4];
      if (scale_map.find(comm_it.kernel) != scale_map.end()){
        exec_time_update *= scale_map[comm_it.kernel].second;
      }
      for (auto j=0; j<n; j++){
        if (scale_map.find(comm_it.kernel) == scale_map.end()){
          for (auto k=0; k<m; k++){
            table[j*m+k] += exec_time_update;
            table[j*m+k] += comm_it.measurements[num_per_process_measures-2];
          }
        }
        else if (scale_map[comm_it.kernel].first != j){
          for (auto k=0; k<m; k++){
            table[j*m+k] += exec_time_update;
            table[j*m+k] += comm_it.measurements[num_per_process_measures-2];
          }
        }
        else{
          for (auto k=0; k<m; k++){
            table[j*m+k] += std::max(exec_time_update*(1.-gradient_scale[k]),0.);
            table[j*m+k] += comm_it.measurements[num_per_process_measures-2];
          }
        }
      }

      if (comm_it.tag == -1){
        // Non-communication event.
      }
      else if (comm_it.tag < 13){
        // Blocking collective or synchronous barrier
        // Note we will first want to use MAXLOC to determine the roots of each gradient, then zero out non-path-root entries and use an Allreduce (so 2-stage)
        PMPI_Allreduce(MPI_IN_PLACE, &table[0], table.size(), MPI_DOUBLE, MPI_MAX,comm_it.comm);
      }
      else if (comm_it.tag < 15){
        // Blocking sendrecv
        PMPI_Sendrecv(&table[0],table.size(),MPI_DOUBLE,comm_it.partner1,internal_tag5,&partner_table[0],partner_table.size(),MPI_DOUBLE,comm_it.partner2,internal_tag5,comm_it.comm,MPI_STATUS_IGNORE);
        for (auto j=0; j<table.size(); j++){
          table[j] = std::max(table[j],partner_table[j]);
        }
      }
      else if ((comm_it.tag < 18) || (comm_it.tag == 32)){
        // Blocking send or recv (various sending protocols)
        if (comm_it.is_eager){
          if (comm_it.is_sender){
            PMPI_Send(&table[0],table.size(),MPI_DOUBLE,comm_it.partner1,internal_tag5,comm_it.comm);
          } else{
            PMPI_Recv(&partner_table[0],partner_table.size(),MPI_DOUBLE,comm_it.partner1,internal_tag5,comm_it.comm,MPI_STATUS_IGNORE);
            for (auto j=0; j<table.size(); j++){
              table[j] = std::max(table[j],partner_table[j]);
            }
          }
        } else{
          PMPI_Sendrecv(&table[0],table.size(),MPI_DOUBLE,comm_it.partner1,internal_tag5,&partner_table[0],partner_table.size(),MPI_DOUBLE,comm_it.partner1,internal_tag5,comm_it.comm,MPI_STATUS_IGNORE);
          for (auto j=0; j<table.size(); j++){
            table[j] = std::max(table[j],partner_table[j]);
          }
        }
      }
      else if (comm_it.tag < 20){
        // Nonblocking send or recv -> no branching on if eager or not
        if (comm_it.is_close){//req_map.find(comm_it.match) != req_map.end()){
          std::vector<MPI_Request> req_vec;//(comm_it.match_size);
          for (auto j=0; j<comm_it.match_size; j++){
            req_vec.push_back(req_map[comm_it.match_vec[j]].first.first);
            if (!comm_it.is_eager){
              req_vec.push_back(req_map[comm_it.match_vec[j]].first.second);
            }
          }
          PMPI_Waitall(req_vec.size(),&req_vec[0],MPI_STATUSES_IGNORE);
          for (auto j=0; j<comm_it.match_size; j++){
            for (auto k=0; k<table.size(); k++){
              table[k] = std::max(table[k],std::max(req_map[comm_it.match_vec[j]].second[k],req_map[comm_it.match_vec[j]].second[table.size()+k]));
            }
            free(req_map[comm_it.match_vec[j]].second);
          }
        } else{
          MPI_Request req1 = MPI_REQUEST_NULL;
          MPI_Request req2 = MPI_REQUEST_NULL;
          // Just allocate 2 instead of 1 assuming all requests use rendezvous protocol
          double* table_nblk = (double*)malloc(2*table.size()*sizeof(double));
          for (auto j=0; j<2*table.size(); j++){
            table_nblk[j] = table[j%table.size()];
          }
          if (comm_it.is_sender){
            PMPI_Isend(&table_nblk[0],table.size(),MPI_DOUBLE,comm_it.partner1,internal_tag5,comm_it.comm,&req1);
            if (!comm_it.is_eager){
              PMPI_Irecv(&table_nblk[table.size()],table.size(),MPI_DOUBLE,comm_it.partner1,internal_tag5,comm_it.comm,&req2);
            }
            req_map[comm_it.match_id] = std::make_pair(std::make_pair(req1,req2),table_nblk);
          }
          else{
            PMPI_Irecv(&table_nblk[0],table.size(),MPI_DOUBLE,comm_it.partner1,internal_tag5,comm_it.comm,&req1);
            if (!comm_it.is_eager){
              PMPI_Isend(&table_nblk[table.size()],table.size(),MPI_DOUBLE,comm_it.partner1,internal_tag5,comm_it.comm,&req2);
            }
            req_map[comm_it.match_id] = std::make_pair(std::make_pair(req1,req2),table_nblk);
          }
        }
      }
      else if (comm_it.tag < 32){
        // Nonblocking collective
      }
    }

    // I think we need one more step.
    PMPI_Allreduce(MPI_IN_PLACE, &table[0], table.size(), MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
    // Find the entry with the max runtime improvement but with the least intensity, and update scale_map.
    for (auto j=0; j<n; j++){
      gradient_save_val[j] = table[j*m];
      gradient_save_index[j] = 0;
      for (auto k=1; k<m; k++){
        if (table[j*m+k] < gradient_save_val[j]){
          gradient_save_val[j] = table[j*m+k];
          gradient_save_index[j] = k;
        }
      }
      if (rank==0){
        std::cout << "Symbol " << symbol_order[j] << " used " << gradient_save_index[j] << " jumps to improve runtime from " << table[j*m] << " to " << gradient_save_val[j] << std::endl;
      }
    }
    if (rank==0) std::cout << "\n";
    // Now iterate over each kernel's representative to identify the best candidate for optimization as we recurse into the next simulation step
    double opt_val = gradient_save_val[0];
    int opt_index = gradient_save_index[0];
    int opt_kernel = 0;
    for (auto j=1; j<n; j++){
      if (gradient_save_val[j] < opt_val){
        opt_val = gradient_save_val[j];
        opt_kernel = j;
        opt_index = gradient_save_index[j];
      }
    }
    // Update scale_map
    scale_map[symbol_order[opt_kernel]].second *= (1.-gradient_scale[opt_index]);
  }
  if (rank==0){
    for (auto j=0; j<n; j++){
      std::cout << "Optimization intensity distribution to kernel " << symbol_order[j] << " is " << scale_map[symbol_order[j]].second << std::endl << std::endl;
    }
  }
  // reset the data structures
  event_list.clear();
  event_list_size=0;
}

}
}
}
