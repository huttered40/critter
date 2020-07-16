#include "path.h"

namespace critter{
namespace internal{
namespace optimization{

void replay(){

  std::map<std::string,double> scale_map;//TODO: Initialize this to all 1s. These will change after each simulation. These are not the same as the gradient scales.
  std::map<int,std::pair<MPI_Request,double*>> req_map;
  // Simulation loop
  for (auto i=0; i<opt_max_iter; i++){
    req_map.clear();
    // iterate over all entries in event_list
    for (auto& comm_it : event_list){
      if (comm_it.tag == -1){
        // Non-communication event. Handle later.
      }
      else if (comm_it.tag < 13){
        // Blocking collective or synchronous barrier
        // Note we will first want to use MAXLOC to determine the roots of each gradient, then zero out non-path-root entries and use an Allreduce (so 2-stage)
        double temp=0;
        PMPI_Allreduce(MPI_IN_PLACE, &temp, 1, MPI_DOUBLE, MPI_SUM,comm_it.comm);
      }
      else if (comm_it.tag < 15){
        // Blocking sendrecv
        double temp1=0; double temp2=0;
        PMPI_Sendrecv(&temp1,1,MPI_DOUBLE,comm_it.partner1,internal_tag5,&temp2,1,MPI_DOUBLE,comm_it.partner2,internal_tag5,comm_it.comm,MPI_STATUS_IGNORE);
      }
      else if ((comm_it.tag < 18) || (comm_it.tag == 32)){
        // Blocking send or recv (various sending protocols)
        double temp1=0; double temp2=0;
        if (comm_it.is_eager){
          if (comm_it.is_sender){
            PMPI_Send(&temp1,1,MPI_DOUBLE,comm_it.partner1,internal_tag5,comm_it.comm);
          } else{
            PMPI_Recv(&temp1,1,MPI_DOUBLE,comm_it.partner1,internal_tag5,comm_it.comm,MPI_STATUS_IGNORE);
          }
        } else{
          PMPI_Sendrecv(&temp1,1,MPI_DOUBLE,comm_it.partner1,internal_tag5,&temp2,1,MPI_DOUBLE,comm_it.partner1,internal_tag5,comm_it.comm,MPI_STATUS_IGNORE);
        }
      }
      else if (comm_it.tag < 20){
        // Nonblocking send or recv -> no branching on if eager or not
        if (comm_it.is_close){//req_map.find(comm_it.match) != req_map.end()){
          std::vector<MPI_Request> req_vec(comm_it.match_size);
          for (auto j=0; j<comm_it.match_size; j++){
            req_vec[j] = req_map[comm_it.match_vec[j]].first;
          }
          PMPI_Waitall(comm_it.match_size,&req_vec[0],MPI_STATUSES_IGNORE);
          for (auto j=0; j<comm_it.match_size; j++){
            free(req_map[comm_it.match_vec[j]].second);
          }
        } else{
          MPI_Request req;
          double* temp = (double*)malloc(1*sizeof(double));
          if (comm_it.is_sender){
            PMPI_Isend(temp,1,MPI_DOUBLE,comm_it.partner1,internal_tag5,comm_it.comm,&req);
            req_map[comm_it.match_id] = std::make_pair(req,temp);
          } else{
            PMPI_Irecv(temp,1,MPI_DOUBLE,comm_it.partner1,internal_tag5,comm_it.comm,&req);
            req_map[comm_it.match_id] = std::make_pair(req,temp);
          }
        }
      }
      else if (comm_it.tag < 32){
        // Nonblocking collective
      }
    }
    // reset the data structures
    event_list.clear();
    event_list_size=0;
  }
}

}
}
}
