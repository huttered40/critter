#ifndef CRITTER__DISCRETIZATION__CONTAINER__COMM_TRACKER_H_
#define CRITTER__DISCRETIZATION__CONTAINER__COMM_TRACKER_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

/* \brief encapsulation of the state of a MPI routine */
class comm_tracker{
  public: 
    /* \brief name of MPI routine */
    std::string name;
    /* \brief integer tag of MPI routine */
    int tag;
    /* \brief function for cost model of MPI routine in bsp cost model, takes (msg_size_in_bytes, number_processors) and returns (latency_cost, bandwidth_cost) */
    std::function< std::pair<double,double>(int64_t,int) > cost_func_bsp;
    /* \brief function for cost model of MPI routine in alpha-beta cost model, takes (msg_size_in_bytes, number_processors) and returns (latency_cost, bandwidth_cost) */
    std::function< std::pair<double,double>(int64_t,int) > cost_func_alphabeta;
    /* \brief duration of computation time for each call made locally, used to save the local computation time between calls to ::start and ::stop variants */
    double comp_time;
    /* \brief time when start() was last called, set to -1.0 initially and after stop() */
    volatile double start_time;
    /* \brief time when start() was last called, set to -1.0 initially and after stop() */
    volatile double synch_time;
    /* \brief save barrier time across start_synch */
    volatile double barrier_time;
    /* \brief cm with which start() was last called */
    MPI_Comm comm;
    /* \brief partner with which start() was last called */
    int partner1;
    /* \brief partner with which start() was last called */
    int partner2;
    /* \brief nbytes with which start() was last called */
    int64_t nbytes;
    /* \brief process count with which start() was last called */
    int comm_size;
    /* \brief is_sender bool with which start() was last called */
    bool is_sender;
    bool should_propagate;
    bool aggregate_comp_kernels;
    bool aggregate_comm_kernels;
    std::vector<comm_kernel_key> save_comm_key;
    std::vector<comp_kernel_key> save_comp_key;
    /** \brief initialization of state called by construtors */
    void init();
    /** */
    void set_header();
};

class blocking : public comm_tracker{
public:
    /**
     * \brief constructor
     * \param[in] name symbol name of MPI routine
     * \param[in] tag integer id of MPI routine
     * \param[in] cost_func_bsp function for simple cost model of MPI routine
     * \param[in] cost_func_alphabeta function for alpha-beta cost model of MPI routine assuming (synchronization-efficient collective communication algorithms)
     */
    blocking(std::string name, int tag,
            std::function< std::pair<double,double>(int64_t,int)> 
              cost_func_bsp = [](int64_t n, int p){ return std::pair<double,double>(1.,n); },
            std::function< std::pair<double,double>(int64_t,int)> 
              cost_func_alphabeta = [](int64_t n, int p){ return std::pair<double,double>(1.,n); }
            );
    /** \brief copy constructor */
    blocking(blocking const& t);
};

class nonblocking : public comm_tracker{
public:
    /**
     * \brief constructor
     * \param[in] name symbol name of MPI routine
     * \param[in] tag integer id of MPI routine
     * \param[in] cost_func_bsp function for simple cost model of MPI routine
     * \param[in] cost_func_alphabeta function for alpha-beta cost model of MPI routine assuming (synchronization-efficient collective communication algorithms)
     */
    nonblocking(std::string name, int tag,
            std::function< std::pair<double,double>(int64_t,int)> 
              cost_func_bsp = [](int64_t n, int p){ return std::pair<double,double>(1.,n); },
            std::function< std::pair<double,double>(int64_t,int)> 
              cost_func_alphabeta = [](int64_t n, int p){ return std::pair<double,double>(1.,n); }
               );
    /** \brief copy constructor */
    nonblocking(nonblocking const& t);
};

struct nonblocking_info{
  nonblocking_info(){}
  nonblocking_info(bool _is_active, bool _is_sender, int _partner, MPI_Comm _comm, double _nbytes, double _num_elems, double _comm_size, int _tag, nonblocking* _track){
    this->is_active = _is_active;
    this->is_sender = _is_sender;
    this->partner = _partner;
    this->comm = _comm;
    this->nbytes = _nbytes;
    this->num_elems = _num_elems;
    this->comm_size = _comm_size;
    this->tag = _tag;
    this->track = _track;
  }
  bool is_active;
  bool is_sender;
  int partner;
  MPI_Comm comm;
  double nbytes;
  double num_elems;
  double comm_size;
  int tag;
  nonblocking* track;
};

extern blocking
         _MPI_Send,
         _MPI_Ssend,
         _MPI_Bsend,
         _MPI_Recv,
         _MPI_Sendrecv,
         _MPI_Sendrecv_replace,
         _MPI_Barrier,
         _MPI_Bcast,
         _MPI_Reduce,
         _MPI_Allreduce,
         _MPI_Gather,
         _MPI_Allgather,
         _MPI_Scatter,
         _MPI_Reduce_scatter,
         _MPI_Alltoall,
         _MPI_Gatherv,
         _MPI_Allgatherv,
         _MPI_Scatterv,
         _MPI_Alltoallv;
extern nonblocking
         _MPI_Isend,
         _MPI_Irecv,
         _MPI_Ibcast,
         _MPI_Iallreduce,
         _MPI_Ireduce,
         _MPI_Igather,
         _MPI_Igatherv,
         _MPI_Iallgather,
         _MPI_Iallgatherv,
         _MPI_Iscatter,
         _MPI_Iscatterv,
         _MPI_Ireduce_scatter,
         _MPI_Ialltoall,
         _MPI_Ialltoallv;
constexpr auto list_size=33;
extern comm_tracker* list[list_size];
extern std::map<MPI_Request,nonblocking_info> nonblocking_internal_info;

}
}
}

#endif /*CRITTER__DISCRETIZATION__CONTAINER__COMM_TRACKER_H_*/
