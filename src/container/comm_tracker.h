#ifndef CRITTER__CONTAINER__COMM_TRACKER_H_
#define CRITTER__CONTAINER__COMM_TRACKER_H_

#include <mpi.h>
#include <unordered_map>
#include <map>
#include <vector>
#include <functional>

namespace internal{

/* \brief encapsulation of the state of a MPI routine */
class comm_tracker{
  public: 
    /* \brief name of MPI routine */
    std::string name;
    /* \brief integer tag of MPI routine */
    int tag;
    /* \brief local duration of synchronization time */
    float* my_synch_time;
    /* \brief local duration of communication time */
    float* my_comm_time;
    /* \brief local comm cost in #messages */
    float* my_msg_count;
    /* \brief local comm cost in #words */
    float* my_wrd_count;
    /* \brief duration of synchronization time along a critical path */
    float* cp_synch_time;
    /* \brief duration of communication time along a critical path */
    float* cp_comm_time;
    /* \brief comm cost in #messages along a critical path */
    float* cp_msg_count;
    /* \brief comm cost in #words along a critical path */
    float* cp_wrd_count;
    /* \brief function for cost model of MPI routine in bsp cost model, takes (msg_size_in_bytes, number_processors) and returns (latency_cost, bandwidth_cost) */
    std::function< std::pair<float,float>(int64_t,int) > cost_func_bsp;
    /* \brief function for cost model of MPI routine in alpha-beta cost model, takes (msg_size_in_bytes, number_processors) and returns (latency_cost, bandwidth_cost) */
    std::function< std::pair<float,float>(int64_t,int) > cost_func_alphabeta;
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
    /** \brief initialization of state called by construtors */
    void init();
    /** */
    void set_header(std::map<std::string,std::vector<float>>& save_info);
    /** */
    void set_cp_costs(std::map<std::string,std::vector<float>>& save_info, size_t idx);
    /** */
    void set_pp_costs(std::map<std::string,std::vector<float>>& save_info, size_t idx);
    /** */
    void set_vol_costs(std::map<std::string,std::vector<float>>& save_info);
    /** \brief sets data members to point into global arrays */
    void set_cost_pointers();
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
            std::function< std::pair<float,float>(int64_t,int)> 
              cost_func_bsp = [](int64_t n, int p){ return std::pair<float,float>(1.,n); },
            std::function< std::pair<float,float>(int64_t,int)> 
              cost_func_alphabeta = [](int64_t n, int p){ return std::pair<float,float>(1.,n); }
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
            std::function< std::pair<float,float>(int64_t,int)> 
              cost_func_bsp = [](int64_t n, int p){ return std::pair<float,float>(1.,n); },
            std::function< std::pair<float,float>(int64_t,int)> 
              cost_func_alphabeta = [](int64_t n, int p){ return std::pair<float,float>(1.,n); }
               );
    /** \brief copy constructor */
    nonblocking(nonblocking const& t);
};

struct nonblocking_info{
  nonblocking_info(){}
  nonblocking_info(float* _path_data, MPI_Request _prop_req,
                   bool _is_sender, int _partner, MPI_Comm _comm, float _nbytes,
                   float _comm_size, nonblocking* _track){
    this->path_data = _path_data;
    //this->barrier_req = _barrier_req;
    this->prop_req = _prop_req;
    this->is_sender = _is_sender;
    this->partner = _partner;
    this->comm = _comm;
    this->nbytes = _nbytes;
    this->comm_size = _comm_size;
    this->track = _track;
  }
  float* path_data;
  //MPI_Request barrier_req;
  MPI_Request prop_req;
  bool is_sender;
  int partner;
  MPI_Comm comm;
  float nbytes;
  float comm_size;
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
         _MPI_Alltoallv,
         _MPI_Comm_split,
         _MPI_Comm_dup;
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
constexpr auto list_size=35;
extern comm_tracker* list[list_size];
extern std::map<MPI_Request,nonblocking_info> nonblocking_internal_info;

}

#endif /*CRITTER__CONTAINER__COMM_TRACKER_H_*/
