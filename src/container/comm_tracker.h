#ifndef CRITTER__CONTAINER__COMM_TRACKER_H_
#define CRITTER__CONTAINER__COMM_TRACKER_H_

#include "../util.h"

namespace critter{
namespace internal{

/* \brief encapsulation of the state of a MPI routine */
class comm_tracker{
  public: 
    /* \brief name of MPI routine */
    std::string name;
    /* \brief integer tag of MPI routine */
    int tag;
    /* \brief local duration of synchronization time */
    double* my_synch_time;
    /* \brief local duration of data mvt time */
    double* my_datamvt_time;
    /* \brief local duration of communication time */
    double* my_comm_time;
    /* \brief local comm cost in #messages */
    double* my_msg_count;
    /* \brief local comm cost in #words */
    double* my_wrd_count;
    /* \brief duration of synchronization time along a critical path */
    double* critical_path_synch_time;
    /* \brief number of time spent moving data about the network along a critical path */
    double* critical_path_datamvt_time;
    /* \brief duration of communication time along a critical path */
    double* critical_path_comm_time;
    /* \brief comm cost in #messages along a critical path */
    double* critical_path_msg_count;
    /* \brief comm cost in #words along a critical path */
    double* critical_path_wrd_count;
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
    /** \brief initialization of state called by construtors */
    void init();
    /** */
    void set_header();
    /** */
    void set_critical_path_costs(size_t idx);
    /** */
    void set_per_process_costs(size_t idx);
    /** */
    void set_volume_costs();
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
            std::function< std::pair<double,double>(int64_t,int)> 
              cost_func_bsp = [](int64_t n, int p){ return std::pair<double,double>(1.,n); },
            std::function< std::pair<double,double>(int64_t,int)> 
              cost_func_alphabeta = [](int64_t n, int p){ return std::pair<double,double>(1.,n); }
            );
    /** \brief copy constructor */
    blocking(blocking const& t);
    /** \brief starts tracking of blocking MPI communication of nelem elements of type t over communicator cm */
    void start(volatile double curTime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, bool is_sender=false, int partner1=-1, int partner2=-1);
    /** \brief starts communication timer for corresponding blocking MPI communication */
    void intermediate();
    /** \brief completes interception of blocking communication protocol */
    void stop();
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
    /** \brief starts tracking of nonblocking MPI communication of nelem elements of type t over communicator cm */
    void start(volatile double curTime, volatile double iTime, int64_t nelem, MPI_Datatype t, MPI_Comm cm, MPI_Request* request, bool is_sender=false, int partner=-1);
    /** \brief completes interception of nonblocking communication protocol */
    void stop(MPI_Request* internal_request, double comp_time, double comm_time);
};

extern blocking
         _MPI_Send,
         _MPI_Ssend,
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
constexpr auto list_size=32;
extern comm_tracker* list[list_size];
extern std::map<MPI_Request,nonblocking*> internal_comm_track;

}
}

#endif /*CRITTER__CONTAINER__COMM_TRACKER_H_*/
