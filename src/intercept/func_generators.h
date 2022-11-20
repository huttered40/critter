#include <tuple>

#include "../profiler.h"
#include "../util.h"

namespace internal{

template<size_t... indx>
class IndexPack{};

inline float flop_generator(int id, float p1=0, float p2=0, float p3=0, float p4=0, float p5=0){
  float flops=0;
  switch (id){
    case 0:
      flops = 2.*p1;
      break;
    case 1:
      flops = p1;
      break;
    case 20:
      flops = 2.*p1*p2 - (p1-p3-1)*(p1-p3) - (p2-p4-1)*(p2-p4);
      break;
    case 21:
      flops = 2.*p1*p2;
      break;
    case 22:
      flops = p1*p1;
      break;
    case 23:
      flops = p1*(4.*p2+2) - 2.*p2*(p2+1);
      break;
    case 24:
      flops = 2.*p1*p1;
      break;
    case 25:
      flops = p1*(p1+1.);
      break;
    case 26:
      flops = 2.*p1*p1 + p1;
      break;
    case 27:
      flops = 2.*p1*p1;
      break;
    case 28:
      flops = p1*(p1+1.);
      break;
    case 29:
      flops = 2.*p1*p1 + p1;
      break;
    case 30:
      flops = p1*p1;
      break;
    case 31:
      flops = p1*p1;
      break;
    case 32:
      flops = p1*p1;
      break;
    case 33:
      flops = p1*p1;
      break;
    case 34:
      flops = p1*(2.*p2+1) - p2*(p2+1);
      break;
    case 35:
      flops = p1*(2.*p2+1) - p2*(p2+1);
      break;
    case 50:
      flops = 2.*p1*p2*p3;
      break;
    case 51:
      flops = p1*p1*p2;
      break;
    case 52:
      flops = p1*p1*p2;
      break;
    case 53:
      flops = p2*p1*(p1+1);
      break;
    case 54:
      flops = 2.*p2*p1*p1 + p1;
      break;
    case 55:
      flops = 2.*p1*p1*p2;
      break;
    case 100:
      flops = p1*p2*p2 - 1./3.*p2*p2*p2 - 1./2.*p2*p2 + 5./6.*p2;
      break;
    case 101:
      flops = 1./3.*p1*p1*p1 + 1./2.*p1*p1 + 1./6.*p1;
      break;
    case 102:
      flops = 1./3.*p1*p1*p1 + 2./3.*p1;
      break;
    case 103:
      flops = 2.*p1*p2*p2 - 2./3.*p2*p2*p2 + p1*p2 + p2*p2 + 14./3.*p2;
      break;
    case 104:
      flops = 4.*p1*p2*p3 - 2.*(p1+p2)*p3*p3 + (4./3.)*p3*p3*p3 + 3.*p2*p3 - p1*p3 - p3*p3 - 4./3.*p3;
      break;
    case 105:
      flops = 4.*p1*p2*p3 - 2.*p2*p3*p3 + 2.*p2*p3 + p1*p3 - 1./2.*p3*p3 + 1./2.*p3;
      break;
    case 106:
      flops = 4./3.*p1*p1*p1 - p1*p1 + 5./3.*p1;
      break;
    case 107:
      flops = 2*p1*p2*p3;// guess
      break;
    case 108:
      flops = 2*p1*p2*p3;// guess
      break;
  }
  return flops;
}

template<typename func_type, typename... t1_types, typename... t2_types, size_t... index_list, typename... arg_types>
inline int engine(size_t id, int guard, std::tuple<t1_types...>&& t1, std::tuple<t2_types...>&& t2,
                             IndexPack<index_list...>, func_type* func, arg_types... args){
  if (mode && guard){
    volatile auto curtime = MPI_Wtime();
    float flops = flop_generator(id,std::get<index_list>(t1)...);
    bool schedule_decision = profiler::initiate_comp(curtime,flops);
    if (schedule_decision) func(args...);//NOTE: Cannot place assert(func(args...))==0 because BLAS is invoked here too (not just LAPACK)
    profiler::complete_comp();
  } else{
    func(args...);//NOTE: Cannot place assert(func(args...))==0 because BLAS is invoked here too (not just LAPACK)
  }
  return 0;
}

}
