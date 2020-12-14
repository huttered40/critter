
#include <tuple>

namespace critter{
namespace internal{

template<size_t... indx>
class IndexPack{};

class reset_matrix_generator{
public:
  static void invoke(){}

  template<typename TupleType, typename... TupleTypes>
  static void invoke(TupleType&& t, TupleTypes... ts){
    auto matrix = std::get<0>(t);
    int m = std::get<1>(t);
    int n = std::get<2>(t);
    int lda = std::get<3>(t);
    auto reset_val = std::get<4>(t);
    for (int i=0; i<n; i++){
      memset(matrix+i*lda,1,m*sizeof(decltype(*matrix)));
    }
    std::get<5>(t)(matrix,m,n,lda);
    invoke(ts...);
  }
};

inline float flop_generator(int id, float p1=0, float p2=0, float p3=0, float p4=0, float p5=0){
  float flops=0;
  switch (id){
    case 100:
      flops = 2.*p1;
      break;
    case 101:
      flops = p1;
      break;
    case 120:
      flops = 2.*p1*p2 - (p1-p3-1)*(p1-p3) - (p2-p4-1)*(p2-p4);
      break;
    case 121:
      flops = 2.*p1*p2;
      break;
    case 122:
      flops = 2.*p1*p2;
      break;
    case 123:
      flops = p1*p1;
      break;
    case 124:
      flops = p1*(4.*p2+2) - 2.*p2*(p2+1);
      break;
    case 125:
      flops = 2.*p1*p1;
      break;
    case 126:
      flops = p1*(p1+1.);
      break;
    case 127:
      flops = 2.*p1*p1 + p1;
      break;
    case 128:
      flops = 2.*p1*p1;
      break;
    case 129:
      flops = p1*(p1+1.);
      break;
    case 130:
      flops = 2.*p1*p1 + p1;
      break;
    case 131:
      flops = p1*p1;
      break;
    case 132:
      flops = p1*p1;
      break;
    case 133:
      flops = p1*p1;
      break;
    case 134:
      flops = p1*p1;
      break;
    case 135:
      flops = p1*(2.*p2+1) - p2*(p2+1);
      break;
    case 136:
      flops = p1*(2.*p2+1) - p2*(p2+1);
      break;
    case 150:
      flops = 2.*p1*p2*p3;
      break;
    case 151:
      flops = p1*p1*p2;
      break;
    case 152:
      flops = p1*p1*p2;
      break;
    case 153:
      flops = p2*p1*(p1+1);
      break;
    case 154:
      flops = 2.*p2*p1*p1 + p1;
      break;
    case 155:
      flops = 2.*p1*p1*p2;
      break;
    case 200:
      flops = p1*p2*p2 - 1./3.*p2*p2*p2 - 1./2.*p2*p2 + 5./6.*p2;
      break;
    case 201:
      flops = 1./3.*p1*p1*p1 + 1./2.*p1*p1 + 1./6.*p1;
      break;
    case 202:
      flops = 1./3.*p1*p1*p1 + 2./3.*p1;
      break;
    case 203:
      flops = 2.*p1*p2*p2 - 2./3.*p2*p2*p2 + p1*p2 + p2*p2 + 14./3.*p2;
      break;
    case 204:
      flops = 4.*p1*p2*p3 - 2.*(p1+p2)*p3*p3 + (4./3.)*p3*p3*p3 + 3.*p2*p3 - p1*p3 - p3*p3 - 4./3.*p3;
      break;
    case 205:
      flops = 4.*p1*p2*p3 - 2.*p2*p3*p3 + 2.*p2*p3 + p1*p3 - 1./2.*p3*p3 + 1./2.*p3;
      break;
    case 206:
      flops = 4./3.*p1*p1*p1 - p1*p1 + 5./3.*p1;
      break;
    case 207:
      break;
    case 208:
      break;
  }
  return flops;
}

template<typename func_type, typename... t1_types, typename... t2_types, size_t... index_list, typename... arg_types>
inline void conditional_blas_engine(int id, int guard, std::tuple<t1_types...>&& t1, std::tuple<t2_types...>&& t2,
                             IndexPack<index_list...>, func_type* func, arg_types... args){
  if (mode && guard){
    volatile float curtime = MPI_Wtime();
    float flops = flop_generator(id,std::get<index_list>(t1)...);
    bool schedule_decision = initiate_comp(id,curtime,flops,std::get<index_list>(t2)...);
    if (schedule_decision) func(args...);
    complete_comp(0,id,flops,std::get<index_list>(t2)...);
  } else{
    func(args...);
  }
}

template<typename func_type, typename... t1_types, typename... t2_types, size_t... index_list1, typename... arg_types, size_t... index_list2, typename... TupleTypes>
inline int conditional_lapack_engine(int id, int guard, std::tuple<t1_types...>&& t1, std::tuple<t2_types...>&& t2,
                              IndexPack<index_list1...>, func_type* func, std::tuple<arg_types...>&& args, IndexPack<index_list2...>,
                              TupleTypes&&... reset_lambdas){
  if (mode && guard){
    volatile float curtime = MPI_Wtime();
    float flops = flop_generator(id,std::get<index_list1>(t1)...);
    float special_time=0;
    bool schedule_decision = initiate_comp(id,curtime,flops,std::get<index_list1>(t2)...);
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(func(std::get<index_list2>(args)...)==0);
      else{
        special_time = MPI_Wtime();
        if (reset_matrix) { reset_matrix_generator::invoke(reset_lambdas...); }
        special_time = MPI_Wtime() - special_time;
        assert(func(std::get<index_list2>(args)...)==0);
      }
    }
    complete_comp(special_time,id,flops,std::get<index_list1>(t2)...);
  } else{
    assert(func(std::get<index_list2>(args)...)==0);
  }
  return 0;// If not 0, an assert would be invoked
}

template<typename T, typename func_type>
int conditional_lapack_engine_tpqrt_(func_type* func, int matrix_layout, int m , int n , int l , int nb , T* a , int lda , T* b , int ldb , T* t , int ldt){
  if (mode && track_lapack){
    volatile float curtime = MPI_Wtime();
    float _m = m; float _n = n; float _l = l;
    float flops = 2.*_m*_n*_l;//Note: this is an educated guess. There is no information on this flop count
    float special_time=0;
    bool schedule_decision = initiate_comp(_LAPACK_tpqrt__id,curtime,flops,m,n,l,nb);
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(func(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,(i+1)*sizeof(T));// Assumes column-major
          memset(a+i*lda+i+1,0,(n-i-1)*sizeof(T));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          memset(b+i*ldb,1,(i+1)*sizeof(T));// Assumes column-major
          memset(b+i*ldb+i+1,0,(std::max(0,(m-l)-i-1))*sizeof(T));// Assumes column-major
          memset(b+i*ldb+(m-l),1,std::min(i+1,l)*sizeof(T));// Assumes column-major
          memset(b+i*ldb+(m-l)+i+1,0,std::max(0,(l-i-1))*sizeof(T));// Assumes column-major
        }
        special_time = MPI_Wtime() - special_time;
        assert(func(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt)==0);
      }
    }
    complete_comp(special_time,_LAPACK_tpqrt__id,flops,m,n,l,nb);
  } else{
    assert(func(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt)==0);
  }
  return 0;// If not 0, an assert would be invoked
}

template<typename T, typename func_type>
int conditional_lapack_engine_tpmqrt_(func_type* func, int matrix_layout, char side , char trans , int m , int n , int k , int l , int nb , const T* v ,
               int ldv , const T* t , int ldt , T* a , int lda , T* b , int ldb){
  if (mode && track_lapack){
    volatile float curtime = MPI_Wtime();
    float _m = m; float _n = n; float _k = k;
    float flops = 2.*_m*_n*_k;//Note: this is an educated guess. There is no information on this flop count
    float special_time=0;
    bool schedule_decision = initiate_comp(_LAPACK_tpmqrt__id,curtime,flops,m,n,k,l,nb);
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(func(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb)==0);
      else{
        special_time = MPI_Wtime();
        T* v_temp = (T*)v;
        T* t_temp = (T*)t;
        if (side == 'L'){
          for (int i=0; i<k; i++){
            memset(v_temp+i*ldv,1,m*sizeof(T));// Assumes column-major
          }
        } else{
          for (int i=0; i<n; i++){
            memset(v_temp+i*ldv,1,n*sizeof(T));// Assumes column-major
          }
        }
        for (int i=0; i<k; i++){
          memset(t_temp+i*ldt,1,nb*sizeof(T));// Assumes column-major
        }
        if (side=='L'){
          for (int i=0; i<n; i++){
            memset(a+i*lda,1,k*sizeof(T));// Assumes column-major
          }
        } else{
          for (int i=0; i<k; i++){
            memset(a+i*lda,1,m*sizeof(T));// Assumes column-major
          }
        }
        for (int i=0; i<n; i++){
          memset(b+i*ldb,1,m*sizeof(T));// Assumes column-major
        }
        special_time = MPI_Wtime() - special_time;
        assert(func(matrix_layout,side,trans,m,n,k,l,nb,v_temp,ldv,t_temp,ldt,a,lda,b,ldb)==0);
      }
    }
    complete_comp(special_time,_LAPACK_tpmqrt__id,flops,m,n,k,l,nb);
  } else{
    assert(func(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb)==0);
  }
  return 0;// If not 0, an assert would be invoked
}

}
}
