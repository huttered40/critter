#ifndef CRITTER__DISCRETIZATION__PARAMETERIZATION__PARAMETERIZATION_H_
#define CRITTER__DISCRETIZATION__PARAMETERIZATION__PARAMETERIZATION_H_

#include <functional>
#include <utility>
#include <vector>

namespace critter{
namespace internal{
namespace discretization{

extern int comm_envelope_param;
extern int comp_envelope_param;
extern int comm_unit_param;
extern int comp_unit_param;
extern int comm_analysis_param;
extern int comp_analysis_param;

// ****************************************************************************************************************************************************
struct pattern_key{};

// ****************************************************************************************************************************************************
struct comm_pattern_key : public pattern_key{

  comm_pattern_key(int _rank=-1, int _pattern_index=-1, int _tag=-1, int _dim_sizes[2]=nullptr, int _dim_strides[2]=nullptr, double _msg_size=-1, int _partner=-1);
  comm_pattern_key(int _pattern_index, int _tag, int _dim_sizes[2], int _dim_strides[2], double _msg_size, int _partner_offset);
  comm_pattern_key(const comm_pattern_key& _copy);
  comm_pattern_key& operator=(const comm_pattern_key& _copy);
  friend bool operator==(const comm_pattern_key& ref1, const comm_pattern_key& ref2);
  friend bool operator<(const comm_pattern_key& ref1, const comm_pattern_key& ref2);

  int tag;
  int dim_sizes[2];// Allow up to 2 dimensions
  int dim_strides[2];// Allow up to 2 dimensions
  int partner_offset;
  int pattern_index;
  double msg_size;
};


// ****************************************************************************************************************************************************
struct comp_pattern_key : public pattern_key{

  comp_pattern_key(int _pattern_index=-1, int _tag=-1, double _flops=-1, int =-1, int _param2=-1, int _param3=-1, int _param4=-1, int _param5=-1);
  comp_pattern_key(const comp_pattern_key& _copy);
  comp_pattern_key& operator=(const comp_pattern_key& _copy);
  friend bool operator==(const comp_pattern_key& ref1, const comp_pattern_key& ref2);
  friend bool operator<(const comp_pattern_key& ref1, const comp_pattern_key& ref2);

  int tag;
  int param1,param2,param3,param4,param5;
  int pattern_index;
  double flops;
};

}
}
}

/*
namespace std{
  template <>
  struct hash<critter::internal::discretization::comm_pattern_key>{
    std::size_t operator()(const critter::internal::discretization::comm_pattern_key& k) const{
      using std::size_t;
      using std::hash;

      // Compute individual hash values for first,
      // second and third and combine them using XOR and bit shifting:
      size_t res = 17;
      res = res*23+hash<int>()(k.tag);
      res = res*23+hash<int>()(k.comm_size);
      res = res*23+hash<int>()(k.comm_color);
      res = res*23+hash<int>()(k.partner_offset);
      res = res*23+hash<double>()(k.msg_size);
      return res;
//
//      return ((hash<int>()(k.tag)
//           ^ (hash<int>()(k.comm_size) << 1)) >> 1)
//           ^ (hash<int>()(k.comm_color) << 1)
//           ^ (hash<int>()(k.partner_offset) << 2)
//           ^ (hash<double>()(k.msg_size) << 3);
//
    }
  };
  template <>
  struct hash<critter::internal::discretization::comp_pattern_key>{
    std::size_t operator()(const critter::internal::discretization::comp_pattern_key& k) const{
      using std::size_t;
      using std::hash;

      // Compute individual hash values for first,
      // second and third and combine them using XOR and bit shifting:
      size_t res = 17;
      res = res*23+hash<int>()(k.tag);
      res = res*23+hash<int>()(k.param1);
      res = res*23+hash<int>()(k.param2);
      res = res*23+hash<int>()(k.param3);
      res = res*23+hash<int>()(k.param4);
      res = res*23+hash<int>()(k.param5);
      res = res*23+hash<double>()(k.flops);
      return res;
//
//      return ((hash<int>()(k.tag)
//           ^ (hash<int>()(k.param1) << 1)) >> 1)
//           ^ (hash<int>()(k.param2) << 1)
//           ^ (hash<int>()(k.param3) << 2)
//           ^ (hash<int>()(k.param4) << 3)
//           ^ (hash<int>()(k.param5) << 4)
//           ^ (hash<double>()(k.flops) << 5);
//
    }
  };
};
*/

#endif /*CRITTER__DISCRETIZATION__PARAMETERIZATION__PARAMETERIZATION_H_*/
