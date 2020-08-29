#ifndef CRITTER__DISCRETIZATION__PARAMETERIZATION__PARAMETERIZATION_H_
#define CRITTER__DISCRETIZATION__PARAMETERIZATION__PARAMETERIZATION_H_

#include <functional>

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

  comm_pattern_key(int _rank=-1, int _pattern_index=-1, int _tag=-1, int _comm_size=-1, int _comm_color=-1, double _msg_size=-1, int _partner=-1);
  comm_pattern_key(int _pattern_index, int _tag, int _comm_size, int _comm_color, double _msg_size, int _partner_offset);
  comm_pattern_key(const comm_pattern_key& _copy);
  comm_pattern_key& operator=(const comm_pattern_key& _copy);
  friend bool operator==(const comm_pattern_key& ref1, const comm_pattern_key& ref2);
  friend bool operator<(const comm_pattern_key& ref1, const comm_pattern_key& ref2);

  int tag;
  int comm_size;
  int comm_color;
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

// ****************************************************************************************************************************************************
struct pattern{
  // If I leverage the kurtosis, I will have to utilize the arithmetic mean.
  //   Note that I'd rather utilize the geometric mean, but I'm not sure how to convert this algorithm
  //     to handle that.
  pattern();
  pattern(const pattern& _copy);
  pattern& operator=(const pattern& _copy);

  int steady_state;
  int global_steady_state;
  int num_schedules;
  int num_non_schedules;
  int num_propagations;
  int num_non_propagations;
  double num_scheduled_units;
  double num_non_scheduled_units;
  double M1,M2;
  double total_exec_time;
};

// ****************************************************************************************************************************************************
struct idle_pattern{
  // If I leverage the kurtosis, I will have to utilize the arithmetic mean.
  //   Note that I'd rather utilize the geometric mean, but I'm not sure how to convert this algorithm
  //     to handle that.
  idle_pattern();
  idle_pattern(const idle_pattern& _copy);
  idle_pattern& operator=(const idle_pattern& _copy);

  int num_schedules;
  int num_non_schedules;
  double M1,M2;
};

// ****************************************************************************************************************************************************
struct pattern_key_id{

  pattern_key_id(bool _is_active=false, int _key_index=0, int _val_index=0, bool _is_updated=false);
  pattern_key_id(const pattern_key_id& _copy);
  pattern_key_id& operator=(const pattern_key_id& _copy);

  // Active just means its still being propogated. It acts as a switch betweeh steady_state arrays and active arrays
  bool is_active;
  bool is_updated;
  int key_index;
  int val_index;
};

}
}
}

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
/*
      return ((hash<int>()(k.tag)
           ^ (hash<int>()(k.comm_size) << 1)) >> 1)
           ^ (hash<int>()(k.comm_color) << 1)
           ^ (hash<int>()(k.partner_offset) << 2)
           ^ (hash<double>()(k.msg_size) << 3);
*/
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
/*
      return ((hash<int>()(k.tag)
           ^ (hash<int>()(k.param1) << 1)) >> 1)
           ^ (hash<int>()(k.param2) << 1)
           ^ (hash<int>()(k.param3) << 2)
           ^ (hash<int>()(k.param4) << 3)
           ^ (hash<int>()(k.param5) << 4)
           ^ (hash<double>()(k.flops) << 5);
*/
    }
  };
};

#endif /*CRITTER__DISCRETIZATION__PARAMETERIZATION__PARAMETERIZATION_H_*/
