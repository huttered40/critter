#ifndef CRITTER__DISCRETIZATION__PARAMETERIZATION__PARAMETERIZATION_H_
#define CRITTER__DISCRETIZATION__PARAMETERIZATION__PARAMETERIZATION_H_

#include <assert.h>

namespace critter{
namespace internal{
namespace discretization{

// ****************************************************************************************************************************************************
struct comm_pattern_key{

  comm_pattern_key(int _pattern_index=0, int _tag=0, int _comm_size=0, int _comm_color=0, double _msg_size=0, int _partner_offset=0);
  comm_pattern_key(const comm_pattern_key& _copy);
  comm_pattern_key& operator=(const comm_pattern_key& _copy);
  //virtual friend bool operator==(const comm_pattern_key& ref1, const comm_pattern_key& ref2){assert(0); return false;}
  //virtual friend bool operator<(const comm_pattern_key& ref1, const comm_pattern_key& ref2){assert(0); return false;}

  virtual bool compare(const comm_pattern_key& key) const{
    assert(0);
  }

  int tag;
  int comm_size;
  int comm_color;
  int partner_offset;
  int pattern_index;
  double msg_size;
};

// ****************************************************************************************************************************************************
struct comm_pattern_key_param1 : public comm_pattern_key{
  comm_pattern_key_param1(int _pattern_index=0, int _tag=0, int _comm_size=0, int _comm_color=0, double _msg_size=0, int _partner_offset=0) : comm_pattern_key(_pattern_index,_tag,_comm_size,_comm_color,_msg_size,_partner_offset){}
  //comm_pattern_key_param1(const comm_pattern_key_param1& _copy);
  //comm_pattern_key_param1& operator=(const comm_pattern_key_param1& _copy);
/*
  friend bool operator==(const comm_pattern_key_param1& ref1, const comm_pattern_key_param1& ref2);
  friend bool operator<(const comm_pattern_key_param1& ref1, const comm_pattern_key_param1& ref2);
*/
  virtual bool compare(const comm_pattern_key_param1& key) const;
};

// ****************************************************************************************************************************************************
struct comm_pattern_key_param2 : public comm_pattern_key{
  comm_pattern_key_param2(int _pattern_index=0, int _tag=0, int _comm_size=0, int _comm_color=0, double _msg_size=0, int _partner_offset=0);
  comm_pattern_key_param2(const comm_pattern_key_param2& _copy);
  comm_pattern_key_param2& operator=(const comm_pattern_key_param2& _copy);
/*
  friend bool operator==(const comm_pattern_key_param1& ref1, const comm_pattern_key_param1& ref2);
  friend bool operator<(const comm_pattern_key_param1& ref1, const comm_pattern_key_param1& ref2);
*/
  virtual bool compare(const comm_pattern_key_param2& key) const;
};

// ****************************************************************************************************************************************************
struct comp_pattern_key{
  comp_pattern_key(int _pattern_index=0, int _tag=0, double _flops=0, int _param1=0, int _param2=0, int _param3=0, int _param4=0, int _param5=0);
  comp_pattern_key(const comp_pattern_key& _copy);
  comp_pattern_key& operator=(const comp_pattern_key& _copy);
  //virtual friend bool operator==(const comp_pattern_key& ref1, const comp_pattern_key& ref2){assert(0); return false;}
  //virtual friend bool operator<(const comp_pattern_key& ref1, const comp_pattern_key& ref2){assert(0); return false;}

  virtual bool compare(const comp_pattern_key& key) const{
    assert(0);
  }

  int tag;
  int param1,param2,param3,param4,param5;
  int pattern_index;
  double flops;
};

// ****************************************************************************************************************************************************
struct comp_pattern_key_param1 : public comp_pattern_key{

  comp_pattern_key_param1(int _pattern_index=0, int _tag=0, double _flops=0, int _param1=0, int _param2=0, int _param3=0, int _param4=0, int _param5=0) : comp_pattern_key(_pattern_index,_tag,_flops,_param1,_param2,_param3,_param4,_param5) {}
  //comp_pattern_key_param1(const comp_pattern_key_param1& _copy);
  //comp_pattern_key_param1& operator=(const comp_pattern_key_param1& _copy);
/*
  friend bool operator==(const comp_pattern_key_param1& ref1, const comp_pattern_key_param1& ref2);
  friend bool operator<(const comp_pattern_key_param1& ref1, const comp_pattern_key_param1& ref2);
*/
  virtual bool compare(const comp_pattern_key_param1& key) const;
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
  double num_scheduled_units;
  double num_non_scheduled_units;
  double M1,M2;
  double total_exec_time;
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

#endif /*CRITTER__DISCRETIZATION__PARAMETERIZATION__PARAMETERIZATION_H_*/
