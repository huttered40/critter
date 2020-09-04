#include "parameterization.h"

#include <stdint.h>
#include <assert.h>
#include <iostream>

namespace critter{
namespace internal{
namespace discretization{

static double truncate(double val, int unit_param){
  if (unit_param==0) return val;
  // returns the next highest power of 2 as a double
  if (val==0) return 0;
  // I'm counting on the int64_t to hold sufficiently large number, especially to store the flop_count exactly.
  int64_t v = val;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v |= v >> 32;
  v++;
  if (unit_param==1){
    return (double)v;
  }
  else{
    int pos=0; int64_t temp = v;
    while (temp){
      pos++; temp>>=1;
    }
    pos--;
    if (pos%unit_param==0) { return (double)v; }
    for (int z = pos%unit_param; z<unit_param; z++){
      v<<=1;
    }
    return (double)v;
  }
}
static double truncate(int v, int unit_param){
  if (unit_param==0) return v;
  // returns the next highest power of 2 as a double
  if (v==0) return 0;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  if (unit_param==1){
    return v;
  }
  else{
    int pos=0; int64_t temp = v;
    while (temp){
      pos++; temp>>=1;
    }
    pos--;
    if (pos%unit_param==0) { return (double)v; }
    for (int z = pos%unit_param; z<unit_param; z++){
      v<<=1;
    }
    return (double)v;
  }
}

// ****************************************************************************************************************************************************
int comm_envelope_param;
int comp_envelope_param;
int comm_unit_param;
int comp_unit_param;
int comm_analysis_param;
int comp_analysis_param;
// ****************************************************************************************************************************************************
comm_pattern_key::comm_pattern_key(int _rank, int _pattern_index, int _tag, int _comm_size, int _comm_color, double _msg_size, int _partner){
  this->pattern_index = _pattern_index;
  // Envelope (non-message-size) parameterization specification
  this->tag = _tag;
  this->comm_size = _comm_size;
  this->comm_color = _comm_color;
  if (comm_envelope_param == 0){
    this->partner_offset = _rank-_partner;
  }
  else if (comm_envelope_param == 1){
    this->partner_offset = abs(_rank-_partner);
  }
  else if (comm_envelope_param == 2){
    this->partner_offset = _partner;
  }
  else if (comm_envelope_param == 3){// Note: this parameter specifies that partner rank will not factor into parameterization (but we still need to give the member a value)
    this->partner_offset = -1;
  }
  // Unit (message-size) parameterization specification
  this->msg_size = truncate(_msg_size,comm_unit_param);
  // Regardless of the specified envelope parameterization, non-p2p communication requires all processes to set partner_offset <- -1
  if (_partner == -1){ this->partner_offset = -1; }
}

// This constructor is used when transferring ownership of kernels following path propagation.
comm_pattern_key::comm_pattern_key(int _pattern_index, int _tag, int _comm_size, int _comm_color, double _msg_size, int _partner_offset){
  this->pattern_index = _pattern_index;
  this->tag = _tag;
  this->comm_size = _comm_size;
  this->comm_color = _comm_color;
  this->partner_offset = _partner_offset;
  this->msg_size = _msg_size;
}

comm_pattern_key::comm_pattern_key(const comm_pattern_key& _copy){
  this->pattern_index = _copy.pattern_index;
  this->tag = _copy.tag;
  this->comm_size = _copy.comm_size;
  this->comm_color = _copy.comm_color;
  this->msg_size = _copy.msg_size;
  this->partner_offset = _copy.partner_offset;
}

comm_pattern_key& comm_pattern_key::operator=(const comm_pattern_key& _copy){
  this->pattern_index = _copy.pattern_index;
  this->tag = _copy.tag;
  this->comm_size = _copy.comm_size;
  this->comm_color = _copy.comm_color;
  this->msg_size = _copy.msg_size;
  this->partner_offset = _copy.partner_offset;
  return *this;
}

bool operator==(const comm_pattern_key& ref1, const comm_pattern_key& ref2){
  // Note that because of how we set the member variables in the constructor based on envlope, unit, and analysis parameterizations, no branching is required here.
  if ((ref1.tag==ref2.tag) && (ref1.comm_size == ref2.comm_size) && (ref1.comm_color == ref2.comm_color) && (ref1.msg_size == ref2.msg_size) && (ref1.partner_offset == ref2.partner_offset)) return true;
  else return false;
}

bool operator<(const comm_pattern_key& ref1, const comm_pattern_key& ref2){
  if (ref1.tag < ref2.tag) return true;
  else if (ref1.tag > ref2.tag) return false;
  if (ref1.comm_size < ref2.comm_size) return true;
  else if (ref1.comm_size > ref2.comm_size) return false;
  if (ref1.comm_color < ref2.comm_color) return true;
  else if (ref1.comm_color > ref2.comm_color) return false;
  if (ref1.msg_size < ref2.msg_size) return true;
  else if (ref1.msg_size > ref2.msg_size) return false;
  if (ref1.partner_offset < ref2.partner_offset) return true;
  else if (ref1.partner_offset > ref2.partner_offset) return false;
  return false;
}

// ****************************************************************************************************************************************************
comp_pattern_key::comp_pattern_key(int _pattern_index, int _tag, double _flops, int _param1, int _param2, int _param3, int _param4, int _param5){
  this->pattern_index = _pattern_index;
  this->tag = _tag;
  this->param1 = truncate(_param1,comp_unit_param);
  this->param2 = truncate(_param2,comp_unit_param);
  this->param3 = truncate(_param3,comp_unit_param);
  this->param4 = truncate(_param4,comp_unit_param);
  this->param5 = truncate(_param5,comp_unit_param);
  this->flops = truncate(_flops,comp_unit_param);
}

comp_pattern_key::comp_pattern_key(const comp_pattern_key& _copy){
  this->pattern_index = _copy.pattern_index;
  this->tag = _copy.tag;
  this->flops = _copy.flops;
  this->param1 = _copy.param1;
  this->param2 = _copy.param2;
  this->param3 = _copy.param3;
  this->param4 = _copy.param4;
  this->param5 = _copy.param5;
}

comp_pattern_key& comp_pattern_key::operator=(const comp_pattern_key& _copy){
  this->pattern_index = _copy.pattern_index;
  this->tag = _copy.tag;
  this->flops = _copy.flops;
  this->param1 = _copy.param1;
  this->param2 = _copy.param2;
  this->param3 = _copy.param3;
  this->param4 = _copy.param4;
  this->param5 = _copy.param5;
  return *this;
}

bool operator==(const comp_pattern_key& ref1, const comp_pattern_key& ref2){
  if ((ref1.tag==ref2.tag) && (ref1.param1 == ref2.param1) && (ref1.param2 == ref2.param2) && (ref1.param3 == ref2.param3) && (ref1.param4 == ref2.param4) && (ref1.param5 == ref2.param5)) return true;
  else return false;
}

bool operator<(const comp_pattern_key& ref1, const comp_pattern_key& ref2){
  if (ref1.tag < ref2.tag) return true;
  else if (ref1.tag > ref2.tag) return false;
  if (ref1.param1 < ref2.param1) return true;
  else if (ref1.param1 > ref2.param1) return false;
  if (ref1.param2 < ref2.param2) return true;
  else if (ref1.param2 > ref2.param2) return false;
  if (ref1.param3 < ref2.param3) return true;
  else if (ref1.param3 > ref2.param3) return false;
  if (ref1.param4 < ref2.param4) return true;
  else if (ref1.param4 > ref2.param4) return false;
  if (ref1.param5 < ref2.param5) return true;
  else if (ref1.param5 > ref2.param5) return false;
  return false;
}

// ****************************************************************************************************************************************************
pattern::pattern(){
  this->total_exec_time = 0;
  this->steady_state=0;
  this->global_steady_state=0;
  this->num_schedules = 0;
  this->num_non_schedules = 0;
  this->num_scheduled_units = 0;
  this->num_non_scheduled_units = 0;
  this->num_propagations=0;
  this->num_non_propagations=0;
  this->M1=0; this->M2=0;
}

pattern::pattern(const pattern& _copy){
  this->total_exec_time = _copy.total_exec_time;
  this->steady_state = _copy.steady_state;
  this->global_steady_state = _copy.global_steady_state;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = _copy.num_non_schedules;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_non_scheduled_units = _copy.num_non_scheduled_units;
  this->num_propagations=_copy.num_propagations;
  this->num_non_propagations=_copy.num_non_propagations;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
}

pattern& pattern::operator=(const pattern& _copy){
  this->total_exec_time = _copy.total_exec_time;
  this->steady_state = _copy.steady_state;
  this->global_steady_state = _copy.global_steady_state;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = _copy.num_non_schedules;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_non_scheduled_units = _copy.num_non_scheduled_units;
  this->num_propagations=_copy.num_propagations;
  this->num_non_propagations=_copy.num_non_propagations;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  return *this;
}

// ****************************************************************************************************************************************************
idle_pattern::idle_pattern(){
  this->M1=0; this->M2=0;
  this->num_schedules = 0;
  this->num_non_schedules = 0;
}

idle_pattern::idle_pattern(const idle_pattern& _copy){
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = _copy.num_non_schedules;
}

idle_pattern& idle_pattern::operator=(const idle_pattern& _copy){
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = _copy.num_non_schedules;
  return *this;
}

// ****************************************************************************************************************************************************
pattern_key_id::pattern_key_id(bool _is_active, int _key_index, int _val_index, bool _is_updated){
  this->is_active = _is_active;
  this->key_index = _key_index;
  this->val_index = _val_index;
  this->is_updated = _is_updated;
}

pattern_key_id::pattern_key_id(const pattern_key_id& _copy){
  this->is_active = _copy.is_active;
  this->key_index = _copy.key_index;
  this->val_index = _copy.val_index;
  this->is_updated = _copy.is_updated;
}

pattern_key_id& pattern_key_id::operator=(const pattern_key_id& _copy){
  this->is_active = _copy.is_active;
  this->key_index = _copy.key_index;
  this->val_index = _copy.val_index;
  this->is_updated = _copy.is_updated;
  return *this;
}


}
}
}
