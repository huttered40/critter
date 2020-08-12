#include "parameterization.h"

namespace critter{
namespace internal{
namespace discretization{

// ****************************************************************************************************************************************************
comm_pattern_key::comm_pattern_key(int _pattern_index, int _tag, int _comm_size, int _comm_color, double _msg_size, int _partner_offset){
  this->pattern_index = _pattern_index;
  this->tag = _tag;
  this->comm_size = _comm_size;
  this->comm_color = _comm_color;
  this->msg_size = _msg_size;
  this->partner_offset = _partner_offset;
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

// ****************************************************************************************************************************************************
/*
bool operator==(const comm_pattern_key_param1& ref1, const comm_pattern_key_param1& ref2){
  if ((ref1.tag==ref2.tag) && (ref1.comm_size == ref2.comm_size) && (ref1.comm_color == ref2.comm_color) && (ref1.msg_size == ref2.msg_size) && (ref1.partner_offset == ref2.partner_offset)) return true;
  else return false;
}

bool operator<(const comm_pattern_key_param1& ref1, const comm_pattern_key_param1& ref2){
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
*/
bool comm_pattern_key_param1::compare(const comm_pattern_key_param1& key) const{
  if (this->tag < key.tag) return true;
  else if (this->tag > key.tag) return false;
  if (this->comm_size < key.comm_size) return true;
  else if (this->comm_size > key.comm_size) return false;
  if (this->comm_color < key.comm_color) return true;
  else if (this->comm_color > key.comm_color) return false;
  if (this->msg_size < key.msg_size) return true;
  else if (this->msg_size > key.msg_size) return false;
  if (this->partner_offset < key.partner_offset) return true;
  else if (this->partner_offset > key.partner_offset) return false;
  return false;
}

// ****************************************************************************************************************************************************
bool comm_pattern_key_param2::compare(const comm_pattern_key_param2& key) const{
  if (this->tag < key.tag) return true;
  else if (this->tag > key.tag) return false;
  if (this->comm_size < key.comm_size) return true;
  else if (this->comm_size > key.comm_size) return false;
  if (this->comm_color < key.comm_color) return true;
  else if (this->comm_color > key.comm_color) return false;
  if (this->msg_size < key.msg_size) return true;
  else if (this->msg_size > key.msg_size) return false;
  return false;
}

// ****************************************************************************************************************************************************
comp_pattern_key::comp_pattern_key(int _pattern_index, int _tag, double _flops, int _param1, int _param2, int _param3, int _param4, int _param5){
  this->pattern_index = _pattern_index;
  this->tag = _tag;
  this->flops = _flops;
  this->param1 = _param1;
  this->param2 = _param2;
  this->param3 = _param3;
  this->param4 = _param4;
  this->param5 = _param5;
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

// ****************************************************************************************************************************************************
/*
bool operator==(const comp_pattern_key_param1& ref1, const comp_pattern_key_param1& ref2){
  if ((ref1.tag==ref2.tag) && (ref1.param1 == ref2.param1) && (ref1.param2 == ref2.param2) && (ref1.param3 == ref2.param3) && (ref1.param4 == ref2.param4) && (ref1.param5 == ref2.param5)) return true;
  else return false;
}

bool operator<(const comp_pattern_key_param1& ref1, const comp_pattern_key_param1& ref2){
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
*/
bool comp_pattern_key_param1::compare(const comp_pattern_key_param1& key) const{
  if (this->tag < key.tag) return true;
  else if (this->tag > key.tag) return false;
  if (this->param1 < key.param1) return true;
  else if (this->param1 > key.param1) return false;
  if (this->param2 < key.param2) return true;
  else if (this->param2 > key.param2) return false;
  if (this->param3 < key.param3) return true;
  else if (this->param3 > key.param3) return false;
  if (this->param4 < key.param4) return true;
  else if (this->param4 > key.param4) return false;
  if (this->param5 < key.param5) return true;
  else if (this->param5 > key.param5) return false;
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
  this->M1=0; this->M2=0;
  //this->M3=0; this->M4=0;
}

pattern::pattern(const pattern& _copy){
  this->total_exec_time = _copy.total_exec_time;
  this->steady_state = _copy.steady_state;
  this->global_steady_state = _copy.global_steady_state;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = _copy.num_non_schedules;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_non_scheduled_units = _copy.num_non_scheduled_units;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  //this->M3=0; this->M4=0;
}

pattern& pattern::operator=(const pattern& _copy){
  this->total_exec_time = _copy.total_exec_time;
  this->steady_state = _copy.steady_state;
  this->global_steady_state = _copy.global_steady_state;
  this->num_schedules = _copy.num_schedules;
  this->num_non_schedules = _copy.num_non_schedules;
  this->num_scheduled_units = _copy.num_scheduled_units;
  this->num_non_scheduled_units = _copy.num_non_scheduled_units;
  this->M1 = _copy.M1;
  this->M2 = _copy.M2;
  //this->M3=0; this->M4=0;
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
