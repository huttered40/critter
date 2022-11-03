#ifndef CRITTER__PROFILE__RECORD__RECORD_H_
#define CRITTER__PROFILE__RECORD__RECORD_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace profile{

class record{
public:
  static void write_file(int variantID, float overhead_time);
  static void print(int variantID, float overhead_time);
};

}
}
}

#endif /*CRITTER__PROFILE__RECORD__RECORD_H_*/
