#ifndef CRITTER__SKELETONIZE__RECORD__RECORD_H_
#define CRITTER__SKELETONIZE__RECORD__RECORD_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace skeletonize{

class record{
public:
  static void write_file(int variantID, int print_mode, float overhead_time);
  static void print(int variantID, int print_mode, float overhead_time);
};

}
}
}

#endif /*CRITTER__SKELETONIZE__RECORD__RECORD_H_*/
