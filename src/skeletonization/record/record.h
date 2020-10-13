#ifndef CRITTER__SKELETONIZATION__RECORD__RECORD_H_
#define CRITTER__SKELETONIZATION__RECORD__RECORD_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace skeletonization{

class record{
public:
  static void write_file(int variantID, int print_mode, double overhead_time);
  static void print(int variantID, int print_mode, double overhead_time);
};

}
}
}

#endif /*CRITTER__SKELETONIZATION__RECORD__RECORD_H_*/
