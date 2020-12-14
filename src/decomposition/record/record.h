#ifndef CRITTER__DECOMPOSITION__RECORD__RECORD_H_
#define CRITTER__DECOMPOSITION__RECORD__RECORD_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace decomposition{

class record{
public:
  static void write_file(int variantID, float overhead_time);
  static void print(int variantID, float overhead_time);
};

}
}
}

#endif /*CRITTER__DECOMPOSITION__RECORD__RECORD_H_*/
