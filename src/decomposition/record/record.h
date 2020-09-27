#ifndef CRITTER__DECOMPOSITION__RECORD__RECORD_H_
#define CRITTER__DECOMPOSITION__RECORD__RECORD_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace decomposition{

class record{
public:
  static void invoke(std::ofstream& Stream, int variantID, double overhead_time);
  static void invoke(std::ostream& Stream, int variantID, double overhead_time);
};

}
}
}

#endif /*CRITTER__DECOMPOSITION__RECORD__RECORD_H_*/
