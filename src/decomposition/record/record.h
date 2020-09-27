#ifndef CRITTER__DECOMPOSITION__RECORD__RECORD_H_
#define CRITTER__DECOMPOSITION__RECORD__RECORD_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace decomposition{

class record{
public:
  static void invoke(std::ofstream& Stream, int variantID);
  static void invoke(std::ostream& Stream, int variantID);
};

}
}
}

#endif /*CRITTER__DECOMPOSITION__RECORD__RECORD_H_*/
