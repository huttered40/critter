#ifndef CRITTER__DECOMPOSITION__RECORD__RECORD_H_
#define CRITTER__DECOMPOSITION__RECORD__RECORD_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace decomposition{

class record{
public:
  static void invoke(std::ofstream& Stream);
  static void invoke(std::ostream& Stream);
};

}
}
}

#endif /*CRITTER__DECOMPOSITION__RECORD__RECORD_H_*/
