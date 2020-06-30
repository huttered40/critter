#ifndef CRITTER__MECHANISM__PER_PROCESS__PER_PROCESS_H_
#define CRITTER__MECHANISM__PER_PROCESS__PER_PROCESS_H_

#include "../../util.h"

namespace critter{
namespace internal{

class per_process{
public:
  static void collect(MPI_Comm comm);
};

}
}

#endif /*CRITTER__MECHANISM__PER_PROCESS__PER_PROCESS_H_*/
