#ifndef CRITTER__SKELETONIZE__VOLUMETRIC__VOLUMETRIC_H_
#define CRITTER__SKELETONIZE__VOLUMETRIC__VOLUMETRIC_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace skeletonize{

class volumetric{
public:
  static void collect(MPI_Comm comm);
};

}
}
}

#endif /*CRITTER__SKELETONIZE__VOLUMETRIC__VOLUMETRIC_H_*/
