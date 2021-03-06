#ifndef CRITTER__SKELETONIZATION__VOLUMETRIC__VOLUMETRIC_H_
#define CRITTER__SKELETONIZATION__VOLUMETRIC__VOLUMETRIC_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace skeletonization{

class volumetric{
public:
  static void collect(MPI_Comm comm);
};

}
}
}

#endif /*CRITTER__SKELETONIZATION__VOLUMETRIC__VOLUMETRIC_H_*/
