#ifndef CRITTER__DISCRETIZATION__VOLUMETRIC__VOLUMETRIC_H_
#define CRITTER__DISCRETIZATION__VOLUMETRIC__VOLUMETRIC_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

class volumetric{
public:
  static void collect(MPI_Comm comm);
};

}
}
}

#endif /*CRITTER__DISCRETIZATION__VOLUMETRIC__VOLUMETRIC_H_*/
