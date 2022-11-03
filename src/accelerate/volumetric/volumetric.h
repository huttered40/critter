#ifndef CRITTER__ACCELERATE__VOLUMETRIC__VOLUMETRIC_H_
#define CRITTER__ACCELERATE__VOLUMETRIC__VOLUMETRIC_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace accelerate{

class volumetric{
public:
  static void collect(MPI_Comm comm);
};

}
}
}

#endif /*CRITTER__ACCELERATE__VOLUMETRIC__VOLUMETRIC_H_*/
