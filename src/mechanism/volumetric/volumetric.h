#ifndef CRITTER__MECHANISM__VOLUMETRIC__VOLUMETRIC_H_
#define CRITTER__MECHANISM__VOLUMETRIC__VOLUMETRIC_H_

#include "../../util.h"

namespace critter{
namespace internal{

class volumetric{
public:
  static void collect(MPI_Comm comm);
};

}
}

#endif /*CRITTER__MECHANISM__VOLUMETRIC__VOLUMETRIC_H_*/
