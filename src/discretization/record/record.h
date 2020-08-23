#ifndef CRITTER__DISCRETIZATION__RECORD__RECORD_H_
#define CRITTER__DISCRETIZATION__RECORD__RECORD_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

class record{
public:
  static void invoke(std::ofstream& Stream, double* data, bool print_statistical_data, bool save_statistical_data);
  static void invoke(std::ostream& Stream, double* data, bool print_statistical_data, bool save_statistical_data);
};

}
}
}

#endif /*CRITTER__DISCRETIZATION__RECORD__RECORD_H_*/
