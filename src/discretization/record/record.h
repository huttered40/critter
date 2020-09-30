#ifndef CRITTER__DISCRETIZATION__RECORD__RECORD_H_
#define CRITTER__DISCRETIZATION__RECORD__RECORD_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

class record{
public:
  static void invoke(std::ofstream& Stream, int variantID, int print_mode, double overhead_time);
  static void invoke(std::ostream& Stream, int variantID, int print_mode, double overhead_time);

private:
  static std::vector<double> set_tuning_statistics(std::ofstream& Stream);
};

}
}
}

#endif /*CRITTER__DISCRETIZATION__RECORD__RECORD_H_*/
